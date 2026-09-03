#pragma once

#include <bwa/core/vector.hpp>
#include <bwa/core/arena.hpp>
#include <span>
#include <cstdint>
#include <cstddef>
#include <array>
#include <algorithm>

namespace bwa::index::detail {

// SA-IS: Linear-time Suffix Array Construction
// Based on: Nong, Zhang, Chan (2009) "Two Efficient Algorithms for Linear Time Suffix Array Construction"
// Implementation reference: https://github.com/algorithm-ninja/sais

namespace sais {

constexpr uint8_t L_TYPE = 0;
constexpr uint8_t S_TYPE = 1;

inline void classify_types(const uint8_t* T, int32_t n, uint8_t* t) noexcept {
    if (n == 0) return;
    t[n - 1] = S_TYPE;
    for (int32_t i = n - 2; i >= 0; --i) {
        if (T[i] < T[i + 1]) {
            t[i] = S_TYPE;
        } else if (T[i] > T[i + 1]) {
            t[i] = L_TYPE;
        } else {
            t[i] = t[i + 1];
        }
    }
}

inline void find_lms_positions(const uint8_t* t, int32_t n, core::Vector<int32_t>& lms) noexcept {
    for (int32_t i = 1; i < n; ++i) {
        if (t[i] == S_TYPE && t[i - 1] == L_TYPE) {
            lms.push_back(i);
        }
    }
}

inline bool lms_substrings_equal(const uint8_t* T, const uint8_t* t,
                                  int32_t i, int32_t j, int32_t n) noexcept {
    if (i == j) return true;
    if (i == n - 1 || j == n - 1) return false;

    int32_t k = 0;
    while (true) {
        bool i_lms = (k > 0) && (t[i + k] == S_TYPE) && (t[i + k - 1] == L_TYPE);
        bool j_lms = (k > 0) && (t[j + k] == S_TYPE) && (t[j + k - 1] == L_TYPE);
        if (k > 0 && i_lms && j_lms) return true;
        if (i_lms != j_lms) return false;
        if (T[i + k] != T[j + k]) return false;
        ++k;
        if (i + k >= n || j + k >= n) return i + k == n && j + k == n;
    }
}

inline void compute_buckets(const uint8_t* T, int32_t n, int32_t sigma,
                             int32_t* bucket_starts, int32_t* bucket_ends) noexcept {
    for (int32_t i = 0; i < sigma; ++i) {
        bucket_starts[i] = 0;
        bucket_ends[i] = 0;
    }
    for (int32_t i = 0; i < n; ++i) {
        ++bucket_ends[T[i]];
    }
    int32_t sum = 0;
    for (int32_t i = 0; i < sigma; ++i) {
        bucket_starts[i] = sum;
        sum += bucket_ends[i];
        bucket_ends[i] = sum;
    }
}

inline void induced_sort_l(const uint8_t* T, const uint8_t* t, int32_t n,
                            int32_t sigma, int32_t* SA,
                            int32_t* bucket_ends,
                            memory::Arena& arena) noexcept {
    auto* processed = static_cast<bool*>(arena.allocate(n * sizeof(bool), alignof(bool)));
    for (int32_t i = 0; i < n; ++i) processed[i] = false;

    for (int32_t i = n - 1; i >= 0; --i) {
        if (SA[i] <= 0) continue;
        int32_t j = SA[i] - 1;
        if (t[j] == L_TYPE && !processed[j]) {
            SA[--bucket_ends[T[j]]] = j;
            processed[j] = true;
        }
    }
    arena.deallocate(processed, n * sizeof(bool), alignof(bool));
}

inline void induced_sort_s(const uint8_t* T, const uint8_t* t, int32_t n,
                            int32_t* SA,
                            const int32_t* bucket_starts) noexcept {
    for (int32_t i = 0; i < n; ++i) {
        if (SA[i] <= 0) continue;
        int32_t j = SA[i] - 1;
        if (t[j] == S_TYPE) {
            SA[bucket_starts[T[j]]++] = j;
        }
    }
}

// Main SA-IS recursive function
// T: input string (with sentinel T[n-1] = 0, T[i] >= 1 for i < n-1)
// n: length of T including sentinel
// SA: output suffix array (size n)
// sigma: alphabet size (max char value + 1)
// Returns: 0 on success
inline int32_t sa_is_recursive(const uint8_t* T, int32_t n, int32_t* SA,
                                int32_t sigma, memory::Arena& arena) noexcept {
    if (n <= 1) {
        if (n == 1) SA[0] = 0;
        return 0;
    }

    core::Vector<uint8_t> t(&arena);
    t.resize(n);
    core::Vector<int32_t> bucket_starts(&arena);
    bucket_starts.resize(sigma);
    core::Vector<int32_t> bucket_ends(&arena);
    bucket_ends.resize(sigma);
    core::Vector<int32_t> lms_positions(&arena);
    lms_positions.reserve(n / 2);

    classify_types(T, n, t.data());
    compute_buckets(T, n, sigma, bucket_starts.data(), bucket_ends.data());

    for (int32_t i = 0; i < n; ++i) SA[i] = -1;
    find_lms_positions(t.data(), n, lms_positions);

    // Place LMS suffixes into SA
    int32_t lms_count = (int32_t)lms_positions.size();
    for (int32_t i = lms_count - 1; i >= 0; --i) {
        int32_t pos = lms_positions[i];
        SA[--bucket_ends[T[pos]]] = pos;
    }

    // Induced sorting
    induced_sort_l(T, t.data(), n, sigma, SA, bucket_ends.data(), arena);
    induced_sort_s(T, t.data(), n, SA, bucket_starts.data());

    // Compact LMS suffixes
    int32_t compacted = 0;
    for (int32_t i = 0; i < n; ++i) {
        int32_t pos = SA[i];
        if (pos > 0 && t[pos] == S_TYPE && t[pos - 1] == L_TYPE) {
            SA[compacted++] = pos;
        }
    }
    for (int32_t i = compacted; i < n; ++i) SA[i] = -1;

    // Name LMS substrings
    int32_t name = 0;
    int32_t prev = -1;
    core::Vector<int32_t> lms_names(&arena);
    lms_names.resize(compacted);
    
    for (int32_t i = 0; i < compacted; ++i) {
        int32_t pos = SA[i];
        if (prev == -1 || !lms_substrings_equal(T, t.data(), pos, prev, n)) {
            ++name;
        }
        lms_names[i] = name - 1;
        prev = pos;
    }

    if (name == compacted) {
        // All LMS substrings are unique
        // name_to_pos maps name -> LMS position
        core::Vector<int32_t> name_to_pos(&arena);
        name_to_pos.resize(name);
        for (int32_t i = 0; i < compacted; ++i) {
            name_to_pos[lms_names[i]] = lms_positions[i];
        }
        // Recompute bucket ends
        compute_buckets(T, n, sigma, bucket_starts.data(), bucket_ends.data());
        for (int32_t i = name - 1; i >= 0; --i) {
            int32_t pos = name_to_pos[i];
            SA[--bucket_ends[T[pos]]] = pos;
        }
    } else {
        // Recurse on reduced string
        core::Vector<uint8_t> T_reduced(&arena);
        T_reduced.resize(compacted + 1);
        T_reduced[compacted] = 0;  // sentinel
        for (int32_t i = 0; i < compacted; ++i) {
            T_reduced[i] = static_cast<uint8_t>(lms_names[i] + 1);
        }

        core::Vector<int32_t> SA_reduced(&arena);
        SA_reduced.resize(compacted + 1);
        sa_is_recursive(T_reduced.data(), compacted + 1, SA_reduced.data(), name + 1, arena);

        // Recompute bucket ends
        compute_buckets(T, n, sigma, bucket_starts.data(), bucket_ends.data());

        // Place LMS suffixes according to SA_reduced order
        for (int32_t i = compacted - 1; i >= 0; --i) {
            int32_t idx = SA_reduced[i + 1];  // SA_reduced[0] is sentinel
            int32_t pos = lms_positions[idx];
            SA[--bucket_ends[T[pos]]] = pos;
        }
    }

    // Final induced sorting
    induced_sort_l(T, t.data(), n, sigma, SA, bucket_ends.data(), arena);
    induced_sort_s(T, t.data(), n, SA, bucket_starts.data());

    return 0;
}

// Wrapper for PackedSequence
inline core::Vector<uint32_t> build_suffix_array(const PackedSequence& seq,
                                                  memory::Arena& arena) noexcept {
    size_t n = seq.size();
    if (n == 0) return {};

    // Convert to byte array with sentinel
    // Sentinel = 0, bases = 1..5 (A=1, C=2, G=3, T=4, N=5)
    core::Vector<uint8_t> T(&arena);
    T.resize(n + 1);
    T[n] = 0;
    for (size_t i = 0; i < n; ++i) {
        uint8_t b = seq.get(i);
        if (b == 4) T[i] = 5;  // N -> 5
        else T[i] = b + 1;     // A=0->1, C=1->2, G=2->3, T=3->4
    }

    core::Vector<int32_t> SA(&arena);
    SA.resize(n + 1);
    sa_is_recursive(T.data(), static_cast<int32_t>(n + 1), SA.data(), 6, arena);

    // Remove sentinel position from SA
    core::Vector<uint32_t> result(&arena);
    result.resize(n);
    int32_t idx = 0;
    for (int32_t i = 0; i <= (int32_t)n; ++i) {
        if (SA[i] != (int32_t)n) {
            result[idx++] = static_cast<uint32_t>(SA[i]);
        }
    }
    return result;
}

} // namespace sais

} // namespace bwa::index::detail