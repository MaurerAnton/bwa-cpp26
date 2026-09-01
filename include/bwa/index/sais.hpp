#pragma once

#include <bwa/core/vector.hpp>
#include <bwa/core/arena.hpp>
#include <span>
#include <cstdint>
#include <cstddef>
#include <array>
#include <algorithm>

// Forward declaration to avoid circular dependency
namespace bwa::index {
    class PackedSequence;
}

namespace bwa::index::detail {

// SA-IS: Linear-time Suffix Array Construction
// Based on: Nong, Zhang, Chan (2009) "Two Efficient Algorithms for Linear Time Suffix Array Construction"
//
// Algorithm overview:
// 1. Classify each position as L-type or S-type using LMS suffix property
// 2. Find LMS substrings (LMS positions followed by S-types)
// 3. Sort LMS substrings using induced sorting
// 4. If LMS substrings are unique, derive SA; otherwise recurse on LMS string
//
// This is O(n) time and O(n) space.

namespace sais {

// Type classification
constexpr uint8_t L_TYPE = 0;
constexpr uint8_t S_TYPE = 1;
constexpr bool LMS_SUFFIX = true;
constexpr bool NOT_LMS_SUFFIX = false;

// Step 1: Classify T[i] as L-type or S-type
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

// Find LMS positions (S-type preceded by L-type)
inline void find_lms_positions(const uint8_t* t, int32_t n, core::Vector<int32_t>& lms) noexcept {
    for (int32_t i = 1; i < n; ++i) {
        if (t[i] == S_TYPE && t[i - 1] == L_TYPE) {
            lms.push_back(i);
        }
    }
}

// Compare two LMS substrings for equality
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

// Induced sort L-type suffixes
inline void induced_sort_l(const uint8_t* T, const uint8_t* t, int32_t n,
                            int32_t sigma, int32_t* SA,
                            const int32_t* bucket_ends,
                            memory::Arena& arena) noexcept {
    auto* processed = static_cast<bool*>(arena.allocate(n * sizeof(bool), alignof(bool)));
    for (int32_t i = 0; i < n; ++i) processed[i] = false;

    for (int32_t i = n - 1; i >= 0; --i) {
        if (SA[i] <= 0) continue;
        int32_t j = SA[i] - 1;
        if (t[j] == L_TYPE && !processed[j]) {
            SA[bucket_ends[T[j]] - 1] = j;
            processed[j] = true;
        }
    }
    arena.deallocate(processed, n * sizeof(bool), alignof(bool));
}

// Induced sort S-type suffixes
inline void induced_sort_s(const uint8_t* T, const uint8_t* t, int32_t n,
                            int32_t* SA,
                            const int32_t* bucket_starts) noexcept {
    for (int32_t i = 0; i < n; ++i) {
        if (SA[i] <= 0) continue;
        int32_t j = SA[i] - 1;
        if (t[j] == S_TYPE) {
            SA[bucket_starts[T[j]]] = j;
        }
    }
}

// Compute bucket boundaries
inline void compute_buckets(const uint8_t* T, int32_t n, int32_t sigma,
                             int32_t* bucket_starts, int32_t* bucket_ends) noexcept {
    int32_t counts[256] = {};
    for (int32_t i = 0; i < n; ++i) {
        ++counts[T[i]];
    }

    int32_t sum = 0;
    for (int32_t i = 0; i < sigma; ++i) {
        bucket_starts[i] = sum;
        sum += counts[i];
        bucket_ends[i] = sum;
    }
}

// Recursive SA-IS
inline int32_t sa_is_recursive(const uint8_t* T, int32_t n, int32_t* SA,
                                int32_t sigma, memory::Arena& arena) {
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
    core::Vector<int32_t> sorted_lms(&arena);
    sorted_lms.reserve(n / 2);

    uint8_t* T_t = t.data();
    int32_t* bs = bucket_starts.data();
    int32_t* be = bucket_ends.data();

    classify_types(T, n, T_t);
    compute_buckets(T, n, sigma, bs, be);

    for (int32_t i = 0; i < n; ++i) SA[i] = -1;
    find_lms_positions(T_t, n, lms_positions);
    for (int32_t i = (int32_t)lms_positions.size() - 1; i >= 0; --i) {
        int32_t pos = lms_positions[i];
        SA[be[T[pos]] - 1] = pos;
        --be[T[pos]];
    }

    induced_sort_l(T, T_t, n, sigma, SA, be, arena);
    induced_sort_s(T, T_t, n, SA, bs);

    int32_t num_lms = (int32_t)lms_positions.size();
    int32_t lms_count = 0;
    for (int32_t i = 0; i < n; ++i) {
        if (SA[i] > 0 && T_t[SA[i]] == S_TYPE && T_t[SA[i] - 1] == L_TYPE) {
            SA[lms_count++] = SA[i];
        }
    }
    for (int32_t i = lms_count; i < n; ++i) SA[i] = -1;

    int32_t name = 0;
    int32_t prev = -1;
    for (int32_t i = 0; i < lms_count; ++i) {
        int32_t pos = SA[i];
        if (prev == -1 || !lms_substrings_equal(T, T_t, pos, prev, n)) {
            ++name;
        }
        sorted_lms.push_back(name - 1);
        prev = pos;
    }

    if (name == lms_count) {
        // All LMS substrings are unique - place them in SA in name order
        // Build inverse mapping: name -> position in lms_positions
        core::Vector<int32_t> name_to_pos(&arena);
        name_to_pos.resize(name);
        for (int32_t i = 0; i < lms_count; ++i) {
            name_to_pos[sorted_lms[i]] = lms_positions[i];
        }
        for (int32_t i = lms_count - 1; i >= 0; --i) {
            int32_t pos = name_to_pos[i];
            SA[be[T[pos]] - 1] = pos;
            --be[T[pos]];
        }
    } else {
        core::Vector<uint8_t> T_reduced(&arena);
        T_reduced.resize(lms_count + 1);
        T_reduced[lms_count] = 0;
        for (int32_t i = 0; i < lms_count; ++i) {
            T_reduced[i] = (uint8_t)(sorted_lms[i] + 1);
        }

        core::Vector<int32_t> SA_reduced(&arena);
        SA_reduced.resize(lms_count + 1);
        sa_is_recursive(T_reduced.data(), lms_count + 1, SA_reduced.data(), name + 1, arena);

        for (int32_t i = lms_count - 1; i >= 0; --i) {
            int32_t idx = SA_reduced[i + 1];
            int32_t pos = lms_positions[idx - 1];
            SA[be[T[pos]] - 1] = pos;
            --be[T[pos]];
        }
    }

    induced_sort_l(T, T_t, n, sigma, SA, be, arena);
    induced_sort_s(T, T_t, n, SA, bs);

    return 0;
}

} // namespace sais

} // namespace bwa::index::detail