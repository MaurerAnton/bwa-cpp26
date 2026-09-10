#include <bwa/index/sais.hpp>
#include <bwa/index/fm_index.hpp>

#include <algorithm>
#include <cstdint>

namespace bwa::index::detail::sais {

using namespace bwa::core;
using namespace bwa::memory;

constexpr uint8_t L_TYPE = 0;
constexpr uint8_t S_TYPE = 1;

// Classify each suffix as L-type (greater than its successor) or S-type.
// The last suffix (sentinel, smallest) is S-type.
inline void classify_types(const int32_t* T, int32_t n, uint8_t* t) noexcept {
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

// Collect LMS positions (S-type preceded by L-type) in text order.
inline void find_lms_positions(const uint8_t* t, int32_t n,
                               core::Vector<int32_t>& lms) noexcept {
    for (int32_t i = 1; i < n; ++i) {
        if (t[i] == S_TYPE && t[i - 1] == L_TYPE) {
            lms.push_back(i);
        }
    }
}

// Compare two LMS substrings: equal chars up to and including the next LMS
// boundary on both sides. Used only to name substrings, never past n.
inline bool lms_substrings_equal(const int32_t* T, const uint8_t* t,
                                 int32_t i, int32_t j, int32_t n) noexcept {
    if (i == j) return true;
    if (i >= n - 1 || j >= n - 1) return i == j;

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

// Bucket boundaries: starts[c] = first index of char c, ends[c] = one past
// the last index. Recomputed fresh before every placement/induced sort so
// no caller depends on consumed bucket state.
inline void compute_buckets(const int32_t* T, int32_t n, int32_t sigma,
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

// Induce L-type suffixes: scan left to right, planting each L-predecessor
// at its bucket head. Placed items to the right are revisited by the scan,
// which is exactly how induction propagates.
inline void induced_sort_l(const int32_t* T, const uint8_t* t, int32_t n,
                           int32_t* SA, int32_t* bucket_heads) noexcept {
    for (int32_t i = 0; i < n; ++i) {
        int32_t s = SA[i];
        if (s <= 0) continue;  // empty cell, or suffix 0 (no predecessor)
        int32_t j = s - 1;
        if (t[j] == L_TYPE) {
            SA[bucket_heads[T[j]]++] = j;
        }
    }
}

// Induce S-type suffixes: scan right to left, planting each S-predecessor
// at its bucket tail.
inline void induced_sort_s(const int32_t* T, const uint8_t* t, int32_t n,
                           int32_t* SA, int32_t* bucket_tails) noexcept {
    for (int32_t i = n - 1; i >= 0; --i) {
        int32_t s = SA[i];
        if (s <= 0) continue;  // empty cell, or suffix 0 (no predecessor)
        int32_t j = s - 1;
        if (t[j] == S_TYPE) {
            SA[--bucket_tails[T[j]]] = j;
        }
    }
}

// Main SA-IS recursive function.
// T: input string, T[n-1] = 0 sentinel, T[i] >= 1 for i < n-1.
// n: length of T including sentinel. SA: output suffix array (size n).
// sigma: alphabet size (max char value + 1).
int32_t sa_is_recursive(const int32_t* T, int32_t n, int32_t* SA,
                        int32_t sigma, Arena& arena) noexcept {
    if (n <= 1) {
        if (n == 1) SA[0] = 0;
        return 0;
    }

    core::Vector<uint8_t> t(&arena);
    t.resize(static_cast<size_t>(n));
    core::Vector<int32_t> bucket_heads(&arena);
    bucket_heads.resize(static_cast<size_t>(sigma));
    core::Vector<int32_t> bucket_tails(&arena);
    bucket_tails.resize(static_cast<size_t>(sigma));
    core::Vector<int32_t> lms_positions(&arena);
    lms_positions.reserve(static_cast<size_t>(n) / 2 + 1);

    classify_types(T, n, t.data());
    find_lms_positions(t.data(), n, lms_positions);
    const int32_t m = static_cast<int32_t>(lms_positions.size());
    // The sentinel forces t[n-2] == L (T[n-2] > T[n-1] == 0), so position
    // n-1 is always LMS: every input with n >= 2 has at least one.
    if (m == 0) return -1;  // unreachable; guards against corrupt output

    for (int32_t i = 0; i < n; ++i) SA[i] = -1;

    // Seed bucket tails with LMS suffixes in text order
    compute_buckets(T, n, sigma, bucket_heads.data(), bucket_tails.data());
    for (int32_t i = 0; i < m; ++i) {
        int32_t pos = lms_positions[i];
        SA[--bucket_tails[T[pos]]] = pos;
    }

    // Induce L then S (fresh buckets each; see compute_buckets docs)
    compute_buckets(T, n, sigma, bucket_heads.data(), bucket_tails.data());
    induced_sort_l(T, t.data(), n, SA, bucket_heads.data());
    compute_buckets(T, n, sigma, bucket_heads.data(), bucket_tails.data());
    induced_sort_s(T, t.data(), n, SA, bucket_tails.data());

    // Compact the now-sorted LMS suffixes to SA[0..m)
    int32_t compacted = 0;
    for (int32_t i = 0; i < n; ++i) {
        int32_t s = SA[i];
        if (s > 0 && t[s] == S_TYPE && t[s - 1] == L_TYPE) {
            SA[compacted++] = s;
        }
    }
    for (int32_t i = compacted; i < n; ++i) SA[i] = -1;
    if (compacted != m) return -1;  // internal error guard

    // Name LMS substrings in sorted order
    core::Vector<int32_t> lms_names(&arena);
    lms_names.resize(static_cast<size_t>(m));
    int32_t name = 0;
    int32_t prev = -1;
    for (int32_t i = 0; i < m; ++i) {
        int32_t pos = SA[i];
        if (prev == -1 || !lms_substrings_equal(T, t.data(), pos, prev, n)) {
            ++name;
        }
        lms_names[i] = name - 1;
        prev = pos;
    }

    // Copy sorted LMS aside: the placement below overwrites SA in place and
    // must not clobber not-yet-read entries.
    core::Vector<int32_t> sorted_lms(&arena);
    sorted_lms.resize(static_cast<size_t>(m));
    for (int32_t i = 0; i < m; ++i) sorted_lms[i] = SA[i];

    // Clear SA: the compacted prefix [0, m) holds stale LMS values that the
    // tail placement below would otherwise leave behind as duplicates.
    for (int32_t i = 0; i < n; ++i) SA[i] = -1;

    compute_buckets(T, n, sigma, bucket_heads.data(), bucket_tails.data());
    if (name == m) {
        // All substrings distinct: compacted order is the sorted LMS order.
        // Seed tails from largest to smallest.
        for (int32_t i = m - 1; i >= 0; --i) {
            int32_t pos = sorted_lms[i];
            SA[--bucket_tails[T[pos]]] = pos;
        }
    } else {
        // Reduced string in TEXT order of LMS positions. lms_positions is
        // text-ordered, so binary-search each sorted position's text rank.
        // (O(m log m), no extra O(n) array.)
        core::Vector<int32_t> reduced(&arena);
        reduced.resize(static_cast<size_t>(m) + 1);
        reduced[m] = 0;  // sentinel
        for (int32_t i = 0; i < m; ++i) {
            const int32_t* begin = lms_positions.begin();
            const int32_t* end = lms_positions.end();
            int32_t text_idx =
                static_cast<int32_t>(std::lower_bound(begin, end, sorted_lms[i]) - begin);
            reduced[text_idx] = lms_names[i] + 1;
        }

        core::Vector<int32_t> reduced_sa(&arena);
        reduced_sa.resize(static_cast<size_t>(m) + 1);
        if (sa_is_recursive(reduced.data(), m + 1, reduced_sa.data(), name + 1, arena) != 0) {
            return -1;
        }

        // Seed tails with LMS suffixes in reduced-suffix order
        // (reduced_sa[0] is the reduced sentinel; skip it)
        for (int32_t i = m - 1; i >= 0; --i) {
            int32_t text_idx = reduced_sa[i + 1];
            int32_t pos = lms_positions[text_idx];
            SA[--bucket_tails[T[pos]]] = pos;
        }
    }

    // Final induced sorting (fresh buckets each)
    compute_buckets(T, n, sigma, bucket_heads.data(), bucket_tails.data());
    induced_sort_l(T, t.data(), n, SA, bucket_heads.data());
    compute_buckets(T, n, sigma, bucket_heads.data(), bucket_tails.data());
    induced_sort_s(T, t.data(), n, SA, bucket_tails.data());

    return 0;
}

// Wrapper for PackedSequence. Contract (see sais.hpp): returns all n+1
// suffixes INCLUDING the empty suffix (value n) sorted first, with codes
// sentinel=0 < N=1 < A=2 < C=3 < G=4 < T=5.
#if defined(__has_feature)
#  if __has_feature(address_sanitizer)
#    define SAIS_NO_SANITIZE __attribute__((no_sanitize("address")))
#  else
#    define SAIS_NO_SANITIZE
#  endif
#else
#  define SAIS_NO_SANITIZE
#endif

SAIS_NO_SANITIZE
core::Vector<uint32_t> build_suffix_array(const PackedSequence& seq,
                                          Arena& arena) noexcept {
    size_t n = seq.size();
    if (n + 1 > static_cast<size_t>(INT32_MAX)) return {};  // caller falls back
    if (n == 0) {
        core::Vector<uint32_t> only(&arena);
        only.resize(1);
        only[0] = 0;
        return only;
    }

    core::Vector<int32_t> T(&arena);
    T.resize(n + 1);
    T[n] = 0;
    for (size_t i = 0; i < n; ++i) {
        uint8_t b = seq.get(i);
        T[i] = (b == 4) ? 1 : (b + 2);
    }

    core::Vector<int32_t> SA(&arena);
    SA.resize(n + 1);
    if (sa_is_recursive(T.data(), static_cast<int32_t>(n + 1), SA.data(), 6, arena) != 0) {
        return {};
    }

    core::Vector<uint32_t> result(&arena);
    result.resize(n + 1);
    for (size_t i = 0; i <= n; ++i) {
        result[i] = static_cast<uint32_t>(SA[i]);
    }
    return result;
}

} // namespace bwa::index::detail::sais
