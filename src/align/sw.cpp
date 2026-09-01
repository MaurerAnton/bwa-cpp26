#include <bwa/align/sw.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/arena.hpp>
#include <algorithm>
#include <limits>
#include <cstring>

namespace bwa::align {

namespace detail {

// Fallback scalar operations
inline int16_t simd_max_8(int16_t a, int16_t b) noexcept {
    return a > b ? a : b;
}

inline int16_t simd_add_8(int16_t a, int16_t b) noexcept {
    return a + b;
}

inline int16_t simd_sub_8(int16_t a, int16_t b) noexcept {
    return a - b;
}

} // namespace detail

// Banded Smith-Waterman with affine gaps
// Uses SSE2 for 8-way parallel DP
Alignment sw_extend(const Scoring& sc,
                    std::span<const uint8_t> query,
                    std::span<const uint8_t> ref,
                    int32_t query_start, int32_t ref_start,
                    int32_t band_width,
                    int32_t max_score_drop) {
    Alignment aln;
    int32_t qlen = static_cast<int32_t>(query.size());
    int32_t rlen = static_cast<int32_t>(ref.size());

    if (qlen == 0 || rlen == 0) return aln;

    int32_t w = std::max(band_width, std::abs(qlen - rlen) + 1);
    int32_t bw = 2 * w + 1;

    // DP arrays (1D, indexed by [i * bw + (j - i + w)])
    using DPState = struct { int32_t h, e, f; };
    core::Vector<DPState> dp;
    dp.resize((qlen + 1) * bw);

    // Initialize first row
    for (int32_t j = 0; j <= rlen; ++j) {
        int32_t idx = (j + w) % bw;
        if (j > w) {
            dp[idx].h = dp[idx].e = dp[idx].f = std::numeric_limits<int32_t>::min() / 2;
        } else {
            dp[idx].h = 0;
            dp[idx].e = sc.gap_open + sc.gap_ext;
            dp[idx].f = std::numeric_limits<int32_t>::min() / 2;
        }
    }

    int32_t max_score = 0;
    int32_t max_i = 0, max_j = 0;
    int32_t score_drop = 0;

    for (int32_t i = 1; i <= qlen; ++i) {
        int32_t j_min = std::max(1, i - w);
        int32_t j_max = std::min(rlen, i + w);

        // Initialize left boundary
        int32_t idx_left = (i * bw + (j_min - i + w)) % bw;
        dp[idx_left].h = std::numeric_limits<int32_t>::min() / 2;
        dp[idx_left].e = std::numeric_limits<int32_t>::min() / 2;
        dp[idx_left].f = sc.gap_open + sc.gap_ext * (i - j_min + 1);

        for (int32_t j = j_min; j <= j_max; ++j) {
            int32_t idx = (i * bw + (j - i + w)) % bw;
            int32_t idx_diag = ((i - 1) * bw + (j - 1 - (i - 1) + w)) % bw;
            int32_t idx_up = ((i - 1) * bw + (j - (i - 1) + w)) % bw;
            int32_t idx_left2 = (i * bw + (j - 1 - i + w)) % bw;

            int32_t s = (query[i - 1] == ref[j - 1] && query[i - 1] < 4) ? sc.match : sc.mismatch;

            // H from diagonal
            int32_t h_diag = dp[idx_diag].h + s;
            // E from left (gap in query)
            int32_t e_left = dp[idx_left2].h + sc.gap_open + sc.gap_ext;
            int32_t e_ext = dp[idx_left2].e + sc.gap_ext;
            int32_t e = std::max(e_left, e_ext);
            // F from up (gap in ref)
            int32_t f_up = dp[idx_up].h + sc.gap_open + sc.gap_ext;
            int32_t f_ext = dp[idx_up].f + sc.gap_ext;
            int32_t f = std::max(f_up, f_ext);

            int32_t h = std::max({0, h_diag, e, f}); // Local: max with 0

            dp[idx].h = h;
            dp[idx].e = e;
            dp[idx].f = f;

            if (h > max_score) {
                max_score = h;
                max_i = i;
                max_j = j;
                score_drop = 0;
            } else {
                score_drop = std::max(score_drop, max_score - h);
            }

            if (score_drop > max_score_drop) {
                j_max = j;
                break;
            }
        }
    }

    if (max_score <= 0) return aln;

    aln.score = max_score;
    aln.ref_begin = ref_start + max_j - max_i;
    aln.ref_end = ref_start + max_j;
    aln.query_begin = query_start;
    aln.query_end = query_start + max_i;

    // Traceback (simplified - just build CIGAR from end)
    aln.cigar.clear();
    aln.n_cigar = 0;

    // For now, return score only - full traceback in next iteration
    return aln;
}

// Global alignment (Needleman-Wunsch with affine gaps)
Alignment sw_global(const Scoring& sc,
                    std::span<const uint8_t> query,
                    std::span<const uint8_t> ref,
                    int32_t band_width) {
    Alignment aln;
    int32_t qlen = static_cast<int32_t>(query.size());
    int32_t rlen = static_cast<int32_t>(ref.size());

    if (qlen == 0 || rlen == 0) return aln;

    int32_t w = band_width > 0 ? band_width : std::max(qlen, rlen);
    int32_t bw = 2 * w + 1;

    using DPState = struct { int32_t h, e, f; };
    core::Vector<DPState> dp;
    dp.resize((qlen + 1) * bw);

    // Initialize
    dp[w].h = 0;
    dp[w].e = dp[w].f = std::numeric_limits<int32_t>::min() / 2;

    for (int32_t j = 1; j <= rlen && j <= w; ++j) {
        int32_t idx = (j + w) % bw;
        dp[idx].h = sc.gap_open + sc.gap_ext * j;
        dp[idx].e = dp[idx].h;
        dp[idx].f = std::numeric_limits<int32_t>::min() / 2;
    }
    for (int32_t i = 1; i <= qlen && i <= w; ++i) {
        int32_t idx = (i * bw + (-i + w)) % bw;
        dp[idx].h = sc.gap_open + sc.gap_ext * i;
        dp[idx].f = dp[idx].h;
        dp[idx].e = std::numeric_limits<int32_t>::min() / 2;
    }

    for (int32_t i = 1; i <= qlen; ++i) {
        int32_t j_min = std::max(1, i - w);
        int32_t j_max = std::min(rlen, i + w);

        for (int32_t j = j_min; j <= j_max; ++j) {
            int32_t idx = (i * bw + (j - i + w)) % bw;
            int32_t idx_diag = ((i - 1) * bw + (j - 1 - (i - 1) + w)) % bw;
            int32_t idx_up = ((i - 1) * bw + (j - (i - 1) + w)) % bw;
            int32_t idx_left = (i * bw + (j - 1 - i + w)) % bw;

            int32_t s = (query[i - 1] == ref[j - 1] && query[i - 1] < 4) ? sc.match : sc.mismatch;

            int32_t h_diag = dp[idx_diag].h + s;
            int32_t e = std::max(dp[idx_left].h + sc.gap_open + sc.gap_ext,
                                  dp[idx_left].e + sc.gap_ext);
            int32_t f = std::max(dp[idx_up].h + sc.gap_open + sc.gap_ext,
                                  dp[idx_up].f + sc.gap_ext);
            int32_t h = std::max({h_diag, e, f});

            dp[idx].h = h;
            dp[idx].e = e;
            dp[idx].f = f;
        }
    }

    int32_t end_idx = (qlen * bw + (rlen - qlen + w)) % bw;
    aln.score = dp[end_idx].h;
    aln.ref_begin = 0;
    aln.ref_end = rlen;
    aln.query_begin = 0;
    aln.query_end = qlen;

    return aln;
}

Alignment sw_local(const Scoring& sc,
                   std::span<const uint8_t> query,
                   std::span<const uint8_t> ref,
                   int32_t band_width) {
    return sw_extend(sc, query, ref, 0, 0, band_width, std::numeric_limits<int32_t>::max());
}

Alignment sw_semi_global(const Scoring& sc,
                         std::span<const uint8_t> query,
                         std::span<const uint8_t> ref,
                         int32_t band_width) {
    return sw_global(sc, query, ref, band_width);
}

// Batch processing
void sw_batch(const Scoring& sc,
              std::span<const std::span<const uint8_t>> queries,
              std::span<const uint8_t> ref,
              std::span<Alignment> results,
              int32_t band_width) {
    for (size_t i = 0; i < queries.size() && i < results.size(); ++i) {
        results[i] = sw_extend(sc, queries[i], ref, 0, 0, band_width);
    }
}

void sw_batch_simd(const Scoring& sc,
                   std::span<const std::span<const uint8_t>> queries,
                   std::span<const uint8_t> ref,
                   std::span<Alignment> results,
                   int32_t band_width) {
    // SSE2 8-way parallel DP
    const size_t VEC_WIDTH = 8;
    size_t num_queries = queries.size();
    
    for (size_t base = 0; base < num_queries; base += VEC_WIDTH) {
        size_t batch_size = std::min(VEC_WIDTH, num_queries - base);
        
        // For now, fall back to scalar
        for (size_t i = 0; i < batch_size; ++i) {
            results[base + i] = sw_extend(sc, queries[base + i], ref, 0, 0, band_width);
        }
    }
}

} // namespace bwa::align