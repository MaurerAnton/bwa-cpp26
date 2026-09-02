#include <bwa/align/sw.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/arena.hpp>
#include <algorithm>
#include <limits>
#include <cstring>
#include <vector>

namespace bwa::align {

namespace detail {

// Traceback direction enum
enum class TraceDir : uint8_t {
    None = 0,
    Match = 1,   // From diagonal (M)
    Insert = 2,  // From left (I - gap in query)
    Delete = 3   // From up (D - gap in ref)
};

// CIGAR operations
enum class CigarOp : uint8_t {
    Match = 0,   // M (alignment match)
    Ins = 1,     // I (insertion to reference)
    Del = 2,     // D (deletion from reference)
    Skip = 3,    // N (skipped region from reference)
    SoftClip = 4, // S (soft clipping)
    HardClip = 5, // H (hard clipping)
    Pad = 6,     // P (padding)
    Equal = 7,   // = (sequence match)
    Diff = 8     // X (sequence mismatch)
};

inline constexpr uint32_t encode_cigar(int len, CigarOp op) noexcept {
    return (static_cast<uint32_t>(len) << 4) | static_cast<uint32_t>(op);
}

} // namespace detail

// Banded Smith-Waterman with affine gaps and full traceback
Alignment sw_extend(const Scoring& sc,
                    std::span<const uint8_t> query,
                    std::span<const uint8_t> ref,
                    int32_t query_start, int32_t ref_start,
                    int32_t band_width,
                    int32_t max_score_drop) {
    int32_t qlen = static_cast<int32_t>(query.size());
    int32_t rlen = static_cast<int32_t>(ref.size());

    if (qlen == 0 || rlen == 0) return Alignment{};

    int32_t w = std::max(band_width, std::abs(qlen - rlen) + 1);
    int32_t bw = 2 * w + 1;

    // DP arrays (full 2D for traceback)
    using DPState = struct { int32_t h, e, f; };
    std::vector<DPState> dp((qlen + 1) * bw);
    std::vector<uint8_t> trace((qlen + 1) * bw); // 0=None, 1=Match, 2=Insert, 3=Delete

    auto idx = [&](int32_t i, int32_t j) -> int32_t {
        return i * bw + (j - i + w);
    };

    // Initialize first row (i=0)
    for (int32_t j = 0; j <= rlen; ++j) {
        int32_t j_rel = j - 0 + w;
        if (j_rel < 0 || j_rel >= bw) continue;
        int32_t id = idx(0, j);
        if (j > w) {
            dp[id].h = dp[id].e = dp[id].f = std::numeric_limits<int32_t>::min() / 2;
        } else {
            dp[id].h = 0;
            dp[id].e = sc.gap_open + sc.gap_ext;
            dp[id].f = std::numeric_limits<int32_t>::min() / 2;
        }
    }

    int32_t max_score = 0;
    int32_t max_i = 0, max_j = 0;
    int32_t score_drop = 0;

    for (int32_t i = 1; i <= qlen; ++i) {
        int32_t j_min = std::max(1, i - w);
        int32_t j_max = std::min(rlen, i + w);

        // Initialize left boundary
        if (j_min <= j_max) {
            int32_t id_left = idx(i, j_min);
            dp[id_left].h = std::numeric_limits<int32_t>::min() / 2;
            dp[id_left].e = std::numeric_limits<int32_t>::min() / 2;
            dp[id_left].f = sc.gap_open + sc.gap_ext * (i - j_min + 1);
        }

        for (int32_t j = j_min; j <= j_max; ++j) {
            int32_t id = idx(i, j);
            int32_t id_diag = idx(i - 1, j - 1);
            int32_t id_up = idx(i - 1, j);
            int32_t id_left = idx(i, j - 1);

            int32_t s = (query[i - 1] == ref[j - 1] && query[i - 1] < 4) ? sc.match : sc.mismatch;

            // H from diagonal (match/mismatch)
            int32_t h_diag = dp[idx(i - 1, j - 1)].h + s;
            // E from left (gap in query / insertion)
            int32_t e_left = dp[idx(i, j - 1)].h + sc.gap_open + sc.gap_ext;
            int32_t e_ext = dp[idx(i, j - 1)].e + sc.gap_ext;
            int32_t e = std::max(e_left, e_ext);
            // F from up (gap in ref / deletion)
            int32_t f_up = dp[idx(i - 1, j)].h + sc.gap_open + sc.gap_ext;
            int32_t f_ext = dp[idx(i - 1, j)].f + sc.gap_ext;
            int32_t f = std::max(f_up, f_ext);

            int32_t h = std::max({0, h_diag, e, f}); // Local: max with 0

            // Determine traceback direction
            if (h > 0) {
                if (h == h_diag) {
                    trace[idx(i, j)] = 1; // Match
                } else if (h == e) {
                    trace[idx(i, j)] = 2; // Insert
                } else {
                    trace[idx(i, j)] = 3; // Delete
                }
            } else {
                trace[idx(i, j)] = 0; // None
            }

            dp[idx(i, j)].h = h;
            dp[idx(i, j)].e = e;
            dp[idx(i, j)].f = f;

            if (h > max_score) {
                max_score = h;
                max_i = i;
                max_j = j;
                score_drop = 0;
            } else {
                score_drop = std::max(score_drop, max_score - h);
            }

            if (score_drop > max_score_drop) {
                break;
            }
        }
    }

    if (max_score <= 0) return Alignment{};

    // Traceback to build CIGAR
    Alignment aln;
    aln.score = max_score;
    aln.ref_begin = ref_start + max_j - max_i;
    aln.ref_end = ref_start + max_j;
    aln.query_begin = query_start;
    aln.query_end = query_start + max_i;

    // Traceback from max position
    int32_t i = max_i, j = max_j;
    uint8_t prev_dir = 0; // 0=None
    int32_t cur_len = 0;

    auto idx_func = [&](int32_t i, int32_t j) -> int32_t {
        return i * bw + (j - i + w);
    };

    while (i > 0 || j > 0) {
        uint8_t dir = trace[idx(i, j)];

        if (dir == 0 || (i == 0 && j == 0)) break;

        if (dir == prev_dir && cur_len < 0xFFF) {
            ++cur_len;
        } else {
            if (cur_len > 0) {
                aln.cigar.push_back(detail::encode_cigar(cur_len, static_cast<detail::CigarOp>(prev_dir)));
                aln.n_cigar++;
            }
            prev_dir = dir;
            cur_len = 1;
        }

        // Move according to traceback direction
        if (dir == 1) { // Match
            --i; --j;
        } else if (dir == 2) { // Insert
            --j;
        } else if (dir == 3) { // Delete
            --i;
        } else {
            break;
        }
    }

    if (cur_len > 0) {
        aln.cigar.push_back(detail::encode_cigar(cur_len, static_cast<detail::CigarOp>(prev_dir)));
        aln.n_cigar++;
    }

    // Reverse CIGAR to get correct order
    std::reverse(aln.cigar.begin(), aln.cigar.end());

    // Count gaps
    aln.n_mismatch = aln.n_gap_open = aln.n_gap_ext = 0;
    for (uint32_t c : aln.cigar) {
        auto op = static_cast<detail::CigarOp>(c & 0xF);
        int len = c >> 4;
        if (op == detail::CigarOp::Ins) {
            aln.n_gap_open++;
            aln.n_gap_ext += len - 1;
        } else if (op == detail::CigarOp::Del) {
            aln.n_gap_open++;
            aln.n_gap_ext += len - 1;
        }
    }

    // Score and positions
    aln.score = max_score;
    aln.ref_begin = ref_start + max_j - max_i;
    aln.ref_end = ref_start + max_j;
    aln.query_begin = query_start;
    aln.query_end = query_start + max_i;

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
    std::vector<DPState> dp((qlen + 1) * bw);

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
    sw_batch(sc, queries, ref, results, band_width);
}

} // namespace bwa::align