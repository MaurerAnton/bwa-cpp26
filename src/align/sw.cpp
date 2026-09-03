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
// Semi-global: query aligned end-to-end, reference can have overhangs
Alignment sw_semi_global_extend(const Scoring& sc,
                                std::span<const uint8_t> query,
                                std::span<const uint8_t> ref,
                                int32_t band_width) {
    int32_t qlen = static_cast<int32_t>(query.size());
    int32_t rlen = static_cast<int32_t>(ref.size());

    if (qlen == 0 || rlen == 0) return Alignment{};

    int32_t w = std::max(band_width, std::abs(qlen - rlen) + 1);
    int32_t bw = 2 * w + 1;

    using DPState = struct { int32_t h, e, f; };

    // Use thread-local static buffers to avoid repeated allocations
    static thread_local std::vector<DPState> dp;
    static thread_local std::vector<uint8_t> trace;

    size_t dp_size = static_cast<size_t>(qlen + 1) * bw;
    if (dp.size() < dp_size) {
        dp.resize(dp_size);
        trace.resize(dp_size);
    } else {
        // Only reset the parts we use
        std::fill_n(dp.data(), dp_size, DPState{0, 0, 0});
        std::fill_n(trace.data(), dp_size, uint8_t(0));
    }

    // Precompute index function parameters
    const int32_t w_val = w;
    const int32_t bw_val = bw;
    auto idx = [w_val, bw_val](int32_t i, int32_t j) -> int32_t {
        return i * bw_val + (j - i + w_val);
    };

    // Initialize first row: semi-global allows free gaps in reference at start
    for (int32_t j = 0; j <= rlen && j <= w; ++j) {
        int32_t id = idx(0, j);
        dp[id].h = 0; // Free gaps in reference
        dp[id].e = 0;
        dp[id].f = std::numeric_limits<int32_t>::min() / 2;
    }

    // Initialize first column: must align query from start (penalize gaps)
    for (int32_t i = 1; i <= qlen && i <= w; ++i) {
        int32_t id = idx(i, 0);
        dp[id].h = sc.gap_open + sc.gap_ext * i;
        dp[id].f = dp[id].h;
        dp[id].e = std::numeric_limits<int32_t>::min() / 2;
        trace[id] = 3; // Delete
    }

    // Fill DP matrix - use direct index arithmetic for speed
    for (int32_t i = 1; i <= qlen; ++i) {
        int32_t j_min = std::max(1, i - w);
        int32_t j_max = std::min(rlen, i + w);

        int32_t base_idx = i * bw;
        int32_t prev_base_idx = (i - 1) * bw;

        for (int32_t j = j_min; j <= j_max; ++j) {
            int32_t cur_idx = base_idx + (j - i + w);

            int32_t s = (query[i - 1] == ref[j - 1] && query[i - 1] < 4) ? sc.match : sc.mismatch;

            int32_t h_diag = dp[prev_base_idx + (j - i + w - 1)].h + s;
            int32_t e_left = dp[base_idx + (j - 1 - i + w)].h + sc.gap_open + sc.gap_ext;
            int32_t e_ext = dp[base_idx + (j - 1 - i + w)].e + sc.gap_ext;
            int32_t e = std::max(e_left, e_ext);
            int32_t f_up = dp[prev_base_idx + (j - i + w)].h + sc.gap_open + sc.gap_ext;
            int32_t f_ext = dp[prev_base_idx + (j - i + w)].f + sc.gap_ext;
            int32_t f = std::max(f_up, f_ext);

            int32_t h = std::max({h_diag, e, f});

            dp[cur_idx].h = h;
            dp[cur_idx].e = e;
            dp[cur_idx].f = f;

            if (h == h_diag) {
                trace[cur_idx] = 1;
            } else if (h == e) {
                trace[cur_idx] = 2;
            } else {
                trace[cur_idx] = 3;
            }
        }
    }

    // Find best end position: must end at i=qlen (query fully aligned)
    // Allow any j in the last row
    int32_t best_score = std::numeric_limits<int32_t>::min() / 2;
    int32_t best_j = qlen;
    for (int32_t j = std::max(1, qlen - w); j <= std::min(rlen, qlen + w); ++j) {
        int32_t id = idx(qlen, j);
        if (dp[id].h > best_score) {
            best_score = dp[id].h;
            best_j = j;
        }
    }

    if (best_score <= 0) return Alignment{};

    // Traceback
    Alignment aln;
    int32_t ti = qlen, tj = best_j;
    int32_t cur_len = 0;
    detail::CigarOp prev_cigar_op = detail::CigarOp::Match;

    while (ti > 0 || tj > 0) {
        if (ti == 0) {
            // Remaining reference is unaligned (soft clip on query)
            break;
        }

        int32_t cur_idx = idx(ti, tj);
        uint8_t dir = trace[cur_idx];
        if (dir == 0) {
            if (ti == 0) break;
            // No trace - query has unaligned prefix
            break;
        }

        detail::CigarOp this_cigar_op;
        if (dir == 1) {
            bool is_match = (query[ti - 1] == ref[tj - 1] && query[ti - 1] < 4);
            this_cigar_op = is_match ? detail::CigarOp::Equal : detail::CigarOp::Diff;
        } else if (dir == 2) {
            this_cigar_op = detail::CigarOp::Ins;
        } else {
            this_cigar_op = detail::CigarOp::Del;
        }

        if (this_cigar_op == prev_cigar_op && cur_len < 0xFFF) {
            ++cur_len;
        } else {
            if (cur_len > 0) {
                aln.cigar.push_back(detail::encode_cigar(cur_len, prev_cigar_op));
            }
            prev_cigar_op = this_cigar_op;
            cur_len = 1;
        }

        if (dir == 1) { --ti; --tj; }
        else if (dir == 2) { --tj; }
        else { --ti; }
    }

    if (cur_len > 0) {
        aln.cigar.push_back(detail::encode_cigar(cur_len, prev_cigar_op));
    }

    // Add soft clip for unaligned query prefix
    if (ti > 0) {
        aln.cigar.push_back(detail::encode_cigar(ti, detail::CigarOp::SoftClip));
    }

    // Reverse CIGAR
    std::reverse(aln.cigar.begin(), aln.cigar.end());

    aln.score = best_score;
    aln.query_begin = 0;
    aln.query_end = qlen;
    aln.ref_begin = best_j - qlen;
    aln.ref_end = best_j;

    // Count mismatches
    aln.n_mismatch = aln.n_gap_open = aln.n_gap_ext = 0;
    for (uint32_t c : aln.cigar) {
        auto op = static_cast<detail::CigarOp>(c & 0xF);
        int len = c >> 4;
        if (op == detail::CigarOp::Ins || op == detail::CigarOp::Del) {
            aln.n_gap_open++;
            aln.n_gap_ext += len - 1;
        } else if (op == detail::CigarOp::Diff) {
            aln.n_mismatch += len;
        }
    }

    return aln;
}

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

    // Use thread-local static buffers to avoid repeated allocations
    static thread_local std::vector<DPState> dp;
    static thread_local std::vector<uint8_t> trace;

    size_t dp_size = static_cast<size_t>(qlen + 1) * bw;
    if (dp.size() < dp_size) {
        dp.resize(dp_size);
        trace.resize(dp_size);
    } else {
        // Only reset the parts we use
        std::fill_n(dp.data(), dp_size, DPState{0, 0, 0});
        std::fill_n(trace.data(), dp_size, uint8_t(0));
    }

    // Precompute index function parameters
    const int32_t w_val = w;
    const int32_t bw_val = bw;
    auto idx = [w_val, bw_val](int32_t i, int32_t j) -> int32_t {
        return i * bw_val + (j - i + w_val);
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

        int32_t base_idx = i * bw;
        int32_t prev_base_idx = (i - 1) * bw;

        for (int32_t j = j_min; j <= j_max; ++j) {
            int32_t cur_idx = base_idx + (j - i + w);

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
                // Map TraceDir to CigarOp: Match->Match, Insert->Ins, Delete->Del
                detail::CigarOp cigar_op;
                switch (prev_dir) {
                    case 1: cigar_op = detail::CigarOp::Match; break;
                    case 2: cigar_op = detail::CigarOp::Ins; break;
                    case 3: cigar_op = detail::CigarOp::Del; break;
                    default: cigar_op = detail::CigarOp::Match;
                }
                aln.cigar.push_back(detail::encode_cigar(cur_len, cigar_op));
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
        detail::CigarOp cigar_op;
        switch (prev_dir) {
            case 1: cigar_op = detail::CigarOp::Match; break;
            case 2: cigar_op = detail::CigarOp::Ins; break;
            case 3: cigar_op = detail::CigarOp::Del; break;
            default: cigar_op = detail::CigarOp::Match;
        }
        aln.cigar.push_back(detail::encode_cigar(cur_len, cigar_op));
        aln.n_cigar++;
    }

    // Reverse CIGAR to get correct order
    std::reverse(aln.cigar.begin(), aln.cigar.end());

    // Second pass: refine Match ops to = (match) or X (mismatch)
    // We need to walk the alignment again to check actual bases
    // Since CIGAR is in forward order (after reversal), but we track from end, iterate in reverse
    int32_t q_pos = max_i - 1; // Current query position (0-based, at end of alignment)
    int32_t r_pos = max_j - 1; // Current ref position (0-based, at end of alignment)
    
    // We'll build a new CIGAR with = and X split (in reverse order, then reverse at end)
    std::vector<uint32_t> new_cigar_rev;
    int32_t new_n_cigar = 0;
    
    // Iterate CIGAR in REVERSE order (from end of alignment to start)
    for (auto it = aln.cigar.rbegin(); it != aln.cigar.rend(); ++it) {
        uint32_t c = *it;
        auto op = static_cast<detail::CigarOp>(c & 0xF);
        int len = c >> 4;
        if (op == detail::CigarOp::Match) {
            // Split this Match run into = and X sub-runs based on actual bases
            int sub_len = 0;
            detail::CigarOp sub_op = detail::CigarOp::Match;
            bool sub_op_set = false;
            
            for (int k = 0; k < len; ++k) {
                // Since we're iterating CIGAR in reverse, we check bases from current position backwards
                bool is_match = (query[query_start + q_pos] == ref[ref_start + r_pos]);
                detail::CigarOp this_op = is_match ? detail::CigarOp::Equal : detail::CigarOp::Diff;
                
                if (!sub_op_set) {
                    sub_op = (q_pos >= 0 && r_pos >= 0 && query[query_start + q_pos] == ref[ref_start + r_pos]) ? detail::CigarOp::Equal : detail::CigarOp::Diff;
                    sub_op_set = true;
                    sub_len = 1;
                } else if (this_op == sub_op) {
                    ++sub_len;
                } else {
                    // Op changed, emit previous sub-run
                    new_cigar_rev.push_back(detail::encode_cigar(sub_len, sub_op));
                    ++new_n_cigar;
                    sub_op = this_op;
                    sub_len = 1;
                }
                
                --q_pos;
                --r_pos;
            }
            
            // Emit last sub-run
            if (sub_op_set) {
                new_cigar_rev.push_back(detail::encode_cigar(sub_len, sub_op));
                ++new_n_cigar;
            }
            
            q_pos -= len;
            r_pos -= len;
        } else if (op == detail::CigarOp::Ins) {
            new_cigar_rev.push_back(c);
            ++new_n_cigar;
            q_pos -= len;
        } else if (op == detail::CigarOp::Del) {
            new_cigar_rev.push_back(c);
            ++new_n_cigar;
            r_pos -= len;
        }
    }
    
    // Reverse to get correct order
    std::reverse(new_cigar_rev.begin(), new_cigar_rev.end());
    aln.cigar = std::move(new_cigar_rev);
    aln.n_cigar = new_n_cigar;

    // Count gaps and mismatches
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
        } else if (op == detail::CigarOp::Diff) {
            aln.n_mismatch += len;
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

    // DP arrays (full 2D for traceback)
    using DPState = struct { int32_t h, e, f; };
    std::vector<DPState> dp((qlen + 1) * bw);
    std::vector<uint8_t> trace((qlen + 1) * bw); // 0=None, 1=Match, 2=Insert, 3=Delete

    auto idx = [&](int32_t i, int32_t j) -> int32_t {
        return i * bw + (j - i + w);
    };

    // Initialize
    dp[w].h = 0;
    dp[w].e = dp[w].f = std::numeric_limits<int32_t>::min() / 2;

    for (int32_t j = 1; j <= rlen && j <= w; ++j) {
        int32_t init_idx = j + w;  // idx(0, j) = j + w
        dp[init_idx].h = sc.gap_open + sc.gap_ext * j;
        dp[init_idx].e = dp[init_idx].h;
        dp[init_idx].f = std::numeric_limits<int32_t>::min() / 2;
        trace[init_idx] = 2; // Insert
    }
    for (int32_t i = 1; i <= qlen && i <= w; ++i) {
        int32_t init_idx = i * bw - i + w;  // idx(i, 0) = i*bw - i + w
        dp[init_idx].h = sc.gap_open + sc.gap_ext * i;
        dp[init_idx].f = dp[init_idx].h;
        dp[init_idx].e = std::numeric_limits<int32_t>::min() / 2;
        trace[init_idx] = 3; // Delete
    }

    for (int32_t i = 1; i <= qlen; ++i) {
        int32_t j_min = std::max(1, i - w);
        int32_t j_max = std::min(rlen, i + w);

        for (int32_t j = j_min; j <= j_max; ++j) {
            int32_t cur_idx = idx(i, j);
            int32_t idx_diag = idx(i - 1, j - 1);
            int32_t idx_up = idx(i - 1, j);
            int32_t idx_left = idx(i, j - 1);

            int32_t s = (query[i - 1] == ref[j - 1] && query[i - 1] < 4) ? sc.match : sc.mismatch;

            int32_t h_diag = dp[idx_diag].h + s;
            int32_t e = std::max(dp[idx(i, j - 1)].h + sc.gap_open + sc.gap_ext,
                                  dp[idx(i, j - 1)].e + sc.gap_ext);
            int32_t f = std::max(dp[idx_up].h + sc.gap_open + sc.gap_ext,
                                  dp[idx(i - 1, j)].f + sc.gap_ext);
            int32_t h = std::max({h_diag, e, f});

            dp[idx(i, j)].h = h;
            dp[idx(i, j)].e = e;
            dp[idx(i, j)].f = f;

            // Determine traceback direction
            if (h == h_diag) {
                trace[cur_idx] = 1; // Match
            } else if (h == e) {
                trace[cur_idx] = 2; // Insert
            } else {
                trace[cur_idx] = 3; // Delete
            }

            dp[cur_idx].h = h;
            dp[cur_idx].e = e;
            dp[cur_idx].f = f;
        }
    }

    int32_t end_idx = qlen * bw + (rlen - qlen + w);
    int32_t max_score = dp[end_idx].h;
    int32_t ti = qlen, tj = rlen;
    uint8_t prev_dir = 0; // 0=None
    int32_t cur_len = 0;

    auto idx_func = [&](int32_t ti, int32_t tj) -> int32_t {
        return ti * bw + (tj - ti + w);
    };

    detail::CigarOp prev_cigar_op = detail::CigarOp::Match;

    while (ti > 0 || tj > 0) {
        uint8_t dir = trace[idx_func(ti, tj)];

        if (dir == 0 || (ti == 0 && tj == 0)) break;

        detail::CigarOp this_cigar_op;
        if (dir == 1) { // Match
            // Check if match or mismatch at this position
            bool is_match = (query[ti - 1] == ref[tj - 1] && query[ti - 1] < 4);
            this_cigar_op = is_match ? detail::CigarOp::Equal : detail::CigarOp::Diff;
        } else if (dir == 2) { // Insert
            this_cigar_op = detail::CigarOp::Ins;
        } else if (dir == 3) { // Delete
            this_cigar_op = detail::CigarOp::Del;
        } else {
            break;
        }

        if (this_cigar_op == prev_cigar_op && cur_len < 0xFFF) {
            ++cur_len;
        } else {
            if (cur_len > 0) {
                aln.cigar.push_back(detail::encode_cigar(cur_len, prev_cigar_op));
                aln.n_cigar++;
            }
            prev_cigar_op = this_cigar_op;
            cur_len = 1;
        }

        // Move according to traceback direction
        if (dir == 1) { // Match
            --ti; --tj;
        } else if (dir == 2) { // Insert
            --tj;
        } else if (dir == 3) { // Delete
            --ti;
        } else {
            break;
        }
    }

    if (cur_len > 0) {
        aln.cigar.push_back(detail::encode_cigar(cur_len, prev_cigar_op));
        aln.n_cigar++;
    }

    // Reverse CIGAR to get correct order
    std::reverse(aln.cigar.begin(), aln.cigar.end());

    // Count gaps and mismatches
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
        } else if (op == detail::CigarOp::Diff) {
            aln.n_mismatch += len;
        }
    }

    aln.score = max_score;
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