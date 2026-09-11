#pragma once

#include <span>
#include <cstdint>
#include <array>
#include <algorithm>
#include <limits>
#include <vector>
#include <numeric>

#include <bwa/index/fm_index.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>

namespace bwa::align {

// Maximal Exact Match (MEM)
struct MEM {
    int32_t query_pos = 0;    // Start position in query (0-based)
    int32_t ref_pos = 0;      // Start position in reference (0-based, local to ref_id)
    int32_t ref_id = 0;       // Reference sequence index (for multi-ref indexes)
    int32_t len = 0;          // Length of match
    int32_t score = 0;        // Score (len * match_score)
    bool is_forward = true;   // Strand

    [[nodiscard]] int32_t query_end() const noexcept { return query_pos + len; }
    [[nodiscard]] int32_t ref_end() const noexcept { return ref_pos + len; }
    [[nodiscard]] bool contains(const MEM& other) const noexcept {
        return ref_id == other.ref_id &&
               query_pos <= other.query_pos && query_end() >= other.query_end() &&
               ref_pos <= other.ref_pos && ref_end() >= other.ref_end();
    }
};

// MEM finder using FM-index
class MEMFinder {
    const index::FMIndex& fm_index_;
    int min_seed_len_ = 19;
    int max_occ_ = 500;        // Skip seeds with > this many occurrences
    int match_score_ = 1;

public:
    MEMFinder(const index::FMIndex& idx,
              int min_seed_len = 19,
              int max_occ = 500,
              int match_score = 1)
        : fm_index_(idx), min_seed_len_(min_seed_len), max_occ_(max_occ), match_score_(match_score) {}

    // Find all MEMs in query (forward strand)
    core::Vector<MEM> find(const std::span<const uint8_t>& query,
                     memory::Arena& arena) const {
        core::Vector<MEM> mems(&arena);
        find(query, mems);
        return mems;
    }

    // Fraction of query covered by exact MEMs of one strand.
    static double strand_coverage(const core::Vector<MEM>& mems, bool is_fwd,
                                  int32_t qlen) {
        if (qlen <= 0) return 0.0;
        // n is small (tens of MEMs); copy spans and sort for a union sweep.
        std::vector<std::pair<int32_t, int32_t>> spans;
        for (const auto& m : mems) {
            if (m.is_forward != is_fwd) continue;
            spans.emplace_back(m.query_pos, m.query_end());
        }
        if (spans.empty()) return 0.0;
        std::sort(spans.begin(), spans.end());
        int32_t total = 0, cur_end = -1;
        for (const auto& [b, e] : spans) {
            if (b > cur_end) {
                total += e - b;
                cur_end = e;
            } else if (e > cur_end) {
                total += e - cur_end;
                cur_end = e;
            }
        }
        return static_cast<double>(total) / static_cast<double>(qlen);
    }

    void find(const std::span<const uint8_t>& query,
              core::Vector<MEM>& mems) const {
        mems.clear();
        if (query.empty()) return;
        int32_t qlen = static_cast<int32_t>(query.size());

        // Forward search
        find_strand(query, mems, true);

        // Inexact (1-mismatch) seeding is ~1000x the cost of exact search
        // per position. BWA only reseeds poorly covered regions, so skip it
        // when exact MEMs already cover most of the query.
        constexpr double kMinExactCoverage = 0.5;
        if (strand_coverage(mems, true, qlen) < kMinExactCoverage) {
            find_inexact_seeds(query, mems, true);
        }

        // Reverse complement search
        core::Vector<uint8_t> rc_query;
        rc_query.resize(query.size());
        for (size_t i = 0; i < query.size(); ++i) {
            uint8_t b = query[query.size() - 1 - i];
            // Complement: A<->T, C<->G, N->N
            if (b == 0) rc_query[i] = 3;      // A -> T
            else if (b == 1) rc_query[i] = 2;  // C -> G
            else if (b == 2) rc_query[i] = 1;  // G -> C
            else if (b == 3) rc_query[i] = 0;  // T -> A
            else rc_query[i] = 4;              // N -> N
        }
        find_strand(std::span<const uint8_t>(rc_query.data(), rc_query.size()),
                    mems, false);
        if (strand_coverage(mems, false, qlen) < kMinExactCoverage) {
            find_inexact_seeds(std::span<const uint8_t>(rc_query.data(), rc_query.size()),
                               mems, false);
        }
    }

    // Rescue scan: find seeds with very short k-mers
    // Used as a last resort when no MEMs are found
    void rescue_scan(const std::span<const uint8_t>& query,
                     core::Vector<MEM>& mems,
                     bool is_fwd) const {
        int32_t qlen = static_cast<int32_t>(query.size());
        if (qlen < 4) return;

        // Try multiple short k values
        const int rescue_ks[] = {5, 4};
        int32_t max_occ_rescue = max_occ_ * 20; // Allow many occurrences for rescue

        for (int rescue_k : rescue_ks) {
            if (qlen < rescue_k) continue;

            for (int32_t i = 0; i + rescue_k <= qlen; ++i) {
                // Check if position has any N
                bool has_n = false;
                for (int k = 0; k < rescue_k; ++k) {
                    if (query[i + k] >= 4) { has_n = true; break; }
                }
                if (has_n) continue;

                // Find exact matches (full row interval [0, rows()))
                size_t l = 0, r = fm_index_.rows();
                for (int k = rescue_k - 1; k >= 0; --k) {
                    auto [nl, nr] = fm_index_.backward_extend(query[i + k], l, r);
                    if (nl >= nr) { l = 0; r = 0; break; }
                    l = nl; r = nr;
                }

                if (r > l && r - l <= static_cast<size_t>(max_occ_rescue)) {
                    emit_mems(mems, i, rescue_k, l, r, is_fwd);
                }
            }
        }
    }

    // Find inexact seeds (1 mismatch) for query positions
    void find_inexact_seeds(const std::span<const uint8_t>& query,
                            core::Vector<MEM>& mems,
                            bool is_fwd) const {
        int32_t qlen = static_cast<int32_t>(query.size());
        if (qlen < min_seed_len_) return;

        // For each position, check if there's already a MEM covering it
        // If not, try 1-mismatch seeds
        std::vector<bool> covered(qlen, false);
        for (const auto& mem : mems) {
            for (int32_t i = mem.query_pos; i < mem.query_end(); ++i) {
                if (i >= 0 && i < qlen) covered[i] = true;
            }
        }

        // Try 1-mismatch seeds at uncovered positions
        for (int32_t i = 0; i + min_seed_len_ <= qlen; ++i) {
            if (covered[i]) continue;

            // Check if position has any N
            bool has_n = false;
            for (int k = 0; k < min_seed_len_; ++k) {
                if (query[i + k] >= 4) { has_n = true; break; }
            }
            if (has_n) continue;

            // Try exact match first (full row interval [0, rows()))
            size_t l = 0, r = fm_index_.rows();
            for (int k = min_seed_len_ - 1; k >= 0; --k) {
                auto [nl, nr] = fm_index_.backward_extend(query[i + k], l, r);
                if (nl >= nr) { l = 0; r = 0; break; }
                l = nl; r = nr;
            }

            if (r > l && r - l <= static_cast<size_t>(max_occ_)) {
                // Exact match found (all occurrences, capped)
                size_t before = mems.size();
                emit_mems(mems, i, min_seed_len_, l, r, is_fwd);
                if (mems.size() > before) covered[i] = true;
            } else {
                // Try 1-mismatch seeds
                for (int mm_pos = 0; mm_pos < min_seed_len_; ++mm_pos) {
                    for (uint8_t alt = 0; alt < 4; ++alt) {
                        if (alt == query[i + mm_pos]) continue;

                        l = 0; r = fm_index_.rows();
                        bool ok = true;
                        for (int k = min_seed_len_ - 1; k >= 0; --k) {
                            uint8_t base = (k == mm_pos) ? alt : query[i + k];
                            auto [nl, nr] = fm_index_.backward_extend(base, l, r);
                            if (nl >= nr) { ok = false; break; }
                            l = nl; r = nr;
                        }
                        if (!ok) continue;
                        if (r - l > static_cast<size_t>(max_occ_)) continue;

                        size_t before = mems.size();
                        emit_mems(mems, i, min_seed_len_, l, r, is_fwd);
                        if (mems.size() > before) {
                            covered[i] = true;
                            break;  // Found matches for this position
                        }
                    }
                    if (covered[i]) break;
                }
            }
        }
    }

    // Max occurrences materialized per seed. BWA retrieves up to max_occ
    // positions; we cap lower because only the count (for MAPQ ambiguity)
    // and a few representatives (for chaining) matter, while the DP cost
    // grows with the square of the MEM count.
    static constexpr int32_t kMaxOccPerSeed = 32;
    // Hard cap on MEMs per find() call to bound chaining cost on
    // pathological repeat reads.
    static constexpr size_t kMaxMemsTotal = 4096;
    // Only seeds at least this multiple of min_seed_len fan out to multiple
    // occurrences. Short (e.g. 9-mer inexact) seeds in repeats would
    // otherwise multiply MEM counts for weak evidence; long multi-occurrence
    // seeds are the meaningful repeat indicators that drive MAPQ.
    static constexpr int32_t kFanoutLenFactor = 2;

    // Emit one MEM per located occurrence (up to kMaxOccPerSeed rows),
    // sharing the same query span. Different loci become separate chains,
    // which is what drives MAPQ down for repeats (BWA sub_n equivalent).
    void emit_mems(core::Vector<MEM>& mems, int32_t query_pos, int32_t len,
                   size_t row_l, size_t row_r, bool is_fwd) const {
        size_t occ = row_r > row_l ? row_r - row_l : 0;
        size_t n = occ;
        if (len < kFanoutLenFactor * min_seed_len_) {
            // Short seed: single representative (avoids fan-out blowup).
            n = occ > 0 ? 1 : 0;
        } else if (n > static_cast<size_t>(kMaxOccPerSeed)) {
            n = static_cast<size_t>(kMaxOccPerSeed);
        }
        for (size_t k = 0; k < n; ++k) {
            if (mems.size() >= kMaxMemsTotal) return;
            auto pos = fm_index_.locate(row_l + k);
            if (!pos) continue;
            MEM mem;
            mem.query_pos = query_pos;
            mem.ref_pos = static_cast<int32_t>(*pos);
            mem.len = len;
            mem.score = len * match_score_;
            mem.is_forward = is_fwd;
            mems.push_back(mem);
        }
    }

    // Find MEMs using SMEM (simple maximal exact matches) algorithm
    void find_strand(const std::span<const uint8_t>& query,
                     core::Vector<MEM>& mems,
                     bool is_fwd) const {        int32_t qlen = static_cast<int32_t>(query.size());
        if (qlen < min_seed_len_) return;

        // Use backward search: process query from right to left
        int32_t i = qlen - 1;
        while (i >= 0) {
            // Skip N's (from right)
            while (i >= 0 && query[i] >= 4) --i;
            if (i < 0) break;

            // Extend backward as long as possible (full row interval)
            int32_t l = 0, r = static_cast<int32_t>(fm_index_.rows());
            int32_t j = i;
            int32_t best_l = l, best_r = r, best_j = j;

            while (j >= 0 && query[j] < 4) {
                uint8_t base = query[j];
                auto [new_l, new_r] = fm_index_.backward_extend(base, l, r);
                if (new_l >= new_r) break;
                l = new_l;
                r = new_r;
                --j;

                // Check if this is a MEM (occurrence count is acceptable)
                int32_t occ = r - l;
                if (occ <= max_occ_) {
                    best_l = l; best_r = r; best_j = j;
                }
            }

            int32_t mem_len = i - best_j;
            if (mem_len >= min_seed_len_) {
                // Retrieve all occurrences (capped): repeat copies become
                // separate MEMs at distinct loci instead of one arbitrary hit.
                emit_mems(mems, best_j + 1, mem_len,
                          static_cast<size_t>(best_l),
                          static_cast<size_t>(best_r), is_fwd);
            }

            // Move to next position (from right)
            if (j == i) --i;
            else i = j;
        }
    }

    // Find MEMs with overlap filtering
    core::Vector<MEM> find_filtered(const std::span<const uint8_t>& query,
                              memory::Arena& arena) const {
        core::Vector<MEM> mems = find(query, arena);
        filter_overlaps(mems);
        return mems;
    }

    // Remove contained MEMs and overlapping MEMs with lower score
    static void filter_overlaps(core::Vector<MEM>& mems) {
        if (mems.size() <= 1) return;

        // Sort by query position, then by length descending (ref_id for determinism)
        core::sort(mems.begin(), mems.end(),
                   [](const MEM& a, const MEM& b) {
                        if (a.query_pos != b.query_pos)
                            return a.query_pos < b.query_pos;
                        if (a.ref_id != b.ref_id)
                            return a.ref_id < b.ref_id;
                        return a.len > b.len;
                    });

        core::Vector<MEM> filtered;
        filtered.reserve(mems.size());

        for (const MEM& mem : mems) {
            bool contained = false;
            for (const MEM& kept : filtered) {
                // Only suppress MEMs on the same reference: MEMs on different
                // references are alternative mappings, not duplicates.
                if (kept.ref_id != mem.ref_id) continue;
                if (kept.contains(mem)) {
                    contained = true;
                    break;
                }
                // Suppress on query overlap only when the two MEMs also
                // overlap in reference (same locus). Same-query-span MEMs at
                // distinct loci are repeat copies: keeping them is what lets
                // chaining report competing placements and MAPQ reflect the
                // ambiguity (BWA sub_n equivalent). Previously the reference
                // check was missing, so repeats collapsed to one arbitrary
                // copy reported with MAPQ 60.
                if (mem.query_pos < kept.query_end() && mem.query_end() > kept.query_pos &&
                    mem.ref_pos < kept.ref_end() && mem.ref_end() > kept.ref_pos) {
                    // Overlap - keep higher score
                    if (mem.score <= kept.score) {
                        contained = true;
                        break;
                    }
                }
            }
            if (!contained) {
                // Remove any kept that are contained by this (same reference only)
                auto it = filtered.begin();
                while (it != filtered.end()) {
                    if (it->ref_id == mem.ref_id && mem.contains(*it)) {
                        it = filtered.erase(it);
                    } else {
                        ++it;
                    }
                }
                filtered.push_back(mem);
            }
        }

        mems.swap(filtered);
    }

    // Chain MEMs colinearly (for chaining step)
    struct Chain {
        core::Vector<MEM> mems;
        int32_t score = 0;
        int32_t ref_id = 0;   // Reference all MEMs in this chain belong to
        int32_t query_begin = 0, query_end = 0;
        int32_t ref_begin = 0, ref_end = 0;
    };

    // DP-optimal chaining (BWA-MEM style colinear chaining).
    // Within each reference, seeds sorted by (query_pos, ref_pos) are linked
    // by longest-path DP: f[i] = score[i] + max(0, max over predecessors j of
    // f[j] - gap_cost), with gap_cost = max(query_gap, ref_gap) in match-score
    // units. A transition j->i requires strict colinearity (non-negative gaps
    // within max_gap); overlapping seeds are left to filter_overlaps. Chains
    // are extracted best-first (traceback stops at seeds already emitted) so
    // chains[0] is the primary candidate and later chains are secondaries.
    // Returned chains are sorted by score descending. O(n^2) per reference
    // group; n is the MEM count per read (tens to low hundreds).
    // MEMs are grouped by reference: chains never span references.
    core::Vector<Chain> chain(const core::Vector<MEM>& mems,
                        int32_t max_gap = 10000,
                        int32_t min_chain_score = 30) const {
        core::Vector<Chain> chains;
        if (mems.empty()) return chains;

        // Copy MEMs to a sortable container
        std::vector<MEM> sorted_mems;
        sorted_mems.reserve(mems.size());
        for (size_t i = 0; i < mems.size(); ++i) {
            sorted_mems.push_back(mems[i]);
        }
        std::sort(sorted_mems.begin(), sorted_mems.end(),
                  [](const MEM& a, const MEM& b) {
                      if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
                      if (a.query_pos != b.query_pos) return a.query_pos < b.query_pos;
                      return a.ref_pos < b.ref_pos;
                  });

        // DP independently per reference group
        size_t g_begin = 0;
        const size_t total = sorted_mems.size();
        // Scratch buffers reused across groups (sized to largest group)
        std::vector<int32_t> best;
        std::vector<int32_t> parent;
        std::vector<char> used;
        while (g_begin < total) {
            size_t g_end = g_begin + 1;
            const int32_t group_ref = sorted_mems[g_begin].ref_id;
            while (g_end < total && sorted_mems[g_end].ref_id == group_ref) {
                ++g_end;
            }
            chain_group(sorted_mems, g_begin, g_end, max_gap, min_chain_score,
                        best, parent, used, chains);
            g_begin = g_end;
        }

        // Best chain first (aligner takes chains[0] as primary candidate).
        // Full tie-break for determinism: repeat copies with equal scores
        // must order identically across runs. stable_sort keeps DP discovery
        // order for full ties. Forward strand preferred (matches user
        // expectation for symmetric repeats; MAPQ still flags ambiguity).
        std::stable_sort(chains.begin(), chains.end(), [](const Chain& a, const Chain& b) {
            if (a.score != b.score) return a.score > b.score;
            if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
            if (a.ref_begin != b.ref_begin) return a.ref_begin < b.ref_begin;
            if (a.query_begin != b.query_begin) return a.query_begin < b.query_begin;
            bool a_fwd = !a.mems.empty() && a.mems.front().is_forward;
            bool b_fwd = !b.mems.empty() && b.mems.front().is_forward;
            return a_fwd > b_fwd;
        });

        return chains;
    }

    // Longest-path DP over one reference group [g_begin, g_end), appending
    // best-first chains to `chains`. Scratch vectors are resized as needed.
    static void chain_group(const std::vector<MEM>& mems,
                            size_t g_begin, size_t g_end,
                            int32_t max_gap, int32_t min_chain_score,
                            std::vector<int32_t>& best,
                            std::vector<int32_t>& parent,
                            std::vector<char>& used,
                            core::Vector<Chain>& chains) {
        const size_t n = g_end - g_begin;
        best.assign(n, 0);
        parent.assign(n, -1);
        used.assign(n, 0);

        for (size_t ii = 0; ii < n; ++ii) {
            const MEM& mem = mems[g_begin + ii];
            int32_t f = mem.score;
            int32_t p = -1;
            for (size_t jj = 0; jj < ii; ++jj) {
                const MEM& prev = mems[g_begin + jj];
                int32_t qgap = mem.query_pos - prev.query_end();
                int32_t rgap = mem.ref_pos - prev.ref_end();
                if (qgap < 0 || rgap < 0 || qgap > max_gap || rgap > max_gap) {
                    continue;
                }
                int32_t gap = qgap > rgap ? qgap : rgap;
                int32_t cand = best[jj] + mem.score - gap;
                if (cand > f) {
                    f = cand;
                    p = static_cast<int32_t>(jj);
                }
            }
            best[ii] = f;
            parent[ii] = p;
        }

        // Best-first extraction: take the highest-scoring unused endpoint,
        // trace back (stopping at already-emitted seeds), rescore the emitted
        // prefix from its members, and keep it if it passes the threshold.
        for (;;) {
            int32_t top = -1;
            for (size_t ii = 0; ii < n; ++ii) {
                if (!used[ii] && (top < 0 || best[ii] > best[static_cast<size_t>(top)])) {
                    top = static_cast<int32_t>(ii);
                }
            }
            if (top < 0 || best[static_cast<size_t>(top)] < min_chain_score) break;

            // Traceback (indices into mems), oldest seed last in path
            std::vector<size_t> path;
            for (int32_t k = top; k >= 0 && !used[static_cast<size_t>(k)];
                 k = parent[static_cast<size_t>(k)]) {
                path.push_back(g_begin + static_cast<size_t>(k));
            }
            if (path.empty()) break;  // cannot happen: top itself is unused
            for (size_t pi = 0; pi < path.size(); ++pi) {
                used[path[pi] - g_begin] = 1;
            }

            Chain ch;
            ch.ref_id = mems[path.back()].ref_id;
            ch.score = 0;
            bool first = true;
            int32_t prev_qend = 0, prev_rend = 0;
            for (size_t pi = path.size(); pi-- > 0;) {
                const MEM& mem = mems[path[pi]];
                ch.mems.push_back(mem);
                if (first) {
                    ch.score = mem.score;
                    ch.query_begin = mem.query_pos;
                    ch.ref_begin = mem.ref_pos;
                    first = false;
                } else {
                    int32_t qgap = mem.query_pos - prev_qend;
                    int32_t rgap = mem.ref_pos - prev_rend;
                    int32_t gap = qgap > rgap ? qgap : rgap;
                    ch.score += mem.score - (gap > 0 ? gap : 0);
                }
                prev_qend = mem.query_end();
                prev_rend = mem.ref_end();
            }
            ch.query_end = prev_qend;
            ch.ref_end = prev_rend;
            if (ch.score >= min_chain_score && !ch.mems.empty()) {
                chains.push_back(std::move(ch));
            }
        }
    }

    // Configuration
    [[nodiscard]] int min_seed_len() const noexcept { return min_seed_len_; }
    void set_min_seed_len(int v) noexcept { min_seed_len_ = v; }

    [[nodiscard]] int max_occ() const noexcept { return max_occ_; }
    void set_max_occ(int v) noexcept { max_occ_ = v; }

    [[nodiscard]] int match_score() const noexcept { return match_score_; }
    void set_match_score(int v) noexcept { match_score_ = v; }
};

// SMEM iterator for streaming MEM finding
// Mirrors MEMFinder::find_strand: FM backward search processes the query
// from right to left, so the iterator extends leftwards from the right end.
class SMEMIterator {
    const index::FMIndex& fm_index_;
    std::span<const uint8_t> query_;
    int min_len_;
    int max_occ_;
    int32_t i_ = -1;      // right end of current search region (inclusive)
    int32_t best_l_ = 0;  // FM interval of longest acceptable match
    int32_t best_j_ = 0;  // left end (inclusive) of longest acceptable match
    int32_t stop_j_ = 0;  // where backward extension stopped (for advancing)
    bool exhausted_ = true;

public:
    SMEMIterator(const index::FMIndex& idx,
                 std::span<const uint8_t> query,
                 int min_len = 19,
                 int max_occ = 500)
        : fm_index_(idx), query_(query), min_len_(min_len), max_occ_(max_occ) {
        if (!query_.empty()) {
            i_ = static_cast<int32_t>(query_.size()) - 1;
            exhausted_ = false;
            advance_to_next();
        }
    }

    [[nodiscard]] bool has_next() const noexcept { return !exhausted_; }

    std::optional<MEM> next() {
        if (exhausted_) return std::nullopt;

        MEM mem;
        mem.query_pos = best_j_ + 1;
        if (auto pos = fm_index_.locate(static_cast<size_t>(best_l_))) {
            mem.ref_pos = static_cast<int32_t>(*pos);
        }
        mem.len = i_ - best_j_;
        mem.score = mem.len;
        mem.is_forward = true;

        // Advance past this MEM (from the right), as in find_strand
        if (stop_j_ == i_) --i_;
        else i_ = stop_j_;
        advance_to_next();
        return mem;
    }

private:
    void advance_to_next() {
        int32_t qlen = static_cast<int32_t>(query_.size());

        while (i_ >= 0) {
            // Skip N's (from right)
            while (i_ >= 0 && query_[static_cast<size_t>(i_)] >= 4) --i_;
            if (i_ < 0) { exhausted_ = true; return; }

            // Extend backward as long as possible (full row interval)
            int32_t l = 0, r = static_cast<int32_t>(fm_index_.rows());
            int32_t j = i_;
            best_l_ = l; best_j_ = j;

            while (j >= 0 && query_[static_cast<size_t>(j)] < 4) {
                auto [new_l, new_r] = fm_index_.backward_extend(
                    query_[static_cast<size_t>(j)],
                    static_cast<size_t>(l), static_cast<size_t>(r));
                if (new_l >= new_r) break;
                l = static_cast<int32_t>(new_l);
                r = static_cast<int32_t>(new_r);
                --j;

                if (r - l <= max_occ_) {
                    best_l_ = l; best_j_ = j;
                }
            }
            stop_j_ = j;

            if (i_ - best_j_ >= min_len_) {
                return; // Found MEM
            }

            if (j == i_) --i_;
            else i_ = j;
            (void)qlen;
        }

        exhausted_ = true;
    }
};

} // namespace bwa::align