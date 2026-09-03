#pragma once

#include <bwa/index/fm_index.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <span>
#include <cstdint>
#include <array>
#include <algorithm>
#include <limits>
#include <vector>
#include <numeric>

namespace bwa::align {

// Maximal Exact Match (MEM)
struct MEM {
    int32_t query_pos = 0;    // Start position in query (0-based)
    int32_t ref_pos = 0;      // Start position in reference (0-based)
    int32_t len = 0;          // Length of match
    int32_t score = 0;        // Score (len * match_score)
    bool is_forward = true;   // Strand

    [[nodiscard]] int32_t query_end() const noexcept { return query_pos + len; }
    [[nodiscard]] int32_t ref_end() const noexcept { return ref_pos + len; }
    [[nodiscard]] bool contains(const MEM& other) const noexcept {
        return query_pos <= other.query_pos && query_end() >= other.query_end() &&
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

    void find(const std::span<const uint8_t>& query,
              core::Vector<MEM>& mems) const {
        mems.clear();
        if (query.empty()) return;

        // Forward search
        find_strand(query, mems, true);

        // Find inexact seeds (1 mismatch) for positions without exact MEMs
        find_inexact_seeds(query, mems, true);

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
        find_inexact_seeds(std::span<const uint8_t>(rc_query.data(), rc_query.size()),
                           mems, false);
    }

    // Rescue scan: find seeds with very short k-mers (k=5)
    // Used as a last resort when no MEMs are found
    void rescue_scan(const std::span<const uint8_t>& query,
                     core::Vector<MEM>& mems,
                     bool is_fwd) const {
        int32_t qlen = static_cast<int32_t>(query.size());
        if (qlen < 5) return;

        // Try k=5 seeds
        const int rescue_k = 5;
        int32_t max_occ_rescue = max_occ_ * 10; // Allow more occurrences for rescue

        for (int32_t i = 0; i + rescue_k <= qlen; ++i) {
            // Check if position has any N
            bool has_n = false;
            for (int k = 0; k < rescue_k; ++k) {
                if (query[i + k] >= 4) { has_n = true; break; }
            }
            if (has_n) continue;

            // Find exact matches
            size_t l = 0, r = fm_index_.length();
            for (int k = rescue_k - 1; k >= 0; --k) {
                auto [nl, nr] = fm_index_.backward_extend(query[i + k], l, r);
                if (nl >= nr) { l = 0; r = 0; break; }
                l = nl; r = nr;
            }

            if (r > l && r - l <= static_cast<size_t>(max_occ_rescue)) {
                if (auto pos = fm_index_.locate(l)) {
                    MEM mem;
                    mem.query_pos = i;
                    mem.ref_pos = *pos;
                    mem.len = rescue_k;
                    mem.score = rescue_k * match_score_;
                    mem.is_forward = is_fwd;
                    mems.push_back(mem);
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

            // Try exact match first
            size_t l = 0, r = fm_index_.length();
            for (int k = min_seed_len_ - 1; k >= 0; --k) {
                auto [nl, nr] = fm_index_.backward_extend(query[i + k], l, r);
                if (nl >= nr) { l = 0; r = 0; break; }
                l = nl; r = nr;
            }

            if (r > l && r - l <= static_cast<size_t>(max_occ_)) {
                // Exact match found
                if (auto pos = fm_index_.locate(l)) {
                    MEM mem;
                    mem.query_pos = i;
                    mem.ref_pos = *pos;
                    mem.len = min_seed_len_;
                    mem.score = min_seed_len_ * match_score_;
                    mem.is_forward = is_fwd;
                    mems.push_back(mem);
                    covered[i] = true;
                }
            } else {
                // Try 1-mismatch seeds
                for (int mm_pos = 0; mm_pos < min_seed_len_; ++mm_pos) {
                    for (uint8_t alt = 0; alt < 4; ++alt) {
                        if (alt == query[i + mm_pos]) continue;

                        l = 0; r = fm_index_.length();
                        bool ok = true;
                        for (int k = min_seed_len_ - 1; k >= 0; --k) {
                            uint8_t base = (k == mm_pos) ? alt : query[i + k];
                            auto [nl, nr] = fm_index_.backward_extend(base, l, r);
                            if (nl >= nr) { ok = false; break; }
                            l = nl; r = nr;
                        }
                        if (!ok) continue;
                        if (r - l > static_cast<size_t>(max_occ_)) continue;

                        if (auto pos = fm_index_.locate(l)) {
                            MEM mem;
                            mem.query_pos = i;
                            mem.ref_pos = *pos;
                            mem.len = min_seed_len_;
                            mem.score = min_seed_len_ * match_score_;
                            mem.is_forward = is_fwd;
                            mems.push_back(mem);
                            covered[i] = true;
                            break;  // Found a match for this position
                        }
                    }
                    if (covered[i]) break;
                }
            }
        }
    }

    // Find MEMs using SMEM (simple maximal exact matches) algorithm
    void find_strand(const std::span<const uint8_t>& query,
                     core::Vector<MEM>& mems,
                     bool is_fwd) const {
        int32_t qlen = static_cast<int32_t>(query.size());
        if (qlen < min_seed_len_) return;

        // Use backward search: process query from right to left
        int32_t i = qlen - 1;
        while (i >= 0) {
            // Skip N's (from right)
            while (i >= 0 && query[i] >= 4) --i;
            if (i < 0) break;

            // Extend backward as long as possible
            int32_t l = 0, r = fm_index_.length();
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
                // Get reference position (sampled SA)
                if (auto pos = fm_index_.locate(best_l)) {
                    MEM mem;
                    mem.query_pos = best_j + 1;
                    mem.ref_pos = *pos;
                    mem.len = mem_len;
                    mem.score = mem_len * match_score_;
                    mem.is_forward = is_fwd;
                    mems.push_back(mem);
                }
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

        // Sort by query position, then by length descending
        core::sort(mems.begin(), mems.end(),
                   [](const MEM& a, const MEM& b) {
                       if (a.query_pos != b.query_pos)
                           return a.query_pos < b.query_pos;
                       return a.len > b.len;
                   });

        core::Vector<MEM> filtered;
        filtered.reserve(mems.size());

        for (const MEM& mem : mems) {
            bool contained = false;
            for (const MEM& kept : filtered) {
                if (kept.contains(mem)) {
                    contained = true;
                    break;
                }
                // Check overlap
                if (mem.query_pos < kept.query_end() && mem.query_end() > kept.query_pos) {
                    // Overlap - keep higher score
                    if (mem.score <= kept.score) {
                        contained = true;
                        break;
                    }
                }
            }
            if (!contained) {
                // Remove any kept that are contained by this
                auto it = filtered.begin();
                while (it != filtered.end()) {
                    if (mem.contains(*it)) {
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
        int32_t query_begin = 0, query_end = 0;
        int32_t ref_begin = 0, ref_end = 0;
    };

    // DP-based optimal chaining (BWA-MEM style)
    // Finds the highest-scoring chain of colinear MEMs
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
                      if (a.query_pos != b.query_pos) return a.query_pos < b.query_pos;
                      return a.ref_pos < b.ref_pos;
                  });

        int32_t n = static_cast<int32_t>(sorted_mems.size());
        if (n == 0) return chains;

        // Simple greedy chaining (avoids memory issues with DP)
        Chain current;
        current.mems.push_back(sorted_mems[0]);
        current.score = sorted_mems[0].score;
        current.query_begin = current.query_end = sorted_mems[0].query_end();
        current.ref_begin = current.ref_end = sorted_mems[0].ref_end();

        for (int32_t i = 1; i < n; ++i) {
            const MEM& mem = sorted_mems[i];
            const MEM& last = current.mems.back();

            int32_t qgap = mem.query_pos - last.query_end();
            int32_t rgap = mem.ref_pos - last.ref_end();

            if (qgap >= 0 && rgap >= 0 && qgap <= max_gap && rgap <= max_gap) {
                current.mems.push_back(mem);
                current.score += mem.score - std::max(qgap, rgap);
                current.query_end = mem.query_end();
                current.ref_end = mem.ref_end();
            } else {
                if (current.score >= min_chain_score) {
                    chains.push_back(std::move(current));
                }
                current = Chain{};
                current.mems.push_back(mem);
                current.score = mem.score;
                current.query_begin = current.query_end = mem.query_end();
                current.ref_begin = current.ref_end = mem.ref_end();
            }
        }

        if (current.score >= min_chain_score) {
            chains.push_back(std::move(current));
        }

        return chains;
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
class SMEMIterator {
    const index::FMIndex& fm_index_;
    std::span<const uint8_t> query_;
    int min_len_;
    int max_occ_;
    int32_t i_ = 0;
    int32_t l_ = 0, r_ = 0;
    int32_t best_l_ = 0, best_r_ = 0, best_j_ = 0;
    bool exhausted_ = false;

public:
    SMEMIterator(const index::FMIndex& idx,
                 std::span<const uint8_t> query,
                 int min_len = 19,
                 int max_occ = 500)
        : fm_index_(idx), query_(query), min_len_(min_len), max_occ_(max_occ) {
        if (!query_.empty()) {
            l_ = 0; r_ = fm_index_.length();
            advance_to_next();
        }
    }

    [[nodiscard]] bool has_next() const noexcept { return !exhausted_; }

    std::optional<MEM> next() {
        if (exhausted_) return std::nullopt;

        int32_t mem_len = best_j_ - i_;
        MEM mem;
        mem.query_pos = i_;
        if (auto pos = fm_index_.locate(best_l_)) {
            mem.ref_pos = *pos;
        }
        mem.len = mem_len;
        mem.score = mem_len;
        mem.is_forward = true;

        advance_to_next();
        return mem;
    }

private:
    void advance_to_next() {
        int32_t qlen = static_cast<int32_t>(query_.size());

        while (i_ < qlen) {
            // Skip N's
            while (i_ < qlen && query_[i_] >= 4) ++i_;
            if (i_ >= qlen) { exhausted_ = true; return; }

            l_ = 0; r_ = fm_index_.length();
            int32_t j = i_;
            best_l_ = l_; best_r_ = r_; best_j_ = j;

            while (j < qlen && query_[j] < 4) {
                auto [new_l, new_r] = fm_index_.backward_extend(query_[j], l_, r_);
                if (new_l >= new_r) break;
                l_ = new_l; r_ = new_r; ++j;

                if (r_ - l_ <= max_occ_) {
                    best_l_ = l_; best_r_ = r_; best_j_ = j;
                }
            }

            if (best_j_ - i_ >= min_len_) {
                return; // Found MEM
            }

            if (j == i_) ++i_;
            else i_ = j;
        }

        exhausted_ = true;
    }
};

} // namespace bwa::align