#include <ranges>  // Must be first to avoid namespace pollution in C++26
#include <bwa/pipeline.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/hash_map.hpp>
#include <bwa/core/sort.hpp>
#include <bwa/index/fm_index.hpp>
#include <bwa/align/mem.hpp>
#include <bwa/align/sw.hpp>
#include <bwa/io/seq_io.hpp>
#include <bwa/io/bam_io.hpp>
#include <bwa/io/bai_io.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <map>
#include <sstream>

namespace bwa {

namespace {

// Reverse complement of 2-bit base code (0=A,1=C,2=G,3=T,4+=N)
inline uint8_t rc_base(uint8_t b) noexcept {
    return b == 0 ? 3 : b == 1 ? 2 : b == 2 ? 1 : b == 3 ? 0 : 4;
}

// Strand vote for a chain: true = reverse (majority of MEM score on RC).
inline bool chain_is_reverse(const bwa::align::MEMFinder::Chain& chain) {
    int32_t fwd = 0, rev = 0;
    for (size_t i = 0; i < chain.mems.size(); ++i) {
        if (chain.mems[i].is_forward) fwd += chain.mems[i].score;
        else rev += chain.mems[i].score;
    }
    return rev > fwd;
}

// Helper: format SAM CIGAR from encoded cigar vector
std::string format_cigar(const std::vector<uint32_t>& cigar) {
    std::string s;
    for (uint32_t c : cigar) {
        int len = bwa::align::cigar_len(c);
        char op = bwa::align::cigar_char(bwa::align::cigar_op(c));
        s += std::to_string(len);
        s += op;
    }
    return s;
}

} // anonymous namespace

void Aligner::align_impl(const bwa::io::SeqRecord& read, AlignmentResult& result) const {
    result.clear();

    // Reset arena for this read
    bwa::memory::reset_tls_arena();

    // Determine effective config based on read length
    Config effective_config = config_;
    int32_t query_len = static_cast<int32_t>(read.seq.size());
    if (query_len > 250 && !config_.read_group.has_value()) {
        // Auto-switch to long-read config for very long reads
        effective_config = Config::long_reads();
    }

    // Pack query sequence
    bwa::index::PackedSequence query;
    query.append(read.seq.data(), read.seq.size());

    if (query.size() == 0) {
        result.mapped = false;
        return;
    }

    // Empty index: nothing to align against
    if (index_.num_references() == 0) {
        result.mapped = false;
        return;
    }

    // Find MEMs on every reference (use TLS arena for memory).
    // Previously only fm_index()[0] was searched, making all other
    // references invisible. Each MEM is tagged with its ref_id.
    bwa::memory::Arena& arena = bwa::memory::get_tls_arena();
    bwa::core::Vector<bwa::align::MEM> mems(&arena);
    const size_t num_refs = index_.num_references();
    std::vector<size_t> ref_mem_counts(num_refs, 0);
    for (size_t ri = 0; ri < num_refs; ++ri) {
        if (ri >= index_.fm_index().num_sequences()) continue;
        if (index_.references()[ri].length == 0) continue;
        bwa::align::MEMFinder finder(index_.fm_index()[ri],
                                     mem_finder_.min_seed_len(),
                                     mem_finder_.max_occ(),
                                     mem_finder_.match_score());
        auto found = finder.find(query.bases(), arena);
        for (size_t i = 0; i < found.size(); ++i) {
            found[i].ref_id = static_cast<int32_t>(ri);
            mems.push_back(found[i]);
        }
        ref_mem_counts[ri] = found.size();
    }
    bwa::align::MEMFinder::filter_overlaps(mems);

    // Rescue scan per reference lacking MEMs
    for (size_t ri = 0; ri < num_refs; ++ri) {
        if (ref_mem_counts[ri] > 0) continue;
        if (ri >= index_.fm_index().num_sequences()) continue;
        if (index_.references()[ri].length == 0) continue;
        bwa::align::MEMFinder finder(index_.fm_index()[ri],
                                     mem_finder_.min_seed_len(),
                                     mem_finder_.max_occ(),
                                     mem_finder_.match_score());
        size_t before = mems.size();
        finder.rescue_scan(query.bases(), mems, true);
        for (size_t i = before; i < mems.size(); ++i) {
            mems[i].ref_id = static_cast<int32_t>(ri);
        }
    }

    if (mems.empty()) {
        result.mapped = false;
        return;
    }

    // Chain MEMs
    auto chains = mem_finder_.chain(mems, effective_config.max_gap, effective_config.min_chain_score);

    if (chains.empty()) {
        result.mapped = false;
        return;
    }

    // Process best chain
    const auto& best_chain = chains[0];

    // Process secondary chains for XS tag (suboptimal alignment score)
    int32_t suboptimal_score = 0;
    if (chains.size() > 1) {
        suboptimal_score = chains[1].score;
    }

    // Helper: reference + global-coordinate window for a chain.
    // MEM ref_pos is local to its reference; extract_ref takes global
    // concatenated coordinates (ref.offset + local).
    auto chain_ref_window = [&](const bwa::align::MEMFinder::Chain& chain,
                                int32_t& ref_id_out,
                                int32_t& ref_begin_out,
                                int32_t& ref_end_out) -> bool {
        ref_id_out = chain.mems.front().ref_id;
        if (ref_id_out < 0 ||
            static_cast<size_t>(ref_id_out) >= index_.num_references()) {
            return false;
        }
        const auto& ref = index_.references()[static_cast<size_t>(ref_id_out)];
        int32_t ref_len = static_cast<int32_t>(ref.length);
        int32_t global_base = static_cast<int32_t>(ref.offset);
        int32_t local_begin = chain.mems.front().ref_pos;
        int32_t padding = effective_config.band_width * 2;
        ref_begin_out = global_base + std::max<int32_t>(0, local_begin - padding);
        ref_end_out = global_base + std::min<int32_t>(ref_len, local_begin + query_len + padding);
        return ref_end_out > ref_begin_out;
    };

    // Pack query for SW (forward + reverse-complement for reverse chains)
    std::vector<uint8_t> query_bytes(query_len);
    std::vector<uint8_t> query_rc(query_len);
    auto qbases = query.bases();
    for (int32_t i = 0; i < query_len; ++i) {
        query_bytes[i] = qbases[i];
        query_rc[query_len - 1 - i] = rc_base(qbases[i]);
    }

    // Adaptive scoring based on read length (BWA-MEM style)
    // For longer reads, increase gap penalties to avoid spurious alignments
    bwa::align::Scoring adaptive_scoring = effective_config.scoring;
    if (query_len > 100) {
        // Scale gap penalties with read length
        adaptive_scoring.gap_open = std::min<int>(-1, static_cast<int>(-0.01 * query_len - 4));
        adaptive_scoring.gap_ext = std::min<int>(-1, static_cast<int>(-0.1 * query_len));
    }

    // Process up to max_secondary alignments (best chain + secondary chains)
    int max_alignments = std::min<int>(chains.size(), 1 + effective_config.max_occ / 100);
    if (max_alignments < 1) max_alignments = 1;
    if (max_alignments > 4) max_alignments = 4; // Cap at 4 total alignments

    bwa::align::Alignment best_swaln;
    int32_t best_chain_idx = -1;
    std::vector<std::pair<int32_t, bwa::align::Alignment>> all_alignments;

    std::vector<char> chain_rev(chains.size(), 0);
    for (size_t ci = 0; ci < chains.size(); ++ci) {
        chain_rev[ci] = chain_is_reverse(chains[ci]) ? 1 : 0;
    }

    for (int32_t ci = 0; ci < max_alignments && ci < (int32_t)chains.size(); ++ci) {
        const auto& chain = chains[ci];
        int32_t ref_id = 0, ref_begin = 0, ref_end = 0;
        if (!chain_ref_window(chain, ref_id, ref_begin, ref_end)) continue;

        auto ref_region = index_.extract_ref(ref_begin, ref_end);
        if (ref_region.empty()) continue;

        const uint8_t* qptr = chain_rev[ci] ? query_rc.data() : query_bytes.data();
        bwa::align::Alignment sw_aln = bwa::align::sw_semi_global_extend(
            adaptive_scoring,
            std::span<const uint8_t>(qptr, query_len),
            std::span<const uint8_t>(ref_region.data(), ref_region.size()),
            effective_config.band_width
        );

        // Accept extensions that explain their chain: absolute floor plus
        // at least half the chain score (BWA drops extensions that keep
        // little of the seeding evidence; free query start/end clipping
        // would otherwise promote short spurious matches).
        if (sw_aln.score >= 10 && sw_aln.score * 2 >= chain.score) {
            all_alignments.push_back({ci, std::move(sw_aln)});
        }
    }

    if (all_alignments.empty()) {
        // Fall back to simple CIGAR from best chain (on its own reference)
        const auto& bmem = best_chain.mems[0];
        const auto& bref = index_.references()[static_cast<size_t>(best_chain.ref_id)];
        int32_t ref_pos = bmem.ref_pos;
        result.mapped = true;
        result.best_score = best_chain.score;
        result.second_best_score = chains.size() > 1 ? chains[1].score : 0;
        result.primary.qname = std::string(read.name.view());
        result.primary.rname = bref.name;
        result.primary.pos = std::max<int32_t>(1, ref_pos + 1);
        result.suboptimal_score = suboptimal_score;
        result.primary.seq = std::string(read.seq.view());
        result.primary.qual = std::string(read.qual.view());

        chain_to_alignment(best_chain, read, bref, result.primary);
        if (chain_is_reverse(best_chain)) {
            result.primary.flag |= AlnRecord::F_REVERSE;
            std::reverse(result.primary.cigar.begin(), result.primary.cigar.end());
        }
        // Aligned spans exclude soft clips (BWA qb/qe convention)
        int32_t fb_qspan = 0, fb_rspan = 0;
        for (uint32_t c : result.primary.cigar) {
            auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
            int len = bwa::align::cigar_len(c);
            if (op == bwa::align::CigarOp::Match || op == bwa::align::CigarOp::Equal ||
                op == bwa::align::CigarOp::Diff || op == bwa::align::CigarOp::Ins) {
                fb_qspan += len;
            }
            if (op == bwa::align::CigarOp::Match || op == bwa::align::CigarOp::Equal ||
                op == bwa::align::CigarOp::Diff || op == bwa::align::CigarOp::Del) {
                fb_rspan += len;
            }
        }
        MapqSignals fb_ms;
        fb_ms.score = best_chain.score;
        fb_ms.sub = result.second_best_score;
        fb_ms.csub = suboptimal_score;
        fb_ms.sub_n = chains.size() > 1 ? static_cast<int32_t>(chains.size()) - 1 : 0;
        fb_ms.seedcov = seed_query_coverage(best_chain);
        fb_ms.query_span = fb_qspan;
        fb_ms.ref_span = fb_rspan;
        fb_ms.match = adaptive_scoring.match;
        fb_ms.mismatch = -adaptive_scoring.mismatch;
        fb_ms.min_seed_len = mem_finder_.min_seed_len();
        result.primary.mapq = approx_mapq_se(fb_ms);
        return;
    }

    // Sort alignments by score (descending)
    std::sort(all_alignments.begin(), all_alignments.end(),
              [](const auto& a, const auto& b) { return a.second.score > b.second.score; });

    // Use the best alignment as primary (positions local to its reference,
    // while extract_ref needs global concatenated coordinates)
    const auto& [primary_chain_idx, primary_swaln] = all_alignments[0];
    const auto& pchain = chains[primary_chain_idx];
    int32_t p_ref_id = 0, p_ref_begin = 0, p_ref_end = 0;
    chain_ref_window(pchain, p_ref_id, p_ref_begin, p_ref_end);
    const auto& pref = index_.references()[static_cast<size_t>(p_ref_id)];
    int32_t p_global_base = static_cast<int32_t>(pref.offset);
    int32_t p_global_end = p_global_base + static_cast<int32_t>(pref.length);
    int32_t primary_local_region_begin = p_ref_begin - p_global_base;
    int32_t final_local_begin = primary_local_region_begin + primary_swaln.ref_begin;
    int32_t final_ref_begin = p_ref_begin + primary_swaln.ref_begin;
    int32_t final_ref_end = p_ref_begin + primary_swaln.ref_end;

    result.mapped = true;
    result.best_score = primary_swaln.score;
    result.second_best_score = all_alignments.size() > 1 ? all_alignments[1].second.score : 0;

    // Build primary alignment
    result.primary.qname = std::string(read.name.view());
    result.primary.rname = pref.name;
    result.primary.pos = final_local_begin + 1;
    MapqSignals ms;
    ms.score = primary_swaln.score;
    ms.sub = result.second_best_score;
    ms.csub = suboptimal_score;
    ms.sub_n = all_alignments.size() > 1 ? static_cast<int32_t>(all_alignments.size()) - 1 : 0;
    ms.seedcov = seed_query_coverage(pchain);
    ms.query_span = primary_swaln.query_end - primary_swaln.query_begin;
    ms.ref_span = primary_swaln.ref_end - primary_swaln.ref_begin;
    ms.match = adaptive_scoring.match;
    ms.mismatch = -adaptive_scoring.mismatch;
    ms.min_seed_len = mem_finder_.min_seed_len();
    result.primary.mapq = approx_mapq_se(ms);
    result.primary.score = primary_swaln.score;
    result.primary.seq = std::string(read.seq.view());
    result.primary.qual = std::string(read.qual.view());
    result.primary.cigar = primary_swaln.cigar;
    if (chain_rev[primary_chain_idx]) {
        result.primary.flag |= AlnRecord::F_REVERSE;
        std::reverse(result.primary.cigar.begin(), result.primary.cigar.end());
    }

    // Add NM and MD tags for primary
    int32_t nm = 0;
    for (uint32_t c : result.primary.cigar) {
        auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
        int len = bwa::align::cigar_len(c);
        if (op == bwa::align::CigarOp::Diff) nm += len;
        else if (op == bwa::align::CigarOp::Ins || op == bwa::align::CigarOp::Del) nm += len;
    }
    result.primary.tags.emplace_back("NM:i", std::to_string(nm));

    // MD tag (positions are global concatenated coordinates here;
    // the leading number is the match run length, per SAM spec)
    std::string md;
    int32_t ref_pos_in_aln = final_ref_begin;
    int32_t run_len = 0;
    bool first = true;
    for (uint32_t c : result.primary.cigar) {
        auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
        int len = bwa::align::cigar_len(c);
        if (op == bwa::align::CigarOp::Equal) {
            run_len += len;
            ref_pos_in_aln += len;
        } else if (op == bwa::align::CigarOp::Diff) {
            md += std::to_string(run_len);
            first = false;
            run_len = 0;
            for (int k = 0; k < len; ++k) {
                if (k > 0) md += "0";
                if (ref_pos_in_aln < p_global_end) {
                    auto ref_base = index_.extract_ref(ref_pos_in_aln, ref_pos_in_aln + 1);
                    md += !ref_base.empty() ? bwa::index::PackedSequence::decode_base(ref_base[0]) : 'N';
                } else { md += 'N'; }
                ref_pos_in_aln++;
            }
        } else if (op == bwa::align::CigarOp::Del) {
            md += std::to_string(run_len);
            first = false;
            run_len = 0;
            md += '^';
            for (int k = 0; k < len; ++k) {
                if (ref_pos_in_aln < p_global_end) {
                    auto ref_base = index_.extract_ref(ref_pos_in_aln, ref_pos_in_aln + 1);
                    md += !ref_base.empty() ? bwa::index::PackedSequence::decode_base(ref_base[0]) : 'N';
                } else { md += 'N'; }
                ref_pos_in_aln++;
            }
        }
    }
    if (first) { md = std::to_string(ref_pos_in_aln - final_ref_begin); }
    else { md += std::to_string(run_len); }
    result.primary.tags.emplace_back("MD:Z", md);

    // Gather alternate alignments (each on its own chain's reference).
    struct AltAln {
        size_t alt_idx;  // index into all_alignments
        int32_t chain_idx;
        int32_t local_pos;  // 1-based leftmost position on its reference
        bool rev;
    };
    std::vector<AltAln> alts;
    for (size_t i = 1; i < all_alignments.size() && i < 4; ++i) {
        const auto& [sec_chain_idx, sec_swaln] = all_alignments[i];
        const auto& schain = chains[sec_chain_idx];
        int32_t s_ref_id = 0, s_ref_begin = 0, s_ref_end = 0;
        if (!chain_ref_window(schain, s_ref_id, s_ref_begin, s_ref_end)) continue;
        const auto& sref = index_.references()[static_cast<size_t>(s_ref_id)];
        int32_t sec_local_begin = (s_ref_begin - static_cast<int32_t>(sref.offset)) +
                                  sec_swaln.ref_begin;
        alts.push_back({i, sec_chain_idx, sec_local_begin + 1,
                        chain_rev[sec_chain_idx] != 0});
    }

    // Split-read detection: an alternate whose aligned query span complements
    // (rather than overlaps) the primary span is a chimeric segment. It is
    // emitted as a hard-clipped supplementary via create_supplementary
    // instead of a secondary (BWA mem_reg2sam split logic).
    int split_alt = -1;
    {
        int32_t p_qb = primary_swaln.query_begin, p_qe = primary_swaln.query_end;
        for (size_t k = 0; k < alts.size(); ++k) {
            const auto& sw = all_alignments[alts[k].alt_idx].second;
            int32_t s_qb = sw.query_begin, s_qe = sw.query_end;
            int32_t overlap =
                std::max(0, std::min(p_qe, s_qe) - std::max(p_qb, s_qb));
            // The alternate must cover query territory the primary does not
            // (left or right extension). Without this, local-SW sub-matches
            // contained inside the primary span were misclassified as
            // chimeric segments, flooding output with supplementaries.
            int32_t ext = std::max(0, p_qb - s_qb) + std::max(0, s_qe - p_qe);
            int32_t combined = (p_qe - p_qb) + (s_qe - s_qb) - overlap;
            if (overlap <= std::max(5, query_len / 10) &&
                ext >= effective_config.min_seed_len &&
                combined >= query_len * 7 / 10 &&
                (s_qe - s_qb) >= effective_config.min_seed_len &&
                sw.score >= config_.min_chain_score) {
                split_alt = static_cast<int>(k);
                break;
            }
        }
    }

    for (size_t k = 0; k < alts.size(); ++k) {
        const auto& alt = alts[k];
        const auto& [sec_chain_idx, sec_swaln] = all_alignments[alt.alt_idx];
        const auto& schain = chains[sec_chain_idx];
        int32_t s_ref_id = 0, s_ref_begin = 0, s_ref_end = 0;
        if (!chain_ref_window(schain, s_ref_id, s_ref_begin, s_ref_end)) continue;
        const auto& sref = index_.references()[static_cast<size_t>(s_ref_id)];

        // Oriented CIGAR (reverse strand: mirrored op order).
        std::vector<uint32_t> oriented = sec_swaln.cigar;
        if (alt.rev) std::reverse(oriented.begin(), oriented.end());

        if (static_cast<int>(k) == split_alt) {
            // Strip end clips, then hard-clip the query outside [s_qb, s_qe).
            int32_t s_qb = sec_swaln.query_begin, s_qe = sec_swaln.query_end;
            size_t core_b = 0, core_e = oriented.size();
            auto is_clip = [](uint32_t c) {
                auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
                return op == bwa::align::CigarOp::SoftClip ||
                       op == bwa::align::CigarOp::HardClip;
            };
            while (core_b < core_e && is_clip(oriented[core_b])) ++core_b;
            while (core_e > core_b && is_clip(oriented[core_e - 1])) --core_e;
            if (core_b >= core_e) continue;  // fully clipped: no segment
            std::vector<uint32_t> hcigar;
            if (s_qb > 0) {
                hcigar.push_back(bwa::align::encode_cigar(
                    s_qb, bwa::align::CigarOp::HardClip));
            }
            hcigar.insert(hcigar.end(), oriented.begin() + core_b,
                          oriented.begin() + core_e);
            if (s_qe < query_len) {
                hcigar.push_back(bwa::align::encode_cigar(
                    query_len - s_qe, bwa::align::CigarOp::HardClip));
            }
            AlnRecord supp = create_supplementary(read, hcigar, sref.name,
                                                  alt.local_pos, alt.rev,
                                                  sec_swaln.score);
            // Hard-clipped bases are not stored in SEQ/QUAL (SAM spec).
            // Query spans are in RC coordinates for reverse segments.
            int32_t seg_b = alt.rev ? query_len - s_qe : s_qb;
            int32_t seg_len = s_qe - s_qb;
            std::string full_seq = std::string(read.seq.view());
            std::string full_qual = std::string(read.qual.view());
            supp.seq = full_seq.substr(static_cast<size_t>(seg_b),
                                       static_cast<size_t>(seg_len));
            supp.qual = full_qual.substr(static_cast<size_t>(seg_b),
                                         static_cast<size_t>(seg_len));
            result.supplementary.push_back(std::move(supp));
            continue;
        }

        AlnRecord sec;
        sec.qname = std::string(read.name.view());
        sec.flag = AlnRecord::F_SECONDARY;
        sec.rname = sref.name;
        sec.pos = alt.local_pos;
        // BWA assigns MAPQ 0 to secondary alignments (mem_reg2aln)
        sec.mapq = 0;
        sec.score = sec_swaln.score;
        sec.seq = std::string(read.seq.view());
        sec.qual = std::string(read.qual.view());
        sec.cigar = std::move(oriented);
        if (alt.rev) sec.flag |= AlnRecord::F_REVERSE;
        int32_t sec_nm = 0;
        for (uint32_t c : sec.cigar) {
            auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
            int len = bwa::align::cigar_len(c);
            if (op == bwa::align::CigarOp::Diff) sec_nm += len;
            else if (op == bwa::align::CigarOp::Ins || op == bwa::align::CigarOp::Del) sec_nm += len;
        }
        sec.tags.emplace_back("NM:i", std::to_string(sec_nm));
        result.secondary.push_back(std::move(sec));
    }

    // Promote cross-reference secondaries to supplementary + attach SA:Z tags.
    detect_supplementary(read, chains, arena, result);
}

// 5' reference position of an alignment (1-based POS is leftmost).
static int32_t five_prime_pos(const AlnRecord& a, int32_t ref_span) {
    bool rev = (a.flag & AlnRecord::F_REVERSE) != 0;
    return rev ? a.pos + ref_span - 1 : a.pos;
}

// Reference span (aligned ref length) from CIGAR, excluding clips.
static int32_t cigar_ref_span(const std::vector<uint32_t>& cigar) {
    int32_t span = 0;
    for (uint32_t c : cigar) {
        auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
        int len = bwa::align::cigar_len(c);
        if (op == bwa::align::CigarOp::Match || op == bwa::align::CigarOp::Equal ||
            op == bwa::align::CigarOp::Diff || op == bwa::align::CigarOp::Del) {
            span += len;
        }
    }
    return span;
}

bool Aligner::rescue_end(const io::SeqRecord& read, const AlnRecord& mate,
                         AlignmentResult& result) const {
    result.clear();
    int32_t qlen = static_cast<int32_t>(read.seq.size());
    if (qlen == 0) return false;

    // Locate mate in global concatenated coordinates.
    const RefSequence* mref = nullptr;
    for (size_t i = 0; i < index_.num_references(); ++i) {
        if (index_.references()[i].name == mate.rname) {
            mref = &index_.references()[i];
            break;
        }
    }
    if (mref == nullptr) return false;
    int32_t mate_global = static_cast<int32_t>(mref->offset) + (mate.pos - 1);

    // Search window around the mate (BWA infers this from the insert-size
    // distribution; we use a fixed window — documented deviation).
    constexpr int32_t kRescueWindow = 2000;
    int32_t w_begin = std::max<int32_t>(static_cast<int32_t>(mref->offset),
                                        mate_global - kRescueWindow);
    int32_t w_end = std::min<int32_t>(static_cast<int32_t>(mref->offset + mref->length),
                                      mate_global + kRescueWindow + qlen);
    if (w_end <= w_begin) return false;
    auto ref_region = index_.extract_ref(w_begin, w_end);
    if (ref_region.empty()) return false;

    // Pack query both strands.
    bwa::index::PackedSequence query;
    query.append(read.seq.data(), read.seq.size());
    auto qbases = query.bases();
    std::vector<uint8_t> qfwd(qlen), qrc(qlen);
    for (int32_t i = 0; i < qlen; ++i) {
        qfwd[i] = qbases[i];
        qrc[qlen - 1 - i] = rc_base(qbases[i]);
    }

    bwa::align::Alignment best;
    bool best_rev = false;
    bool have = false;
    for (int s = 0; s < 2; ++s) {
        const uint8_t* qptr = s == 0 ? qfwd.data() : qrc.data();
        auto aln = bwa::align::sw_semi_global_extend(
            config_.scoring,
            std::span<const uint8_t>(qptr, qlen),
            std::span<const uint8_t>(ref_region.data(), ref_region.size()),
            config_.band_width);
        if (!have || aln.score > best.score) {
            best = std::move(aln);
            best_rev = (s == 1);
            have = true;
        }
    }
    if (!have || best.score < config_.min_chain_score) return false;

    AlnRecord rec;
    rec.qname = std::string(read.name.view());
    if (best_rev) {
        rec.flag |= AlnRecord::F_REVERSE;
        std::reverse(best.cigar.begin(), best.cigar.end());
    }
    rec.rname = mref->name;
    rec.pos = (w_begin - static_cast<int32_t>(mref->offset)) + best.ref_begin + 1;
    rec.mapq = 0;  // re-evaluated by the pairing step
    rec.cigar = best.cigar;
    rec.seq = std::string(read.seq.view());
    rec.qual = std::string(read.qual.view());
    rec.score = best.score;
    rec.is_primary = true;
    int32_t nm = 0;
    for (uint32_t c : rec.cigar) {
        auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
        int len = bwa::align::cigar_len(c);
        if (op == bwa::align::CigarOp::Diff) nm += len;
        else if (op == bwa::align::CigarOp::Ins || op == bwa::align::CigarOp::Del) nm += len;
    }
    rec.tags.emplace_back("NM:i", std::to_string(nm));
    result.mapped = true;
    result.best_score = best.score;
    result.primary = std::move(rec);
    return true;
}

void Aligner::mark_pair(AlignmentResult& res, const io::SeqRecord& read,
                        bool is_read1, const AlnRecord* mate_primary,
                        bool mate_mapped, bool proper, int32_t tlen) const {
    uint32_t base = AlnRecord::F_PAIRED | (is_read1 ? AlnRecord::F_READ1 : AlnRecord::F_READ2);
    auto mark_alt = [&](AlnRecord& a) {
        a.flag |= base;
        if (proper) a.flag |= AlnRecord::F_PROPER_PAIR;
        // Non-primary records carry no mate coordinates (BWA convention).
        a.rnext = "*";
        a.pnext = 0;
        a.tlen = 0;
    };
    for (auto& s : res.secondary) mark_alt(s);
    for (auto& s : res.supplementary) mark_alt(s);

    AlnRecord& p = res.primary;
    if (!res.mapped) {
        p.qname = std::string(read.name.view());
        p.seq = std::string(read.seq.view());
        p.qual = std::string(read.qual.view());
        p.rname = "*";
        p.pos = 0;
        p.mapq = 0;
        p.cigar.clear();
        p.score = 0;
        p.is_primary = true;
        p.flag = base | AlnRecord::F_UNMAP;
        if (mate_mapped && mate_primary != nullptr) {
            p.rnext = mate_primary->rname;
            p.pnext = mate_primary->pos;
        } else {
            p.flag |= AlnRecord::F_MUNMAP;
            p.rnext = "*";
            p.pnext = 0;
        }
        p.tlen = 0;
        return;
    }
    // Preserve strand/secondary/supplementary bits; replace read-index bits.
    p.flag &= ~(AlnRecord::F_READ1 | AlnRecord::F_READ2);
    p.flag |= base;
    if (proper) p.flag |= AlnRecord::F_PROPER_PAIR;
    if (mate_mapped && mate_primary != nullptr) {
        p.rnext = (mate_primary->rname == p.rname) ? "=" : mate_primary->rname;
        p.pnext = mate_primary->pos;
        p.tlen = tlen;
    } else {
        p.flag |= AlnRecord::F_MUNMAP;
        p.rnext = "*";
        p.pnext = 0;
        p.tlen = 0;
    }
}

void Aligner::align_pair_impl(const bwa::io::SeqRecord& read1,
                              const bwa::io::SeqRecord& read2,
                              AlignmentResult& res1,
                              AlignmentResult& res2) const {
    align_impl(read1, res1);
    align_impl(read2, res2);

    // Mate rescue for unmapped ends (BWA mem_mate_rescue equivalent).
    bool rescued1 = false, rescued2 = false;
    if (!res1.mapped && res2.mapped) rescued1 = rescue_end(read1, res2.primary, res1);
    if (!res2.mapped && res1.mapped) rescued2 = rescue_end(read2, res1.primary, res2);

    // Proper-pair evaluation: same reference, FR orientation, insert size
    // within [0, max_gap]. BWA infers the insert distribution from the
    // library; we use a fixed bound — documented deviation.
    bool proper = false;
    int32_t tlen1 = 0, tlen2 = 0;
    if (res1.mapped && res2.mapped && res1.primary.rname == res2.primary.rname) {
        int32_t span1 = cigar_ref_span(res1.primary.cigar);
        int32_t span2 = cigar_ref_span(res2.primary.cigar);
        int32_t five1 = five_prime_pos(res1.primary, span1);
        int32_t five2 = five_prime_pos(res2.primary, span2);
        bool fwd1 = (res1.primary.flag & AlnRecord::F_REVERSE) == 0;
        bool fwd2 = (res2.primary.flag & AlnRecord::F_REVERSE) == 0;
        bool fr = (fwd1 && !fwd2 && five1 <= five2) ||
                  (!fwd1 && fwd2 && five2 <= five1);
        int32_t insert = (five1 <= five2) ? (five2 - five1) : (five1 - five2);
        if (fr && insert <= config_.max_gap) {
            proper = true;
            // TLEN sign: positive when this end holds the leftmost 5'.
            tlen1 = (five1 <= five2) ? insert : -insert;
            tlen2 = -tlen1;
        }
    }

    // Rescued ends that form a proper pair earn a modest MAPQ; BWA derives
    // this from the pairing score, we use a fixed credit — documented.
    constexpr uint8_t kRescueMapq = 20;
    if (proper) {
        if (rescued1) res1.primary.mapq = std::min<uint8_t>(kRescueMapq, res2.primary.mapq);
        if (rescued2) res2.primary.mapq = std::min<uint8_t>(kRescueMapq, res1.primary.mapq);
    }

    const AlnRecord* mate_of_1 = res2.mapped ? &res2.primary : nullptr;
    const AlnRecord* mate_of_2 = res1.mapped ? &res1.primary : nullptr;
    mark_pair(res1, read1, true, mate_of_1, res2.mapped, proper, tlen1);
    mark_pair(res2, read2, false, mate_of_2, res1.mapped, proper, tlen2);
}

void Aligner::chain_to_alignment(const bwa::align::MEMFinder::Chain& chain,
                                 const bwa::io::SeqRecord& read,
                                 const RefSequence& ref,
                                 AlnRecord& aln) const {
    aln.cigar.clear();
    int32_t prev_query_end = 0;
    int32_t query_len = static_cast<int32_t>(read.seq.size());

    for (size_t i = 0; i < chain.mems.size(); ++i) {
        const auto& mem = chain.mems[i];

        // Add soft clip for gap before this MEM
        if (mem.query_pos > prev_query_end) {
            int32_t clip_len = mem.query_pos - prev_query_end;
            aln.cigar.push_back(
                bwa::align::encode_cigar(clip_len, bwa::align::CigarOp::SoftClip));
        }

        // Add the MEM as a match
        aln.cigar.push_back(
            bwa::align::encode_cigar(mem.len, bwa::align::CigarOp::Match));

        prev_query_end = mem.query_end();
    }

    // Add trailing soft clip
    if (prev_query_end < query_len) {
        int32_t clip_len = query_len - prev_query_end;
        aln.cigar.push_back(
            bwa::align::encode_cigar(clip_len, bwa::align::CigarOp::SoftClip));
    }

    aln.score = chain.score;
}

uint8_t Aligner::approx_mapq_se(const MapqSignals& s) noexcept {
    // Port of BWA mem_approx_mapq_se() (bwamem.c). Constants from
    // mem_opt_init(): MEM_MAPQ_COEF 30.0, MEM_MAPQ_MAX 60, mapQ_coef_len 50,
    // mapQ_coef_fac log(50); 6.02 = 10*log(4)/log(10), 4.343 = 10/log(10).
    int sub = (s.sub != 0) ? s.sub : s.min_seed_len * s.match;
    if (s.csub > sub) sub = s.csub;
    if (sub >= s.score) return 0;
    int32_t alen = s.query_span > s.ref_span ? s.query_span : s.ref_span;
    if (alen <= 0 || s.match <= 0 || s.mismatch < 0) return 0;  // hardening
    double len = static_cast<double>(alen);
    double match = static_cast<double>(s.match);
    double mismatch = static_cast<double>(s.mismatch);
    double identity = 1.0 - (len * match - s.score) / (match + mismatch) / len;
    if (s.score == 0) {
        return 0;
    }
    int mapq;
    constexpr double coef_len = 50.0;
    constexpr double coef_fac = 3.912023005428146;  // log(50)
    if (coef_len > 0) {
        double tmp = len < coef_len ? 1.0 : coef_fac / std::log(len);
        tmp *= identity * identity;
        mapq = static_cast<int>(6.02 * (s.score - sub) / match * tmp * tmp + .499);
    } else {
        if (s.seedcov < 1) return 0;
        mapq = static_cast<int>(30.0 * (1.0 - static_cast<double>(sub) / s.score) *
                                std::log(static_cast<double>(s.seedcov)) +
                                .499);
        if (identity < 0.95) {
            mapq = static_cast<int>(mapq * identity * identity + .499);
        }
    }
    if (s.sub_n > 0) {
        mapq -= static_cast<int>(4.343 * std::log(static_cast<double>(s.sub_n) + 1.0) + .499);
    }
    if (mapq > 60) mapq = 60;
    if (mapq < 0) mapq = 0;
    return static_cast<uint8_t>(mapq);
}

int32_t Aligner::seed_query_coverage(const align::MEMFinder::Chain& chain) {
    // Union length of member MEM query spans (members are query-ordered)
    int32_t total = 0;
    int32_t cur_end = -1;
    for (size_t i = 0; i < chain.mems.size(); ++i) {
        int32_t b = chain.mems[i].query_pos;
        int32_t e = chain.mems[i].query_end();
        if (b > cur_end) {
            total += e - b;
            cur_end = e;
        } else if (e > cur_end) {
            total += e - cur_end;
            cur_end = e;
        }
    }
    return total;
}



// Detect and create supplementary alignments for chimeric reads.
// Any secondary alignment on a different reference than the primary is a
// chimeric candidate: move it from secondary to supplementary, flag it
// F_SUPPLEMENTARY (clearing F_SECONDARY), and attach SA:Z tags linking
// primary <-> supplementaries (BWA mem_mark_secondary / mem_reg2sam logic).
void Aligner::detect_supplementary(const io::SeqRecord& read,
                                   const bwa::core::Vector<bwa::align::MEMFinder::Chain>& chains,
                                   const bwa::memory::Arena& arena,
                                   AlignmentResult& result) const {
    (void)read;
    (void)chains;
    (void)arena;
    if (!config_.output_supplementary) return;
    if (!result.mapped) return;
    if (result.secondary.empty() && result.supplementary.empty()) return;

    // Move secondaries on a different reference than primary to supplementary.
    std::vector<AlnRecord> kept_secondary;
    kept_secondary.reserve(result.secondary.size());
    for (auto& sec : result.secondary) {
        if (sec.flag & AlnRecord::F_UNMAP) {
            kept_secondary.push_back(std::move(sec));
            continue;
        }
        if (sec.rname != result.primary.rname &&
            sec.score >= config_.min_chain_score &&
            result.supplementary.size() < 4) {
            AlnRecord supp = std::move(sec);
            supp.is_supplementary = true;
            supp.is_primary = false;
            supp.flag |= AlnRecord::F_SUPPLEMENTARY;
            supp.flag &= ~AlnRecord::F_SECONDARY;  // supplementary, not secondary
            supp.mapq = 0;  // BWA assigns MAPQ 0 to supplementary
            result.supplementary.push_back(std::move(supp));
        } else {
            kept_secondary.push_back(std::move(sec));
        }
    }
    result.secondary = std::move(kept_secondary);
    if (result.supplementary.empty()) return;

    // Build SA:Z tags once: primary lists all supplementaries, each
    // supplementary lists primary + other supplementaries. Format per entry:
    // RNAME,POS,STRAND,CIGAR,MAPQ,NM;
    for (const auto& t : result.primary.tags) {
        if (t.first == "SA:Z") return;  // already tagged
    }
    auto entry = [this](const AlnRecord& a) { return make_sa_tag(a); };
    std::string primary_sa;
    for (const auto& supp : result.supplementary) primary_sa += entry(supp);
    if (!primary_sa.empty()) {
        result.primary.tags.emplace_back("SA:Z", primary_sa);
        for (size_t i = 0; i < result.supplementary.size(); ++i) {
            std::string sa = entry(result.primary);
            for (size_t j = 0; j < result.supplementary.size(); ++j) {
                if (j == i) continue;
                sa += entry(result.supplementary[j]);
            }
            result.supplementary[i].tags.emplace_back("SA:Z", sa);
        }
    }
}

// Create supplementary alignment for a split-read segment.
AlnRecord Aligner::create_supplementary(const io::SeqRecord& read,
                                        const std::vector<uint32_t>& cigar,
                                        const std::string& rname,
                                        int32_t pos_1based,
                                        bool is_reverse,
                                        int32_t score) const {
    AlnRecord supp;
    supp.qname = std::string(read.name.view());
    supp.flag = AlnRecord::F_SUPPLEMENTARY;
    if (is_reverse) supp.flag |= AlnRecord::F_REVERSE;
    supp.rname = rname;
    supp.pos = pos_1based;
    supp.mapq = 0;  // BWA assigns MAPQ 0 to supplementary alignments
    supp.cigar = cigar;
    supp.seq = std::string(read.seq.view());
    supp.qual = std::string(read.qual.view());
    supp.score = score;
    supp.is_supplementary = true;
    supp.is_primary = false;

    // NM over aligned ops only (H/S clips excluded by op filter).
    int32_t nm = 0;
    for (uint32_t c : supp.cigar) {
        auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
        int len = bwa::align::cigar_len(c);
        if (op == bwa::align::CigarOp::Diff) nm += len;
        else if (op == bwa::align::CigarOp::Ins || op == bwa::align::CigarOp::Del) nm += len;
    }
    supp.tags.emplace_back("NM:i", std::to_string(nm));

    return supp;
}

// Compute single SA:Z entry for an alignment: RNAME,POS,STRAND,CIGAR,MAPQ,NM;
std::string Aligner::make_sa_tag(const AlnRecord& aln) const {
    std::string sa;
    sa += aln.rname + ',';
    sa += std::to_string(aln.pos) + ',';
    sa += ((aln.flag & AlnRecord::F_REVERSE) ? "-" : "+");
    sa += ',';
    sa += format_cigar(aln.cigar) + ',';
    sa += std::to_string(static_cast<int>(aln.mapq)) + ',';
    // NM value: tags use typed keys ("NM:i"); accept legacy "NM" too.
    for (const auto& tag : aln.tags) {
        if (tag.first == "NM:i" || tag.first == "NM") {
            sa += tag.second;
            break;
        }
    }
    sa += ';';
    return sa;
}

void Pipeline::write_header(std::ostream& out) const {
    out << "@HD\tVN:1.6\tSO:coordinate\n";

    for (size_t i = 0; i < index_.num_references(); ++i) {
        const auto& ref = index_.get_ref(i);
        out << "@SQ\tSN:" << ref.name
            << "\tLN:" << ref.length;
        if (!ref.md5.empty()) {
            out << "\tM5:" << ref.md5;
        }
        out << "\n";
    }

    // @RG header (if read group is specified)
    if (config_.read_group.has_value()) {
        const auto& rg = config_.read_group.value();
        out << "@RG\tID:" << rg.id;
        if (!rg.sample.empty()) out << "\tSM:" << rg.sample;
        if (!rg.library.empty()) out << "\tLB:" << rg.library;
        if (!rg.platform.empty()) out << "\tPL:" << rg.platform;
        if (!rg.platform_unit.empty()) out << "\tPU:" << rg.platform_unit;
        out << "\n";
    }

    // @PG header (program info)
    out << "@PG\tID:" << config_.program_name
        << "\tPN:" << config_.program_name
        << "\tVN:" << config_.program_version
        << "\tCL:" << config_.program_command << "\n";
}

void Pipeline::write_alignment(std::ostream& out, const AlignmentResult& result) const {
    if (!result.mapped) {
        return;
    }

    // Write primary alignment with suboptimal score
    AlnRecord primary = result.primary;
    primary.suboptimal_score = result.suboptimal_score;
    write_sam_record(out, primary);

    // Write secondary alignments (BWA suppresses these unless -a)
    if (config_.output_secondary) {
        for (const auto& sec : result.secondary) {
            write_sam_record(out, sec);
        }
    }

    // Write supplementary alignments (already carry SA:Z tags)
    if (config_.output_supplementary) {
        for (const auto& supp : result.supplementary) {
            write_sam_record(out, supp);
        }
    }
}

// Paired-end output: both ends adjacent (read1 records, then read2
// records), primaries emitted even when unmapped per SAM.
void Pipeline::write_pair(std::ostream& out, const AlignmentResult& res1,
                          const AlignmentResult& res2) const {
    auto write_end = [&](const AlignmentResult& res) {
        AlnRecord primary = res.primary;
        if (res.mapped) primary.suboptimal_score = res.suboptimal_score;
        write_sam_record(out, primary);
        if (!res.mapped) return;
        if (config_.output_secondary) {
            for (const auto& sec : res.secondary) write_sam_record(out, sec);
        }
        if (config_.output_supplementary) {
            for (const auto& supp : res.supplementary) write_sam_record(out, supp);
        }
    };
    write_end(res1);
    write_end(res2);
}

void Pipeline::write_sam_record(std::ostream& out, const AlnRecord& aln) const {
    out << aln.qname << '\t';
    out << aln.flag << '\t';
    out << aln.rname << '\t';
    out << aln.pos << '\t';
    out << static_cast<int>(aln.mapq) << '\t';
    if (aln.cigar.empty()) out << "*\t";
    else out << format_cigar(aln.cigar) << '\t';
    out << aln.rnext << '\t'
        << aln.pnext << '\t'
        << aln.tlen << '\t';
    out << aln.seq << '\t';
    out << aln.qual;

    // Optional tags
    for (const auto& tag : aln.tags) {
        out << '\t' << tag.first << ':' << tag.second;
    }
    // XS tag: suboptimal alignment score
    if (aln.suboptimal_score > 0) {
        out << "\tXS:i:" << aln.suboptimal_score;
    }
    out << '\n';
}

// Simple MD5 implementation (sufficient for reference checksums)
struct MD5 {
    uint32_t state[4] = {0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476};
    uint64_t count = 0;
    std::array<uint8_t, 64> buffer{};
    size_t buffer_len = 0;

    static constexpr uint32_t S[] = {
        7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
        5,  9, 14, 20, 5,  9, 14, 20, 5,  9, 14, 20, 5,  9, 14, 20,
        4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
        6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21
    };

    static constexpr uint32_t K[] = {
        0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee,
        0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
        0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
        0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
        0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa,
        0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
        0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed,
        0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
        0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
        0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
        0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05,
        0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
        0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039,
        0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
        0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
        0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391
    };

    void transform(const uint8_t* block) {
        uint32_t a = state[0], b = state[1], c = state[2], d = state[3];
        uint32_t M[16];
        std::memcpy(M, block, 64);

        for (int i = 0; i < 64; ++i) {
            uint32_t f, g;
            if (i < 16) { f = (b & c) | (~b & d); g = i; }
            else if (i < 32) { f = (d & b) | (~d & c); g = (5 * i + 1) % 16; }
            else if (i < 48) { f = b ^ c ^ d; g = (3 * i + 5) % 16; }
            else { f = c ^ (b | ~d); g = (7 * i) % 16; }
            uint32_t temp = d;
            d = c; c = b;
            b += ((a + f + K[i] + M[g]) << S[i]) | (a + f + K[i] + M[g]) >> (32 - S[i]);
            a = temp;
        }
        state[0] += a; state[1] += b; state[2] += c; state[3] += d;
    }

    void update(const uint8_t* data, size_t len) {
        count += len;
        while (len > 0) {
            size_t to_copy = std::min(len, size_t(64) - buffer_len);
            std::memcpy(buffer.data() + buffer_len, data, to_copy);
            buffer_len += to_copy;
            data += to_copy;
            len -= to_copy;
            if (buffer_len == 64) {
                transform(buffer.data());
                buffer_len = 0;
            }
        }
    }

    std::array<uint8_t, 16> finalize() {
        uint64_t bit_count = count * 8;
        update(reinterpret_cast<const uint8_t*>("\x80"), 1);
        while (buffer_len != 56) {
            update(reinterpret_cast<const uint8_t*>("\0"), 1);
        }
        for (int i = 0; i < 8; ++i) {
            update(reinterpret_cast<const uint8_t*>(&bit_count) + i, 1);
        }
        std::array<uint8_t, 16> digest;
        for (int i = 0; i < 4; ++i) {
            digest[i*4] = state[i] & 0xFF;
            digest[i*4+1] = (state[i] >> 8) & 0xFF;
            digest[i*4+2] = (state[i] >> 16) & 0xFF;
            digest[i*4+3] = (state[i] >> 24) & 0xFF;
        }
        return digest;
    }
};

std::string compute_md5(const std::string& data) {
    MD5 md5;
    md5.update(reinterpret_cast<const uint8_t*>(data.data()), data.size());
    auto digest = md5.finalize();
    std::ostringstream oss;
    for (auto b : digest) {
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b);
    }
    return oss.str();
}

void Index::build_impl(const char* fasta_path, const Config& cfg) {
    bwa::io::SeqReader reader(fasta_path);
    if (!reader.is_open()) {
        throw std::runtime_error("Failed to open FASTA file");
    }

    bwa::io::SeqRecord rec;
    size_t total_len = 0;

    auto read_result = reader.read(rec);
    while (read_result && *read_result) {
        if (rec.is_fasta() && !rec.seq.empty()) {
            RefSequence ref;
            ref.name = std::string(rec.name.view());
            ref.length = rec.seq.size();
            ref.offset = total_len;

            // Compute MD5 of the sequence (uppercase, no newlines)
            {
                std::string upper_seq;
                upper_seq.reserve(rec.seq.size());
                for (size_t i = 0; i < rec.seq.size(); ++i) {
                    char c = rec.seq.view()[i];
                    upper_seq += (c >= 'a' && c <= 'z') ? (c - 32) : c;
                }
                ref.md5 = compute_md5(upper_seq);
            }

            // Pack sequence
            bwa::index::PackedSequence seq;
            seq.append(rec.seq.data(), rec.seq.size());

            // Store packed reference for later extraction
            packed_refs_.push_back(seq);

            fm_index_.add_sequence(seq, bwa::memory::get_tls_arena());

            refs_.push_back(std::move(ref));
            total_len += rec.seq.size();
        }
        rec.clear();
        read_result = reader.read(rec);
    }

    // Build metadata string (v3: 1-byte BWT codes with sentinel row,
    // 6-symbol occ tables, N-mask in .pac)
    std::ostringstream oss;
    oss << "BWA-CPP26-METADATA-v3\n";
    oss << refs_.size() << "\n";
    oss << total_len << "\n";
    for (const auto& ref : refs_) {
        oss << ref.name << "\t" << ref.length << "\t" << ref.offset << "\t" << ref.md5 << "\n";
    }
    meta_ = oss.str();
}

void Index::save_impl(const char* prefix) const {
    // Save metadata
    std::string meta_path = std::string(prefix) + ".meta";
    std::ofstream meta_out(meta_path, std::ios::binary);
    if (!meta_out) {
        throw std::runtime_error("Cannot open metadata file for writing");
    }
    if (!meta_.empty()) {
        meta_out.write(meta_.data(), static_cast<std::streamsize>(meta_.size()));
    }
    meta_out.close();

    // Save BWT: v3 format stores one byte per row (sort codes $=0..N=5),
    // length+1 rows per reference (sentinel row included). N is a distinct
    // code, so no mask section is needed. (v2 stored 2-bit words + mask and
    // v1 dropped the mask entirely.)
    std::string bwt_path = std::string(prefix) + ".bwt";
    std::ofstream bwt_out(bwt_path, std::ios::binary);
    for (const auto& fm : fm_index_) {
        const auto& bwt_vec = fm.bwt();
        bwt_out.write(reinterpret_cast<const char*>(bwt_vec.data()),
                      static_cast<std::streamsize>(bwt_vec.size()));
    }
    bwt_out.close();

    // Save SA samples
    std::string sa_path = std::string(prefix) + ".sa";
    std::ofstream sa_out(sa_path, std::ios::binary);
    for (const auto& fm : fm_index_) {
        sa_out.write(reinterpret_cast<const char*>(fm.sa_samples().data()),
                      static_cast<std::streamsize>(fm.sa_samples().size() * sizeof(uint32_t)));
    }
    sa_out.close();

    // Save occ
    std::string occ_path = std::string(prefix) + ".occ";
    std::ofstream occ_out(occ_path, std::ios::binary);
    for (const auto& fm : fm_index_) {
        occ_out.write(reinterpret_cast<const char*>(fm.occ_table().data()),
                      static_cast<std::streamsize>(fm.occ_table().size() * sizeof(uint32_t)));
    }
    occ_out.close();

    // Save packed reference: data words + N-mask words per reference.
    // (The mask restore keeps N vs ACGT exact in extract_ref for SW/MD.)
    std::string pac_path = std::string(prefix) + ".pac";
    std::ofstream pac_out(pac_path, std::ios::binary);
    if (!pac_out) {
        throw std::runtime_error("Cannot open packed reference file for writing");
    }
    for (const auto& packed : packed_refs_) {
        pac_out.write(reinterpret_cast<const char*>(packed.words().data()),
                      static_cast<std::streamsize>(packed.words().size() * sizeof(uint64_t)));
        pac_out.write(reinterpret_cast<const char*>(packed.n_masks().data()),
                      static_cast<std::streamsize>(packed.n_masks().size() * sizeof(uint64_t)));
    }
    pac_out.close();
}

void Index::load_impl(const char* prefix) {
    // Load metadata
    std::string meta_path = std::string(prefix) + ".meta";
    std::ifstream meta_in(meta_path, std::ios::binary);
    if (!meta_in) throw std::runtime_error("Cannot open index metadata");

    meta_in.seekg(0, std::ios::end);
    size_t meta_size = meta_in.tellg();
    meta_in.seekg(0);
    meta_.resize(meta_size);
    meta_in.read(meta_.data(), meta_size);
    meta_in.close();

    // Parse metadata (v3 required: 1-byte BWT codes with sentinel row,
    // 6-symbol occ tables, N-mask in .pac; older files misload silently,
    // so refuse them explicitly)
    std::string_view mv = meta_;
    if (!mv.starts_with("BWA-CPP26-METADATA-v3\n")) {
        throw std::runtime_error(
            "Unsupported index version (need v3, found incompatible metadata). "
            "Rebuild the index with this bwa-cpp26 version.");
    }
    size_t pos = mv.find('\n') + 1;
    size_t num_refs = 0, total_len = 0;
    sscanf(mv.data() + pos, "%zu\n%zu\n", &num_refs, &total_len);

    // After sscanf, pos is at the start of the first ref line
    // (sscanf doesn't update pos, so we need to find it)
    pos = mv.find('\n', pos) + 1;
    pos = mv.find('\n', pos) + 1;

    refs_.reserve(num_refs);
    for (size_t i = 0; i < num_refs; ++i) {
        RefSequence ref;
        size_t end = mv.find('\t', pos);
        ref.name = std::string(mv.substr(pos, end - pos));
        pos = end + 1;
        end = mv.find('\t', pos);
        ref.length = std::stoull(std::string(mv.substr(pos, end - pos)));
        pos = end + 1;
        end = mv.find('\t', pos);
        ref.offset = std::stoull(std::string(mv.substr(pos, end - pos)));
        pos = end + 1;
        end = mv.find('\n', pos);
        if (end != std::string_view::npos) {
            ref.md5 = std::string(mv.substr(pos, end - pos));
        }
        pos = end + 1;
        refs_.push_back(std::move(ref));
    }

    // Load BWT using memory-mapped file for large genome support
    std::string bwt_path = std::string(prefix) + ".bwt";
    bwa::io::MmapFile bwt_mmap;
    bwt_mmap.open(bwt_path.c_str());
    if (!bwt_mmap.is_open()) throw std::runtime_error("Cannot open BWT file");

    size_t bwt_size = bwt_mmap.size();
    const uint8_t* bwt_data = static_cast<const uint8_t*>(bwt_mmap.data());
    size_t bwt_bytes = bwt_size;

    // Reconstruct one FM-index per reference, keeping fm_index_[ri] aligned
    // with refs_[ri]. save_impl writes each index's BWT bytes (length+1 rows
    // of sort codes), SA and occ arrays consecutively, so split the blobs
    // using per-ref sizes from metadata.
    // (A single combined index would be invalid: BWT(concat) != concat(BWTs).)
    {
        size_t bwt_byte_off = 0;  // in bytes
        for (size_t ri = 0; ri < num_refs; ++ri) {
            const size_t len = refs_[ri].length;
            const size_t rows = len + 1;

            bwa::index::FMIndex idx;
            idx.set_size(len);

            // BWT slice: rows bytes of sort codes ($=0..N=5)
            std::vector<uint8_t> bwt_vec(rows, 0);
            for (size_t i = 0; i < rows; ++i) {
                if (bwt_byte_off < bwt_bytes) {
                    bwt_vec[i] = bwt_data[bwt_byte_off];
                }
                bwt_byte_off++;
            }
            idx.set_bwt(bwt_vec);

            fm_index_.add_index(std::move(idx));
        }
    }

    // Load SA samples (one slice per reference)
    {
        std::string sa_path = std::string(prefix) + ".sa";
        bwa::io::MmapFile sa_mmap;
        sa_mmap.open(sa_path.c_str());
        if (!sa_mmap.is_open()) throw std::runtime_error("Cannot open SA file");

        const uint32_t* sa_data = static_cast<const uint32_t*>(sa_mmap.data());
        size_t sa_avail = sa_mmap.size() / sizeof(uint32_t);
        size_t off = 0;
        for (size_t ri = 0; ri < num_refs; ++ri) {
            // SA has length+1 rows (sentinel included)
            size_t want = (refs_[ri].length + 1 + 31) / 32;
            size_t have = (off < sa_avail) ? std::min(want, sa_avail - off) : 0;
            fm_index_[ri].set_sa_samples(
                std::vector<uint32_t>(sa_data + off, sa_data + off + have));
            off += want;
        }
    }

    // Load occ tables + rebuild count tables per reference
    {
        std::string occ_path = std::string(prefix) + ".occ";
        bwa::io::MmapFile occ_mmap;
        occ_mmap.open(occ_path.c_str());
        if (!occ_mmap.is_open()) throw std::runtime_error("Cannot open OCC file");

        const uint32_t* occ_data = static_cast<const uint32_t*>(occ_mmap.data());
        size_t occ_avail = occ_mmap.size() / sizeof(uint32_t);
        size_t off = 0;
        for (size_t ri = 0; ri < num_refs; ++ri) {
            // Rows = length+1 (sentinel row included)
            size_t nrows = refs_[ri].length + 1;
            size_t num_intv = (nrows + 128 - 1) / 128 + 1;
            size_t want = num_intv * bwa::index::FMIndex::ALPHABET;
            size_t have = (off < occ_avail) ? std::min(want, occ_avail - off) : 0;
            std::vector<uint32_t> occ(occ_data + off, occ_data + off + have);
            if (occ.size() < want) occ.resize(want, 0);
            off += want;
            fm_index_[ri].set_occ_table(occ);

            // Rebuild count table from this ref's final interval
            std::vector<uint32_t> cnt(bwa::index::FMIndex::ALPHABET + 1, 0);
            size_t last_intv = (nrows + 128 - 1) / 128;
            for (int b = 0; b < bwa::index::FMIndex::ALPHABET; ++b) {
                cnt[static_cast<size_t>(b) + 1] =
                    cnt[static_cast<size_t>(b)] +
                    occ[last_intv * bwa::index::FMIndex::ALPHABET + static_cast<size_t>(b)];
            }
            fm_index_[ri].set_count_table(cnt);
        }
    }

    // Load packed reference
    std::string pac_path = std::string(prefix) + ".pac";
    bwa::io::MmapFile pac_mmap;
    pac_mmap.open(pac_path.c_str());
    if (!pac_mmap.is_open()) {
        // No packed reference - extract_ref will return empty
        return;
    }
    const uint64_t* pac_data = static_cast<const uint64_t*>(pac_mmap.data());
    size_t pac_words = pac_mmap.size() / sizeof(uint64_t);

    // Distribute packed data among references. save_impl writes each ref's
    // whole 2-bit words followed by its N-mask words, so advance by whole
    // words per ref (a flat base cursor would drift when a length isn't a
    // multiple of 32). Restoring the mask keeps N exact in SW/MD.
    packed_refs_.resize(refs_.size());
    size_t word_off = 0;
    for (size_t i = 0; i < refs_.size(); ++i) {
        packed_refs_[i].resize(refs_[i].length);
        size_t data_words = (refs_[i].length + 31) / 32;
        for (size_t j = 0; j < refs_[i].length; ++j) {
            size_t wi = word_off + j / 32;
            size_t bi = (j % 32) * 2;
            uint8_t val = 0;
            if (wi < pac_words) {
                val = static_cast<uint8_t>((pac_data[wi] >> bi) & 0x3);
            }
            packed_refs_[i].set(j, val);
            size_t mi = word_off + data_words + j / 64;
            if (mi < pac_words) {
                packed_refs_[i].set_n(j, ((pac_data[mi] >> (j % 64)) & 1) != 0);
            }
        }
        word_off += data_words + (refs_[i].length + 63) / 64;
    }
}

} // namespace bwa
