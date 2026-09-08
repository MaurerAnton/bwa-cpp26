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
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>

namespace bwa {

namespace {

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

// Helper: extract a subregion of the packed reference as uint8_t span
std::vector<uint8_t> extract_ref_region(const bwa::index::FMIndex& fm,
                                        size_t start, size_t end) {
    std::vector<uint8_t> result;
    if (end > fm.length()) end = fm.length();
    if (start >= end) return result;
    result.reserve(end - start);
    for (size_t i = start; i < end; ++i) {
        // The BWT[i] gives us the character. To get the actual reference
        // sequence, we need to use the SA. For simplicity, we'll use
        // backward search to extract a substring.
        // Actually, for extension, we need the original reference text.
        // The FMIndex doesn't directly store the original reference.
        // We can reconstruct it by repeatedly applying LF from the primary.
        // For now, just get the BWT character.
        result.push_back(static_cast<uint8_t>(fm.bwt().get(i)));
    }
    return result;
}

// Extract the actual reference sequence at [ref_start, ref_end)
// by walking through the SA from a sampled position
std::vector<uint8_t> extract_ref_segment(size_t ref_start, size_t ref_len) {
    // Placeholder - in a real implementation, we'd use the packed reference
    // stored alongside the FM-index. For now, we can only return N's.
    // The actual implementation requires saving the original reference
    // alongside the index.
    return std::vector<uint8_t>(ref_len, 0); // A's as placeholder
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

    // Find MEMs (use TLS arena for memory)
    bwa::memory::Arena& arena = bwa::memory::get_tls_arena();
    auto mems = mem_finder_.find(query.bases(), arena);
    mem_finder_.filter_overlaps(mems);

    // Rescue scan if no MEMs found
    if (mems.empty()) {
        mem_finder_.rescue_scan(query.bases(), mems, true);
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

    // Get reference sequence info
    if (index_.num_references() == 0) {
        result.mapped = false;
        return;
    }

    const auto& ref_seq = index_.references()[0];
    int32_t ref_len = static_cast<int32_t>(ref_seq.length);

    // Pack query for SW
    std::vector<uint8_t> query_bytes(query_len);
    auto qbases = query.bases();
    for (int32_t i = 0; i < query_len; ++i) {
        query_bytes[i] = qbases[i];
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

    for (int32_t ci = 0; ci < max_alignments && ci < (int32_t)chains.size(); ++ci) {
        const auto& chain = chains[ci];
        int32_t chain_ref_begin = chain.mems.front().ref_pos;

        int32_t padding = effective_config.band_width * 2;
        int32_t ref_begin = std::max<int32_t>(0, chain_ref_begin - padding);
        int32_t ref_end = std::min<int32_t>(ref_len, chain_ref_begin + query_len + padding);

        if (ref_end <= ref_begin) continue;

        auto ref_region = index_.extract_ref(ref_begin, ref_end);
        if (ref_region.empty()) continue;

        bwa::align::Alignment sw_aln = bwa::align::sw_semi_global_extend(
            adaptive_scoring,
            std::span<const uint8_t>(query_bytes.data(), query_len),
            std::span<const uint8_t>(ref_region.data(), ref_region.size()),
            effective_config.band_width
        );

        if (sw_aln.score > 0) {
            all_alignments.push_back({ci, std::move(sw_aln)});
        }
    }

    if (all_alignments.empty()) {
        // Fall back to simple CIGAR from best chain
        int32_t ref_pos = best_chain.mems[0].ref_pos;
        result.mapped = true;
        result.best_score = best_chain.score;
        result.best_score = best_chain.score;
        result.second_best_score = chains.size() > 1 ? chains[1].score : 0;
        result.primary.qname = std::string(read.name.view());
        result.primary.rname = ref_seq.name;
        result.primary.pos = std::max<int32_t>(1, ref_pos + 1);
        result.primary.mapq = compute_mapq(result.best_score, result.second_best_score);
        result.primary.score = best_chain.score;
        result.suboptimal_score = suboptimal_score;
        result.primary.seq = std::string(read.seq.view());
        result.primary.qual = std::string(read.qual.view());

        int32_t prev_query_end = 0;
        for (const auto& mem : best_chain.mems) {
            if (mem.query_pos > prev_query_end) {
                result.primary.cigar.push_back(
                    bwa::align::encode_cigar(mem.query_pos - prev_query_end, bwa::align::CigarOp::SoftClip));
            }
            result.primary.cigar.push_back(
                bwa::align::encode_cigar(mem.len, bwa::align::CigarOp::Match));
            prev_query_end = mem.query_end();
        }
        if (prev_query_end < query_len) {
            result.primary.cigar.push_back(
                bwa::align::encode_cigar(query_len - prev_query_end, bwa::align::CigarOp::SoftClip));
        }
        return;
    }

    // Sort alignments by score (descending)
    std::sort(all_alignments.begin(), all_alignments.end(),
              [](const auto& a, const auto& b) { return a.second.score > b.second.score; });

    // Use the best alignment as primary
    const auto& [primary_chain_idx, primary_swaln] = all_alignments[0];
    int32_t primary_ref_begin = chains[primary_chain_idx].mems.front().ref_pos;
    int32_t padding = effective_config.band_width * 2;
    int32_t primary_ref_region_begin = std::max<int32_t>(0, primary_ref_begin - padding);
    int32_t final_ref_begin = primary_ref_region_begin + primary_swaln.ref_begin;
    int32_t final_ref_end = primary_ref_region_begin + primary_swaln.ref_end;

    result.mapped = true;
    result.best_score = primary_swaln.score;
    result.second_best_score = all_alignments.size() > 1 ? all_alignments[1].second.score : 0;

    // Build primary alignment
    result.primary.qname = std::string(read.name.view());
    result.primary.rname = ref_seq.name;
    result.primary.pos = final_ref_begin + 1;
    result.primary.mapq = compute_mapq(result.best_score, result.second_best_score);
    result.primary.score = primary_swaln.score;
    result.primary.seq = std::string(read.seq.view());
    result.primary.qual = std::string(read.qual.view());
    result.primary.cigar = primary_swaln.cigar;

    // Add NM and MD tags for primary
    int32_t nm = 0;
    for (uint32_t c : result.primary.cigar) {
        auto op = static_cast<bwa::align::CigarOp>(c & 0xF);
        int len = bwa::align::cigar_len(c);
        if (op == bwa::align::CigarOp::Diff) nm += len;
        else if (op == bwa::align::CigarOp::Ins || op == bwa::align::CigarOp::Del) nm += len;
    }
    result.primary.tags.push_back({"NM", std::to_string(nm)});

    // MD tag
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
            if (first) { md += std::to_string(ref_pos_in_aln); first = false; }
            else { md += std::to_string(run_len); }
            run_len = 0;
            for (int k = 0; k < len; ++k) {
                if (k > 0) md += "0";
                if (ref_pos_in_aln < ref_len) {
                    auto ref_base = index_.extract_ref(ref_pos_in_aln, ref_pos_in_aln + 1);
                    md += !ref_base.empty() ? bwa::index::PackedSequence::decode_base(ref_base[0]) : 'N';
                } else { md += 'N'; }
                ref_pos_in_aln++;
            }
        } else if (op == bwa::align::CigarOp::Del) {
            if (first) { md += std::to_string(ref_pos_in_aln); first = false; }
            else { md += std::to_string(run_len); }
            run_len = 0;
            md += '^';
            for (int k = 0; k < len; ++k) {
                if (ref_pos_in_aln < ref_len) {
                    auto ref_base = index_.extract_ref(ref_pos_in_aln, ref_pos_in_aln + 1);
                    md += !ref_base.empty() ? bwa::index::PackedSequence::decode_base(ref_base[0]) : 'N';
                } else { md += 'N'; }
                ref_pos_in_aln++;
            }
        }
    }
    if (first) { md = std::to_string(ref_pos_in_aln - final_ref_begin); }
    else { md += std::to_string(run_len); }
    result.primary.tags.push_back({"MD", md});

    // Add secondary alignments
    for (size_t i = 1; i < all_alignments.size() && i < 4; ++i) {
        const auto& [sec_chain_idx, sec_swaln] = all_alignments[i];
        int32_t sec_ref_begin = chains[sec_chain_idx].mems.front().ref_pos;
        int32_t sec_ref_region_begin = std::max<int32_t>(0, sec_ref_begin - padding);
        int32_t sec_final_begin = sec_ref_region_begin + sec_swaln.ref_begin;
        int32_t sec_final_end = sec_ref_region_begin + sec_swaln.ref_end;

        AlnRecord sec;
        sec.qname = std::string(read.name.view());
        sec.flag = AlnRecord::F_SECONDARY;
        sec.rname = ref_seq.name;
        sec.pos = sec_final_begin + 1;
        sec.mapq = compute_mapq(result.best_score, sec_swaln.score);
        sec.score = sec_swaln.score;
        sec.seq = std::string(read.seq.view());
        sec.qual = std::string(read.qual.view());
        sec.cigar = sec_swaln.cigar;
        result.secondary.push_back(std::move(sec));
    }
}

void Aligner::align_pair_impl(const bwa::io::SeqRecord& read1,
                              const bwa::io::SeqRecord& read2,
                              AlignmentResult& result) const {
    // Align both reads independently
    AlignmentResult result1, result2;
    align_impl(read1, result1);
    align_impl(read2, result2);

    // Use the better alignment as primary
    AlnRecord* primary;
    AlnRecord* mate;

    if (result1.best_score >= result2.best_score) {
        result = std::move(result1);
        primary = &result.primary;
        if (result2.mapped) {
            result.secondary.push_back(std::move(result2.primary));
            mate = &result.secondary.back();
        } else {
            mate = nullptr;
        }
    } else {
        result = std::move(result2);
        primary = &result.primary;
        if (result1.mapped) {
            result.secondary.push_back(std::move(result1.primary));
            mate = &result.secondary.back();
        } else {
            mate = nullptr;
        }
    }

    // Set paired-end FLAG bits and mate information
    if (result.mapped && mate) {
        // Set paired flags
        primary->flag = AlnRecord::F_PAIRED | AlnRecord::F_READ1;
        mate->flag = AlnRecord::F_PAIRED | AlnRecord::F_READ2;

        // Set mate reference and position
        primary->rnext = mate->rname;
        primary->pnext = mate->pos;
        mate->rnext = primary->rname;
        mate->pnext = primary->pos;

        // Check for proper pair
        if (primary->rname == mate->rname) {
            // Both mapped to same reference
            // Check orientation: FR (forward-reverse) is the typical case
            // For now, always mark as proper pair if on same reference
            primary->flag |= AlnRecord::F_PROPER_PAIR;
            mate->flag |= AlnRecord::F_PROPER_PAIR;

            // Calculate template length (approximate)
            int32_t tlen = mate->pos - primary->pos + 20; // +20 for read length
            primary->tlen = tlen;
            mate->tlen = -tlen;
        }
    } else if (result.mapped) {
        // Mate didn't map
        primary->flag = AlnRecord::F_PAIRED | AlnRecord::F_READ1 | AlnRecord::F_MUNMAP;
    }
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

uint8_t Aligner::compute_mapq(int32_t best, int32_t second_best) const noexcept {
    if (best <= 0) return 0;

    // BWA-MEM style MAPQ calculation
    // Based on score difference and uniqueness
    if (second_best <= 0) {
        // Unique alignment
        return 60;
    }

    int diff = best - second_best;

    // Scale by score difference
    if (diff >= best / 2) return 60;      // Very confident
    if (diff >= best / 3) return 40;      // Confident
    if (diff >= best / 4) return 30;      // Moderate
    if (diff >= 10) return 20;            // Low
    if (diff >= 5) return 10;             // Very low

    return 0;                              // Uncertain
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

    // Write secondary alignments
    for (const auto& sec : result.secondary) {
        write_sam_record(out, sec);
    }
}

void Pipeline::write_sam_record(std::ostream& out, const AlnRecord& aln) const {
    out << aln.qname << '\t';
    out << aln.flag << '\t';
    out << aln.rname << '\t';
    out << aln.pos << '\t';
    out << static_cast<int>(aln.mapq) << '\t';
    out << format_cigar(aln.cigar) << '\t';
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

    // Build metadata string
    std::ostringstream oss;
    oss << "BWA-CPP26-METADATA-v1\n";
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

    // Save FM-index (simplified - just BWT for now)
    std::string bwt_path = std::string(prefix) + ".bwt";
    std::ofstream bwt_out(bwt_path, std::ios::binary);
    for (const auto& fm : fm_index_) {
        bwt_out.write(reinterpret_cast<const char*>(fm.bwt().words().data()),
                      static_cast<std::streamsize>(fm.bwt().words().size() * sizeof(uint64_t)));
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

    // Save packed reference (for subsequence extraction)
    std::string pac_path = std::string(prefix) + ".pac";
    std::ofstream pac_out(pac_path, std::ios::binary);
    if (!pac_out) {
        throw std::runtime_error("Cannot open packed reference file for writing");
    }
    for (const auto& packed : packed_refs_) {
        pac_out.write(reinterpret_cast<const char*>(packed.words().data()),
                      static_cast<std::streamsize>(packed.words().size() * sizeof(uint64_t)));
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

    // Parse metadata
    std::string_view mv = meta_;
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
    const uint64_t* bwt_data = static_cast<const uint64_t*>(bwt_mmap.data());
    size_t bwt_words = bwt_size / sizeof(uint64_t);

    // Create FM index from loaded data (single index for now)
    {
        bwa::index::FMIndex idx;

        // Set BWT - directly use mmapped data
        bwa::index::PackedSequence bwt;
        bwt.resize(total_len);
        for (size_t i = 0; i < total_len; ++i) {
            // Unpack the 2-bit values from the uint64_t array
            size_t wi = i / 32;
            size_t bi = (i % 32) * 2;
            uint8_t val = (bwt_data[wi] >> bi) & 0x3;
            bwt.set(i, val);
        }
        idx.set_bwt(bwt);
        idx.set_size(total_len);

        // Load SA samples
        std::string sa_path = std::string(prefix) + ".sa";
        bwa::io::MmapFile sa_mmap;
        sa_mmap.open(sa_path.c_str());
        if (!sa_mmap.is_open()) throw std::runtime_error("Cannot open SA file");

        const uint32_t* sa_data = static_cast<const uint32_t*>(sa_mmap.data());
        size_t sa_words = sa_mmap.size() / sizeof(uint32_t);
        idx.set_sa_samples(std::vector<uint32_t>(sa_data, sa_data + sa_words));

        // Load occ table
        std::string occ_path = std::string(prefix) + ".occ";
        bwa::io::MmapFile occ_mmap;
        occ_mmap.open(occ_path.c_str());
        if (!occ_mmap.is_open()) throw std::runtime_error("Cannot open OCC file");

        const uint32_t* occ_data = static_cast<const uint32_t*>(occ_mmap.data());
        size_t occ_words = occ_mmap.size() / sizeof(uint32_t);
        idx.set_occ_table(std::vector<uint32_t>(occ_data, occ_data + occ_words));

        // Build count table from occ
        std::vector<uint32_t> cnt(5, 0);
        for (int b = 0; b < 4; ++b) {
            // Last interval has final count
            size_t last_intv = (total_len + 128 - 1) / 128;
            cnt[b + 1] = cnt[b] + occ_data[last_intv * 4 + b];
        }
        idx.set_count_table(cnt);

        // Add to fm_index_
        fm_index_.add_index(std::move(idx));
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

    // Distribute packed data among references
    packed_refs_.resize(refs_.size());
    size_t offset = 0;
    for (size_t i = 0; i < refs_.size(); ++i) {
        packed_refs_[i].resize(refs_[i].length);
        for (size_t j = 0; j < refs_[i].length; ++j) {
            size_t wi = offset / 32;
            size_t bi = (offset % 32) * 2;
            uint8_t val = (pac_data[wi] >> bi) & 0x3;
            packed_refs_[i].set(j, val);
            ++offset;
        }
    }
}

} // namespace bwa
