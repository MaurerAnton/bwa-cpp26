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
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <algorithm>

using namespace bwa;

namespace {

// Helper: format SAM CIGAR from encoded cigar vector
std::string format_cigar(const std::vector<uint32_t>& cigar) {
    std::string s;
    for (uint32_t c : cigar) {
        int len = align::cigar_len(c);
        char op = align::cigar_char(align::cigar_op(c));
        s += std::to_string(len);
        s += op;
    }
    return s;
}

// Helper: extract a subregion of the packed reference as uint8_t span
std::vector<uint8_t> extract_ref_region(const index::FMIndex& fm,
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

void Aligner::align_impl(const io::SeqRecord& read, AlignmentResult& result) const {
    result.clear();

    // Reset arena for this read
    memory::reset_tls_arena();

    // Pack query sequence
    index::PackedSequence query;
    query.append(read.seq.data(), read.seq.size());

    if (query.size() == 0) {
        result.mapped = false;
        return;
    }

    // Find MEMs (use TLS arena for memory)
    memory::Arena& arena = memory::get_tls_arena();
    auto mems = mem_finder_.find(query.bases(), arena);
    mem_finder_.filter_overlaps(mems);

    if (mems.empty()) {
        result.mapped = false;
        return;
    }

    // Chain MEMs
    auto chains = mem_finder_.chain(mems, config_.max_gap, config_.min_chain_score);

    if (chains.empty()) {
        result.mapped = false;
        return;
    }

    // Process best chain
    const auto& best_chain = chains[0];

    // Get reference sequence info
    if (index_.num_references() == 0) {
        result.mapped = false;
        return;
    }

    const auto& ref_seq = index_.references()[0];
    int32_t ref_len = static_cast<int32_t>(ref_seq.length);
    int32_t query_len = static_cast<int32_t>(read.seq.size());

    // Pack query for SW
    std::vector<uint8_t> query_bytes(query_len);
    auto qbases = query.bases();
    for (int32_t i = 0; i < query_len; ++i) {
        query_bytes[i] = qbases[i];
    }

    // Adaptive scoring based on read length (BWA-MEM style)
    // For longer reads, increase gap penalties to avoid spurious alignments
    align::Scoring adaptive_scoring = config_.scoring;
    if (query_len > 100) {
        // Scale gap penalties with read length
        adaptive_scoring.gap_open = std::min<int>(-1, static_cast<int>(-0.01 * query_len - 4));
        adaptive_scoring.gap_ext = std::min<int>(-1, static_cast<int>(-0.1 * query_len));
    }

    // Process up to max_secondary alignments (best chain + secondary chains)
    int max_alignments = std::min<int>(chains.size(), 1 + config_.max_occ / 100);
    if (max_alignments < 1) max_alignments = 1;
    if (max_alignments > 4) max_alignments = 4; // Cap at 4 total alignments

    align::Alignment best_swaln;
    int32_t best_chain_idx = -1;
    std::vector<std::pair<int32_t, align::Alignment>> all_alignments;

    for (int32_t ci = 0; ci < max_alignments && ci < (int32_t)chains.size(); ++ci) {
        const auto& chain = chains[ci];
        int32_t chain_ref_begin = chain.mems.front().ref_pos;

        int32_t padding = config_.band_width * 2;
        int32_t ref_begin = std::max<int32_t>(0, chain_ref_begin - padding);
        int32_t ref_end = std::min<int32_t>(ref_len, chain_ref_begin + query_len + padding);

        if (ref_end <= ref_begin) continue;

        auto ref_region = index_.extract_ref(ref_begin, ref_end);
        if (ref_region.empty()) continue;

        align::Alignment sw_aln = align::sw_semi_global_extend(
            adaptive_scoring,
            std::span<const uint8_t>(query_bytes.data(), query_len),
            std::span<const uint8_t>(ref_region.data(), ref_region.size()),
            config_.band_width
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
        result.primary.seq = std::string(read.seq.view());
        result.primary.qual = std::string(read.qual.view());

        int32_t prev_query_end = 0;
        for (const auto& mem : best_chain.mems) {
            if (mem.query_pos > prev_query_end) {
                result.primary.cigar.push_back(
                    align::encode_cigar(mem.query_pos - prev_query_end, align::CigarOp::SoftClip));
            }
            result.primary.cigar.push_back(
                align::encode_cigar(mem.len, align::CigarOp::Match));
            prev_query_end = mem.query_end();
        }
        if (prev_query_end < query_len) {
            result.primary.cigar.push_back(
                align::encode_cigar(query_len - prev_query_end, align::CigarOp::SoftClip));
        }
        return;
    }

    // Sort alignments by score (descending)
    std::sort(all_alignments.begin(), all_alignments.end(),
              [](const auto& a, const auto& b) { return a.second.score > b.second.score; });

    // Use the best alignment as primary
    const auto& [primary_chain_idx, primary_swaln] = all_alignments[0];
    int32_t primary_ref_begin = chains[primary_chain_idx].mems.front().ref_pos;
    int32_t padding = config_.band_width * 2;
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
        auto op = static_cast<align::CigarOp>(c & 0xF);
        int len = align::cigar_len(c);
        if (op == align::CigarOp::Diff) nm += len;
        else if (op == align::CigarOp::Ins || op == align::CigarOp::Del) nm += len;
    }
    result.primary.tags.push_back({"NM", std::to_string(nm)});

    // MD tag
    std::string md;
    int32_t ref_pos_in_aln = final_ref_begin;
    int32_t run_len = 0;
    bool first = true;
    for (uint32_t c : result.primary.cigar) {
        auto op = static_cast<align::CigarOp>(c & 0xF);
        int len = align::cigar_len(c);
        if (op == align::CigarOp::Equal) {
            run_len += len;
            ref_pos_in_aln += len;
        } else if (op == align::CigarOp::Diff) {
            if (first) { md += std::to_string(ref_pos_in_aln); first = false; }
            else { md += std::to_string(run_len); }
            run_len = 0;
            for (int k = 0; k < len; ++k) {
                if (k > 0) md += "0";
                if (ref_pos_in_aln < ref_len) {
                    auto ref_base = index_.extract_ref(ref_pos_in_aln, ref_pos_in_aln + 1);
                    md += !ref_base.empty() ? index::PackedSequence::decode_base(ref_base[0]) : 'N';
                } else { md += 'N'; }
                ref_pos_in_aln++;
            }
        } else if (op == align::CigarOp::Del) {
            if (first) { md += std::to_string(ref_pos_in_aln); first = false; }
            else { md += std::to_string(run_len); }
            run_len = 0;
            md += '^';
            for (int k = 0; k < len; ++k) {
                if (ref_pos_in_aln < ref_len) {
                    auto ref_base = index_.extract_ref(ref_pos_in_aln, ref_pos_in_aln + 1);
                    md += !ref_base.empty() ? index::PackedSequence::decode_base(ref_base[0]) : 'N';
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

void Aligner::align_pair_impl(const io::SeqRecord& read1,
                              const io::SeqRecord& read2,
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

void Aligner::chain_to_alignment(const align::MEMFinder::Chain& chain,
                                 const io::SeqRecord& read,
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
                align::encode_cigar(clip_len, align::CigarOp::SoftClip));
        }

        // Add the MEM as a match
        aln.cigar.push_back(
            align::encode_cigar(mem.len, align::CigarOp::Match));

        prev_query_end = mem.query_end();
    }

    // Add trailing soft clip
    if (prev_query_end < query_len) {
        int32_t clip_len = query_len - prev_query_end;
        aln.cigar.push_back(
            align::encode_cigar(clip_len, align::CigarOp::SoftClip));
    }

    aln.score = chain.score;
}

uint8_t Aligner::compute_mapq(int32_t best, int32_t second_best) const noexcept {
    if (second_best <= 0) return 60;
    int diff = best - second_best;
    if (diff >= 40) return 60;
    return static_cast<uint8_t>(std::max(1, diff * 60 / 40));
}

void Pipeline::write_header(std::ostream& out) const {
    out << "@HD\tVN:1.6\tSO:coordinate\n";

    for (size_t i = 0; i < index_.num_references(); ++i) {
        const auto& ref = index_.get_ref(i);
        out << "@SQ\tSN:" << ref.name
            << "\tLN:" << ref.length << "\n";
    }
}

void Pipeline::write_alignment(std::ostream& out, const AlignmentResult& result) const {
    if (!result.mapped) {
        return;
    }

    // Write primary alignment
    write_sam_record(out, result.primary);

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
    out << '\n';
}

void Index::build_impl(const char* fasta_path, const Config& cfg) {
    io::SeqReader reader(fasta_path);
    if (!reader.is_open()) {
        throw std::runtime_error("Failed to open FASTA file");
    }

    io::SeqRecord rec;
    size_t total_len = 0;

    auto read_result = reader.read(rec);
    while (read_result && *read_result) {
        if (rec.is_fasta() && !rec.seq.empty()) {
            RefSequence ref;
            ref.name = std::string(rec.name.view());
            ref.length = rec.seq.size();
            ref.offset = total_len;

            // Pack sequence
            index::PackedSequence seq;
            seq.append(rec.seq.data(), rec.seq.size());

            // Store packed reference for later extraction
            packed_refs_.push_back(seq);

            fm_index_.add_sequence(seq, memory::get_tls_arena());

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
        oss << ref.name << "\t" << ref.length << "\t" << ref.offset << "\n";
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
        end = mv.find('\n', pos);
        ref.offset = std::stoull(std::string(mv.substr(pos, end - pos)));
        pos = end + 1;
        refs_.push_back(std::move(ref));
    }

    // Load BWT
    std::string bwt_path = std::string(prefix) + ".bwt";
    std::ifstream bwt_in(bwt_path, std::ios::binary);
    if (!bwt_in) throw std::runtime_error("Cannot open BWT file");

    // Determine size and read into buffer
    bwt_in.seekg(0, std::ios::end);
    size_t bwt_size = bwt_in.tellg();
    bwt_in.seekg(0);
    std::vector<uint64_t> bwt_data(bwt_size / sizeof(uint64_t));
    bwt_in.read(reinterpret_cast<char*>(bwt_data.data()), bwt_size);
    bwt_in.close();

    // Create FM index from loaded data (single index for now)
    {
        index::FMIndex idx;

        // Set BWT
        index::PackedSequence bwt;
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
        std::ifstream sa_in(sa_path, std::ios::binary);
        sa_in.seekg(0, std::ios::end);
        size_t sa_size = sa_in.tellg();
        sa_in.seekg(0);
        std::vector<uint32_t> sa_data(sa_size / sizeof(uint32_t));
        sa_in.read(reinterpret_cast<char*>(sa_data.data()), sa_size);
        sa_in.close();
        idx.set_sa_samples(sa_data);

        // Load occ table
        std::string occ_path = std::string(prefix) + ".occ";
        std::ifstream occ_in(occ_path, std::ios::binary);
        occ_in.seekg(0, std::ios::end);
        size_t occ_size = occ_in.tellg();
        occ_in.seekg(0);
        std::vector<uint32_t> occ_data(occ_size / sizeof(uint32_t));
        occ_in.read(reinterpret_cast<char*>(occ_data.data()), occ_size);
        occ_in.close();
        idx.set_occ_table(occ_data);

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
    std::ifstream pac_in(pac_path, std::ios::binary);
    if (!pac_in) {
        // No packed reference - extract_ref will return empty
        return;
    }
    pac_in.seekg(0, std::ios::end);
    size_t pac_size = pac_in.tellg();
    pac_in.seekg(0);
    std::vector<uint64_t> pac_data(pac_size / sizeof(uint64_t));
    pac_in.read(reinterpret_cast<char*>(pac_data.data()), pac_size);
    pac_in.close();

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
