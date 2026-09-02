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

    // Get reference position from the first MEM
    int32_t ref_pos = best_chain.mems[0].ref_pos;
    int32_t chain_query_begin = best_chain.mems.front().query_pos;
    int32_t chain_ref_begin = best_chain.mems.front().ref_pos;

    // Calculate alignment region with padding
    int32_t padding = config_.band_width * 2;
    int32_t ref_begin = std::max<int32_t>(0, chain_ref_begin - padding);
    int32_t ref_end = std::min<int32_t>(ref_len, chain_ref_begin + query_len + padding);

    if (ref_end <= ref_begin) {
        result.mapped = false;
        return;
    }

    // Extract reference region
    auto ref_region = index_.extract_ref(ref_begin, ref_end);
    if (ref_region.empty()) {
        // No packed reference available - fall back to simple CIGAR from chain
        result.mapped = true;
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
                int32_t clip_len = mem.query_pos - prev_query_end;
                result.primary.cigar.push_back(
                    align::encode_cigar(clip_len, align::CigarOp::SoftClip));
            }
            result.primary.cigar.push_back(
                align::encode_cigar(mem.len, align::CigarOp::Match));
            prev_query_end = mem.query_end();
        }
        if (prev_query_end < query_len) {
            int32_t clip_len = query_len - prev_query_end;
            result.primary.cigar.push_back(
                align::encode_cigar(clip_len, align::CigarOp::SoftClip));
        }
        return;
    }

    // Convert query to span for SW
    std::vector<uint8_t> query_bytes(query_len);
    auto qbases = query.bases();
    for (int32_t i = 0; i < query_len; ++i) {
        query_bytes[i] = qbases[i];
    }

    // Run banded SW extension
    align::Alignment sw_aln = align::sw_extend(
        config_.scoring,
        std::span<const uint8_t>(query_bytes.data(), query_len),
        std::span<const uint8_t>(ref_region.data(), ref_region.size()),
        0, 0, config_.band_width, config_.max_score_drop
    );

    // Calculate final alignment position
    int32_t final_ref_begin = ref_begin + sw_aln.ref_begin;
    int32_t final_ref_end = ref_begin + sw_aln.ref_end;

    result.mapped = sw_aln.score > 0;
    result.best_score = sw_aln.score;
    result.second_best_score = chains.size() > 1 ? chains[1].score : 0;

    // Build primary alignment
    result.primary.qname = std::string(read.name.view());
    result.primary.rname = ref_seq.name;
    result.primary.pos = final_ref_begin + 1; // SAM is 1-based
    result.primary.mapq = compute_mapq(result.best_score, result.second_best_score);
    result.primary.score = sw_aln.score;
    result.primary.seq = std::string(read.seq.view());
    result.primary.qual = std::string(read.qual.view());

    // Copy CIGAR from SW alignment
    result.primary.cigar = sw_aln.cigar;
}

void Aligner::align_pair_impl(const io::SeqRecord& read1,
                              const io::SeqRecord& read2,
                              AlignmentResult& result) const {
    // Align first read
    align_impl(read1, result);

    // For paired-end, we would also align read2 and check insert size
    // For now, this is a stub
    AlignmentResult result2;
    align_impl(read2, result2);
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
    if (!result.mapped) return;

    // Build SAM line
    out << result.primary.qname << '\t';
    out << result.primary.flag << '\t';
    out << result.primary.rname << '\t';
    out << result.primary.pos << '\t';
    out << static_cast<int>(result.primary.mapq) << '\t';
    out << format_cigar(result.primary.cigar) << '\t';
    out << result.primary.rnext << '\t'
        << result.primary.pnext << '\t'
        << result.primary.tlen << '\t';
    out << result.primary.seq << '\t';
    out << result.primary.qual << '\n';
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
