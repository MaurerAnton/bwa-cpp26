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
#include <chrono>
#include <cstdlib>
#include <algorithm>

using namespace bwa;

void Aligner::align_impl(const io::SeqRecord& read, AlignmentResult& result) const {
    // Pack query sequence
    index::PackedSequence query;
    query.append(read.seq.data(), read.seq.size());

    // Find MEMs
    auto mems = mem_finder_.find(query.bases(), arena_);
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

    // Extend best chain with Smith-Waterman
    const auto& best_chain = chains[0];

    // Get reference sequence
    const auto& ref_seq = index_.references()[0]; // Single reference for now

    // For now, just score based on chain
    result.mapped = true;
    result.best_score = best_chain.score;
    result.second_best_score = chains.size() > 1 ? chains[1].score : 0;

    // Build primary alignment
    result.primary.qname = read.name;
    result.primary.rname = ref_seq.name;
    result.primary.pos = best_chain.ref_begin + 1; // SAM is 1-based
    result.primary.mapq = compute_mapq(result.best_score, result.second_best_score);
    result.primary.score = best_chain.score;

    // Create simple CIGAR (placeholder)
    result.primary.cigar.push_back(align::encode_cigar(best_chain.mems[0].len, align::CigarOp::Match));
}

void Aligner::align_pair_impl(const io::SeqRecord& read1,
                              const io::SeqRecord& read2,
                              AlignmentResult& result) const {
    // Align first read
    align_impl(read1, result);
    
    // Align second read
    AlignmentResult result2;
    align_impl(read2, result2);
    
    // Combine (simplified)
    if (result.mapped && result2.mapped) {
        // Paired-end logic
    }
}

void Aligner::chain_to_alignment(const align::MEMFinder::Chain& chain,
                                 const io::SeqRecord& read,
                                 const RefSequence& ref,
                                 AlnRecord& aln) const {
    // Build CIGAR from chain
    aln.cigar.clear();
    for (const auto& mem : chain.mems) {
        aln.cigar.push_back(align::encode_cigar(mem.len, align::CigarOp::Match));
    }
}

uint8_t Aligner::compute_mapq(int32_t best, int32_t second_best) const noexcept {
    if (second_best <= 0) return 60;
    int diff = best - second_best;
    if (diff >= 40) return 60;
    return static_cast<uint8_t>(std::max(1, diff * 60 / 40));
}

void Pipeline::write_header(const io::SeqWriter& writer) const {
    // Write SAM header
    core::PmrString header;
    header.kputs("@HD\tVN:1.6\tSO:coordinate\n");
    
    for (size_t i = 0; i < index_.num_references(); ++i) {
        const auto& ref = index_.get_ref(i);
        header.kputs("@SQ\tSN:");
        header.kputsn(ref.name.data(), ref.name.size());
        header.kputs("\tLN:");
        header.kputl(static_cast<int64_t>(ref.length));
        header.kputc('\n');
    }
    
    // Write to output (placeholder - would use writer)
}

void Pipeline::write_alignment(const io::SeqWriter& writer, const AlignmentResult& result) const {
    if (!result.mapped) return;
    
    // Write SAM line for primary
    core::PmrString line;
    line.kputsn(result.primary.qname.data(), result.primary.qname.size());
    line.kputc('\t');
    line.kputw(result.primary.flag);
    line.kputc('\t');
    line.kputsn(result.primary.rname.data(), result.primary.rname.size());
    line.kputc('\t');
    line.kputl(result.primary.pos);
    line.kputc('\t');
    line.kputw(result.primary.mapq);
    line.kputc('\t');
    
    // CIGAR
    for (size_t i = 0; i < result.primary.cigar.size(); ++i) {
        if (i > 0) line.kputc('\t');
        uint32_t c = result.primary.cigar[i];
        line.kputw(align::cigar_len(c));
        line.kputc(align::cigar_char(align::cigar_op(c)));
    }
    line.kputc('\t');
    
    line.kputs("*\t0\t0\t"); // rnext, pnext, tlen
    line.kputsn(result.primary.seq.data(), result.primary.seq.size());
    line.kputc('\t');
    line.kputsn(result.primary.qual.data(), result.primary.qual.size());
    line.kputc('\n');
    
    // Write to output (placeholder)
}

void Index::build_impl(const char* fasta_path, const Config& cfg, memory::Arena& arena) {
    io::SeqReader reader(fasta_path);
    if (!reader.is_open()) {
        throw std::runtime_error("Failed to open FASTA file");
    }

    io::SeqRecord rec;
    size_t total_len = 0;

    while (reader.read(rec)) {
        if (rec.is_fasta() && !rec.seq.empty()) {
            RefSequence ref;
            ref.name = std::move(rec.name);
            ref.length = rec.seq.size();
            ref.offset = total_len;

            // Pack sequence
            index::PackedSequence seq;
            seq.append(rec.seq.data(), rec.seq.size());
            fm_index_.add_sequence(seq, arena);

            refs_.push_back(std::move(ref));
            total_len += rec.seq.size();
        }
    }

    // Build metadata string for saving
    core::PmrString meta;
    meta.kputs("BWA-CPP26_INDEX\n");
    meta.kputl(static_cast<int64_t>(refs_.size()));
    meta.kputc('\n');
    meta.kputl(static_cast<int64_t>(total_len));
    meta.kputc('\n');
    for (const auto& ref : refs_) {
        meta.kputsn(ref.name.data(), ref.name.size());
        meta.kputc('\t');
        meta.kputl(static_cast<int64_t>(ref.length));
        meta.kputc('\t');
        meta.kputl(static_cast<int64_t>(ref.offset));
        meta.kputc('\n');
    }
    meta_ = std::move(meta);
}

void Index::save_impl(const char* prefix) const {
    // Save metadata
    core::PmrString meta_path;
    meta_path.kputs(prefix);
    meta_path.kputs(".meta");
    std::ofstream meta_out(meta_path.data(), std::ios::binary);
    meta_out.write(meta_.data(), meta_.size());
    meta_out.close();

    // Save FM-index (simplified - just BWT for now)
    core::PmrString bwt_path;
    bwt_path.kputs(prefix);
    bwt_path.kputs(".bwt");
    std::ofstream bwt_out(bwt_path.data(), std::ios::binary);
    for (const auto& fm : fm_index_) {
        bwt_out.write(reinterpret_cast<const char*>(fm.bwt().words().data()),
                      fm.bwt().words().size() * sizeof(uint64_t));
    }
    bwt_out.close();

    // Save SA samples
    core::PmrString sa_path;
    sa_path.kputs(prefix);
    sa_path.kputs(".sa");
    std::ofstream sa_out(sa_path.data(), std::ios::binary);
    for (const auto& fm : fm_index_) {
        sa_out.write(reinterpret_cast<const char*>(fm.sa_samples().data()),
                     fm.sa_samples().size() * sizeof(uint32_t));
    }
    sa_out.close();

    // Save occ
    core::PmrString occ_path;
    occ_path.kputs(prefix);
    occ_path.kputs(".occ");
    std::ofstream occ_out(occ_path.data(), std::ios::binary);
    for (const auto& fm : fm_index_) {
        occ_out.write(reinterpret_cast<const char*>(fm.occ_table().data()),
                      fm.occ_table().size() * sizeof(uint32_t));
    }
    occ_out.close();
}

void Index::load_impl(const char* prefix) {
    // Load metadata
    core::PmrString meta_path;
    meta_path.kputs(prefix);
    meta_path.kputs(".meta");
    std::ifstream meta_in(meta_path.data(), std::ios::binary);
    if (!meta_in) throw std::runtime_error("Cannot open index metadata");

    meta_in.seekg(0, std::ios::end);
    size_t meta_size = meta_in.tellg();
    meta_in.seekg(0);
    meta_.resize(meta_size);
    meta_in.read(meta_.data(), meta_size);
    meta_in.close();

    // Parse metadata (simplified)
    std::string_view mv = meta_.view();
    // Skip header
    size_t pos = mv.find('\n') + 1;
    size_t num_refs = 0, total_len = 0;
    sscanf(mv.data() + pos, "%zu\n%zu\n", &num_refs, &total_len);

    refs_.reserve(num_refs);
    for (size_t i = 0; i < num_refs; ++i) {
        RefSequence ref;
        pos = mv.find('\n', pos) + 1;
        size_t end = mv.find('\t', pos);
        ref.name.assign(mv.substr(pos, end - pos));
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
    core::PmrString bwt_path;
    bwt_path.kputs(prefix);
    bwt_path.kputs(".bwt");
    std::ifstream bwt_in(bwt_path.data(), std::ios::binary);
    // ... load each FM-index BWT
    bwt_in.close();

    // Similar for SA, occ
}