#pragma once

#include <bwa/index/fm_index.hpp>
#include <bwa/align/mem.hpp>
#include <bwa/align/sw.hpp>
#include <bwa/io/seq_io.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/hash_map.hpp>
#include <span>
#include <string_view>
#include <optional>
#include <functional>
#include <chrono>

namespace bwa {

// Main configuration
struct Config {
    // Algorithm parameters
    int min_seed_len = 19;        // Minimum seed length
    int max_occ = 500;            // Max occurrences for a seed
    int max_gap = 10000;          // Max gap in chaining
    int min_chain_score = 30;     // Minimum chain score
    int max_chain_gap = 10000;    // Max gap between chains
    int band_width = 32;          // Band width for SW extension
    int max_score_drop = 100;     // Max score drop for early termination

    // Scoring
    align::Scoring scoring = align::Scoring::bwa_mem_default();

    // Threading
    int num_threads = 1;

    // Output
    bool output_sam = true;
    bool output_secondary = true;
    bool output_supplementary = true;
    int min_mapq = 1;             // Minimum MAPQ

    // Memory
    size_t batch_size = 1000;     // Reads per batch

    [[nodiscard]] static Config default_mem() noexcept {
        return Config{};
    }

    [[nodiscard]] static Config fast() noexcept {
        Config c;
        c.min_seed_len = 23;
        c.max_occ = 200;
        c.band_width = 16;
        return c;
    }

    [[nodiscard]] static Config accurate() noexcept {
        Config c;
        c.min_seed_len = 15;
        c.max_occ = 1000;
        c.band_width = 64;
        c.max_score_drop = 200;
        return c;
    }
};

// Read group / sample info
struct ReadGroup {
    core::PmrString id;
    core::PmrString sample;
    core::PmrString library;
    core::PmrString platform;
    core::PmrString platform_unit;
};

// Alignment record (SAM-compatible)
struct AlnRecord {
    core::PmrString qname;          // Query name
    uint32_t flag = 0;        // SAM flag
    core::PmrString rname;          // Reference name
    int32_t pos = 0;          // 1-based position
    uint8_t mapq = 0;         // MAPQ
    core::Vector<uint32_t> cigar;   // CIGAR
    core::PmrString rnext = "*";    // Mate reference
    int32_t pnext = 0;        // Mate position
    int32_t tlen = 0;         // Template length
    core::PmrString seq;            // Query sequence
    core::PmrString qual;           // Query quality
    core::Vector<std::pair<core::PmrString, core::PmrString>> tags; // Optional tags

    // Computed fields
    int32_t score = 0;
    int32_t secondary_score = 0;
    bool is_primary = true;
    bool is_supplementary = false;

    void clear() noexcept {
        qname.clear(); flag = 0; rname.clear(); pos = 0; mapq = 0;
        cigar.clear(); rnext = "*"; pnext = 0; tlen = 0;
        seq.clear(); qual.clear(); tags.clear();
        score = secondary_score = 0;
        is_primary = true; is_supplementary = false;
    }

    // SAM flag constants
    static constexpr uint32_t F_PAIRED = 0x1;
    static constexpr uint32_t F_PROPER_PAIR = 0x2;
    static constexpr uint32_t F_UNMAP = 0x4;
    static constexpr uint32_t F_MUNMAP = 0x8;
    static constexpr uint32_t F_REVERSE = 0x10;
    static constexpr uint32_t F_MREVERSE = 0x20;
    static constexpr uint32_t F_READ1 = 0x40;
    static constexpr uint32_t F_READ2 = 0x80;
    static constexpr uint32_t F_SECONDARY = 0x100;
    static constexpr uint32_t F_QCFAIL = 0x200;
    static constexpr uint32_t F_DUP = 0x400;
    static constexpr uint32_t F_SUPPLEMENTARY = 0x800;
};

// Reference sequence info
struct RefSequence {
    std::string name;   // Heap-allocated, persistent
    std::string md5;
    size_t length = 0;
    size_t offset = 0; // Offset in concatenated reference
};

// Main BWA index
// Uses heap allocators for metadata to avoid arena lifetime issues.
// FMIndex data is stored in arenas owned by the FMIndex objects themselves.
class Index {
    index::MultiFMIndex fm_index_;
    std::vector<RefSequence> refs_;  // Heap-allocated, persistent
    std::string meta_;               // Heap-allocated, persistent

public:
    Index() = default;

    // Build from FASTA files
    static Index build(const char* fasta_path,
                       const Config& cfg = Config::default_mem()) {
        Index idx;
        idx.build_impl(fasta_path, cfg);
        return idx;
    }

    // Load from disk
    static Index load(const char* prefix) {
        Index idx;
        idx.load_impl(prefix);
        return idx;
    }

    void save(const char* prefix) const { save_impl(prefix); }

    // Accessors
    [[nodiscard]] const index::MultiFMIndex& fm_index() const noexcept { return fm_index_; }
    [[nodiscard]] index::MultiFMIndex& fm_index() noexcept { return fm_index_; }
    [[nodiscard]] const std::vector<RefSequence>& references() const noexcept { return refs_; }
    [[nodiscard]] size_t num_references() const noexcept { return refs_.size(); }
    [[nodiscard]] size_t total_length() const noexcept { return fm_index_.total_length(); }

    [[nodiscard]] const RefSequence& get_ref(size_t i) const noexcept { return refs_[i]; }
    [[nodiscard]] std::optional<size_t> find_ref(std::string_view name) const noexcept {
        for (size_t i = 0; i < refs_.size(); ++i) {
            if (refs_[i].name == name) return i;
        }
        return std::nullopt;
    }

private:
    void build_impl(const char* fasta_path, const Config& cfg);
    void load_impl(const char* prefix);
    void save_impl(const char* prefix) const;
};

// Alignment result for a single read
struct AlignmentResult {
    AlnRecord primary;
    core::Vector<AlnRecord> secondary;
    core::Vector<AlnRecord> supplementary;
    bool mapped = false;
    int32_t best_score = 0;
    int32_t second_best_score = 0;

    void clear() noexcept {
        primary.clear(); secondary.clear(); supplementary.clear();
        mapped = false; best_score = second_best_score = 0;
    }
};

// Main aligner class
class Aligner {
    const Index& index_;
    Config config_;
    align::MEMFinder mem_finder_;
    mutable memory::Arena arena_;

public:
    Aligner(const Index& idx, const Config& cfg = Config::default_mem())
        : index_(idx), config_(cfg),
          mem_finder_(idx.fm_index()[0], cfg.min_seed_len, cfg.max_occ, cfg.scoring.match),
          arena_(64 * 1024) {}

    // Align single read
    AlignmentResult align(const io::SeqRecord& read) const {
        AlignmentResult result;
        align_impl(read, result);
        return result;
    }

    // Align paired-end reads
    AlignmentResult align_pair(const io::SeqRecord& read1,
                               const io::SeqRecord& read2) const {
        AlignmentResult result;
        align_pair_impl(read1, read2, result);
        return result;
    }

    // Align batch of reads
    void align_batch(std::span<const io::SeqRecord> reads,
                     std::span<AlignmentResult> results) const {
        for (size_t i = 0; i < reads.size() && i < results.size(); ++i) {
            results[i].clear();
            align_impl(reads[i], results[i]);
        }
    }

    // Streaming align (callback per result)
    template <typename Callback>
    void align_stream(io::SeqReader& reader, Callback&& cb) const {
        io::SeqRecord read;
        while (reader.read(read)) {
            AlignmentResult result;
            align_impl(read, result);
            cb(result);
            memory::reset_tls_arena();
        }
    }

    [[nodiscard]] const Config& config() const noexcept { return config_; }
    [[nodiscard]] const Index& index() const noexcept { return index_; }

private:
    void align_impl(const io::SeqRecord& read, AlignmentResult& result) const;
    void align_pair_impl(const io::SeqRecord& read1,
                         const io::SeqRecord& read2,
                         AlignmentResult& result) const;

    // Convert MEM chain to alignment
    void chain_to_alignment(const align::MEMFinder::Chain& chain,
                            const io::SeqRecord& read,
                            const RefSequence& ref,
                            AlnRecord& aln) const;

    // Compute MAPQ from scores
    uint8_t compute_mapq(int32_t best, int32_t second_best) const noexcept;
};

// High-level pipeline
class Pipeline {
    Index index_;
    Config config_;
    mutable Aligner aligner_;

public:
    Pipeline(const char* index_prefix, const Config& cfg = Config::default_mem())
        : index_(Index::load(index_prefix)), config_(cfg), aligner_(index_, config_) {}

    static Pipeline build(const char* fasta_path,
                          const char* index_prefix,
                          const Config& cfg = Config::default_mem()) {
        Index idx = Index::build(fasta_path, cfg);
        idx.save(index_prefix);
        return Pipeline(index_prefix, cfg);
    }

    // Align FASTQ file to SAM output
    void align_file(const char* fastq_path, const char* sam_path = "-") const {
        io::SeqReader reader(fastq_path);
        io::SeqWriter writer(sam_path, "wb");
        write_header(writer);

        aligner_.align_stream(reader, [&](const AlignmentResult& result) {
            write_alignment(writer, result);
        });
    }

    // Align paired FASTQ files
    void align_pair(const char* fastq1, const char* fastq2, const char* sam_path = "-") const {
        io::SeqReader r1(fastq1), r2(fastq2);
        io::SeqWriter writer(sam_path, "wb");
        write_header(writer);

        io::SeqRecord read1, read2;
        while (r1.read(read1) && r2.read(read2)) {
            AlignmentResult result = aligner_.align_pair(read1, read2);
            write_alignment(writer, result);
            memory::reset_tls_arena();
        }
    }

private:
    void write_header(const io::SeqWriter& writer) const;
    void write_alignment(const io::SeqWriter& writer, const AlignmentResult& result) const;
};

} // namespace bwa