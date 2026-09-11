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
#include <utility>
#include <chrono>
#include <iostream>
#include <thread>
#include <mutex>
#include <atomic>
#include <queue>
#include <condition_variable>
#include <sstream>
#include <zlib.h>
#include <bwa/io/bam_io.hpp>
#include <bwa/io/bai_io.hpp>

namespace bwa {

// Read group / sample info
struct ReadGroup {
    std::string id;
    std::string sample;
    std::string library;
    std::string platform;
    std::string platform_unit;
};

// Insert-size distribution estimated from a sample of read pairs
// (BWA mem_pestat equivalent). Used for proper-pair classification and
// mate-rescue windows; falls back to config bounds when invalid.
struct InsertSizeStats {
    double mean = 0.0;
    double std = 0.0;
    bool valid = false;
};

// Main configuration
struct Config {
    // Algorithm parameters
    int min_seed_len = 9;         // Minimum seed length (lowered for sensitivity)
    int max_occ = 500;            // Max occurrences for a seed
    int max_gap = 10000;          // Max gap in chaining
    int min_chain_score = 10;     // Minimum chain score
    int max_chain_gap = 10000;    // Max gap between chains
    int band_width = 32;          // Band width for SW extension
    int max_score_drop = 100;     // Max score drop for early termination

    // Scoring
    align::Scoring scoring = align::Scoring::bwa_mem_default();

    // Threading
    int num_threads = 1;

    // Output
    bool output_sam = true;
    // Secondary alignments are suppressed by default (BWA mem behavior);
    // the -a flag enables them.
    bool output_secondary = false;
    bool output_supplementary = true;
    int min_mapq = 1;             // Minimum MAPQ
    // Mark split (supplementary) hits as secondary instead (-M), for
    // pipelines that require Picard-compatible flags.
    bool mark_split_secondary = false;
    // Minimum alignment score to output (-T); 0 disables the filter.
    int min_output_score = 0;

    // Pairing controls (BWA -S / -P)
    bool skip_mate_rescue = false;
    bool skip_pairing = false;
    // -p: consecutive records in one FASTQ are mates (interleaved input).
    bool smart_pairing = false;
    // -Y: soft-clip (not hard-clip) split/supplementary segments.
    bool soft_clip_supplementary = false;
    // -5: for split alignments, the smallest-coordinate hit is primary.
    bool smallest_coord_primary = false;

    // Read group (optional, written as @RG header)
    std::optional<ReadGroup> read_group;

    // Program info (written as @PG header)
    std::string program_name = "bwa-cpp26";
    std::string program_version = "0.1.0";
    std::string program_command = "bwa-cpp26 mem";

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

    // Preset for long reads (ONT/PacBio)
    [[nodiscard]] static Config long_reads() noexcept {
        Config c;
        c.min_seed_len = 9;         // Shorter seeds for higher sensitivity
        c.max_occ = 10000;         // Allow more occurrences
        c.max_gap = 50000;         // Larger gaps
        c.min_chain_score = 20;    // Lower threshold
        c.band_width = 128;        // Wider band for more indels
        c.max_score_drop = 500;    // Allow more score drop
        // More permissive scoring for long reads with higher error rates
        c.scoring.match = 1;
        c.scoring.mismatch = -1;   // Less penalty for mismatches
        c.scoring.gap_open = -2;   // Less penalty for gaps
        c.scoring.gap_ext = -1;
        return c;
    }

    // Apply a BWA -x read-type preset in place. Later options override the
    // preset values, so call this as the options are parsed in order.
    // Returns false for unknown preset names.
    bool apply_preset(std::string_view name) noexcept {
        if (name == "ont2d" || name == "ont" || name == "nanopore") {
            min_seed_len = 14;
            scoring = {1, -1, -1, -1, 0, scoring.unpaired};
            return true;
        }
        if (name == "pacbio" || name == "pb" || name == "ccs") {
            min_seed_len = 17;
            scoring = {1, -1, -1, -1, 0, scoring.unpaired};
            return true;
        }
        if (name == "intractable" || name == "sr" || name == "short") {
            min_seed_len = 19;
            scoring = align::Scoring::bwa_mem_default();
            return true;
        }
        return false;
    }
};

// Alignment record (SAM-compatible)
struct AlnRecord {
    std::string qname;          // Query name
    uint32_t flag = 0;        // SAM flag
    std::string rname;          // Reference name
    int32_t pos = 0;          // 1-based position
    uint8_t mapq = 0;         // MAPQ
    std::vector<uint32_t> cigar;   // CIGAR
    std::string rnext = "*";    // Mate reference
    int32_t pnext = 0;        // Mate position
    int32_t tlen = 0;         // Template length
    std::string seq;            // Query sequence
    std::string qual;           // Query quality
    std::vector<std::pair<std::string, std::string>> tags; // Optional tags
    int32_t suboptimal_score = 0;  // XS tag: suboptimal alignment score

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
        suboptimal_score = 0;
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
};// Main BWA index
// Uses heap allocators for metadata to avoid arena lifetime issues.
// FMIndex data is stored in arenas owned by the FMIndex objects themselves.
class Index {
    index::MultiFMIndex fm_index_;
    std::vector<RefSequence> refs_;  // Heap-allocated, persistent
    std::string meta_;               // Heap-allocated, persistent
    std::vector<index::PackedSequence> packed_refs_;  // Packed reference sequences (persistent)
    std::vector<uint64_t> ref_concat_;  // Concatenated reference as packed bits (persistent)

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

    // Extract a packed reference subsequence [start, end) as uint8_t vector
    [[nodiscard]] std::vector<uint8_t> extract_ref(size_t start, size_t end) const {
        std::vector<uint8_t> result;
        if (packed_refs_.empty() || start >= end) return result;

        // Find which reference contains this range
        for (size_t ri = 0; ri < refs_.size(); ++ri) {
            const auto& ref = refs_[ri];
            if (start >= ref.offset && start < ref.offset + ref.length) {
                size_t local_start = start - ref.offset;
                size_t local_end = std::min(end - ref.offset, ref.length);
                const auto& packed = packed_refs_[ri];
                result.reserve(local_end - local_start);
                for (size_t i = local_start; i < local_end; ++i) {
                    result.push_back(packed.get(i));
                }
                return result;
            }
        }
        return result;
    }

    // Get total packed reference length
    [[nodiscard]] size_t ref_length() const noexcept {
        size_t total = 0;
        for (const auto& ref : refs_) total += ref.length;
        return total;
    }

private:
    void build_impl(const char* fasta_path, const Config& cfg);
    void load_impl(const char* prefix);
    void save_impl(const char* prefix) const;
};

// Alignment result for a single read
struct AlignmentResult {
    AlnRecord primary;
    std::vector<AlnRecord> secondary;
    std::vector<AlnRecord> supplementary;
    bool mapped = false;
    int32_t best_score = 0;
    int32_t second_best_score = 0;
    int32_t suboptimal_score = 0;  // XS tag: suboptimal alignment score

    void clear() noexcept {
        primary.clear(); secondary.clear(); supplementary.clear();
        mapped = false; best_score = second_best_score = suboptimal_score = 0;
    }
};

// Main aligner class
class Aligner {
    const Index& index_;
    Config config_;
    align::MEMFinder mem_finder_;
    mutable memory::Arena arena_;
    InsertSizeStats insert_stats_;

public:
    Aligner(const Index& idx, const Config& cfg = Config::default_mem())
        : index_(idx), config_(cfg),
          mem_finder_(idx.num_references() > 0 ? idx.fm_index()[0] : empty_fm_index(),
                      cfg.min_seed_len, cfg.max_occ, cfg.scoring.match),
          arena_(64 * 1024) {}

    // Insert-size stats (set from a sample before the main paired pass).
    void set_insert_stats(const InsertSizeStats& s) noexcept { insert_stats_ = s; }
    [[nodiscard]] const InsertSizeStats& insert_stats() const noexcept {
        return insert_stats_;
    }

    // Align single read
    AlignmentResult align(const io::SeqRecord& read) const {
        AlignmentResult result;
        align_impl(read, result);
        return result;
    }

    // Align paired-end reads (mate rescue + proper-pair marking).
    // Returns one result per end; unmapped ends are returned with
    // mapped=false but a populated unmapped primary record for SAM output.
    std::pair<AlignmentResult, AlignmentResult> align_pair(
        const io::SeqRecord& read1, const io::SeqRecord& read2) const {
        AlignmentResult res1, res2;
        align_pair_impl(read1, read2, res1, res2);
        return {std::move(res1), std::move(res2)};
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
    // Fallback FM-index for empty references (avoids UB from fm_index()[0])
    static const index::FMIndex& empty_fm_index() {
        static const index::FMIndex empty;
        return empty;
    }

    void align_impl(const io::SeqRecord& read, AlignmentResult& result) const;
    void align_pair_impl(const io::SeqRecord& read1,
                         const io::SeqRecord& read2,
                         AlignmentResult& res1,
                         AlignmentResult& res2) const;
    // Targeted SW for an unmapped end in a window around its mapped mate.
    bool rescue_end(const io::SeqRecord& read, const AlnRecord& mate,
                    AlignmentResult& result) const;
    // Set PAIRED/READ1/READ2/PROPER/UNMAP/MUNMAP + RNEXT/PNEXT/TLEN.
    void mark_pair(AlignmentResult& res, const io::SeqRecord& read,
                   bool is_read1, const AlnRecord* mate_primary,
                   bool mate_mapped, bool proper, int32_t tlen) const;

    // Convert MEM chain to alignment
    void chain_to_alignment(const align::MEMFinder::Chain& chain,
                            const io::SeqRecord& read,
                            const RefSequence& ref,
                            AlnRecord& aln) const;

public:
    // Input signals for mapping-quality estimation, mirroring the fields of
    // BWA's mem_alnreg_t used by mem_approx_mapq_se (bwamem.c).
    struct MapqSignals {
        int32_t score = 0;       // best alignment score (a->score)
        int32_t sub = 0;         // best suboptimal alignment score (a->sub), 0 if none
        int32_t csub = 0;        // best suboptimal chain score (a->csub), 0 if none
        int32_t sub_n = 0;       // # suboptimal hits (a->sub_n)
        int32_t seedcov = 0;     // seed coverage of primary (a->seedcov)
        int32_t query_span = 0;  // aligned query length (a->qe - a->qb)
        int32_t ref_span = 0;    // aligned ref length (a->re - a->rb)
        int32_t match = 1;       // match score (opt->a)
        int32_t mismatch = 4;    // mismatch penalty magnitude (opt->b)
        int32_t min_seed_len = 19;
    };

    // Detect and create supplementary alignments for chimeric reads
    void detect_supplementary(const io::SeqRecord& read,
                              const bwa::core::Vector<bwa::align::MEMFinder::Chain>& chains,
                              const bwa::memory::Arena& arena,
                              AlignmentResult& result) const;

    // Create supplementary alignment for a split-read segment. Takes the
    // oriented CIGAR (hard clips included by the caller), reference name,
    // 1-based leftmost position, strand and SW score.
    AlnRecord create_supplementary(const io::SeqRecord& read,
                                   const std::vector<uint32_t>& cigar,
                                   const std::string& rname,
                                   int32_t pos_1based,
                                   bool is_reverse,
                                   int32_t score) const;

    // Compute SA tag entry for a single alignment (RNAME,POS,STRAND,CIGAR,MAPQ,NM;)
    std::string make_sa_tag(const AlnRecord& aln) const;

    // Faithful port of BWA's mem_approx_mapq_se(): 6.02-scaled score gap
    // discounted by alignment length (mapQ_coef_len/fac) and identity, minus
    // 4.343*log(sub_n+1) per extra hit, clamped to [0, 60]. Public static so
    // unit tests can pin the model. Deviations from BWA, all documented at
    // the call sites: frac_rep scaling is skipped (repetitiveness untracked;
    // sub_n/second-best already depress repeats), and sub_n counts distinct
    // emitted chains rather than query-overlapping hits.
    static uint8_t approx_mapq_se(const MapqSignals& s) noexcept;

private:
    // Query-span union length of a chain's MEMs (seed coverage estimate)
    static int32_t seed_query_coverage(const align::MEMFinder::Chain& chain);
};

// High-level pipeline
class Pipeline {
    Index index_;
    Config config_;
    mutable Aligner aligner_;
    mutable InsertSizeStats insert_stats_;

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
        if (config_.num_threads <= 1) {
            align_file_serial(fastq_path, sam_path);
        } else {
            align_file_parallel(fastq_path, sam_path, config_.num_threads);
        }
    }

    // Serial alignment (single-threaded)
    void align_file_serial(const char* fastq_path, const char* sam_path) const {
        io::SeqReader reader(fastq_path);
        std::ofstream sam_file;
        std::ostream* out = &std::cout;

        if (std::string_view(sam_path) != "-") {
            // Check if output should be gzipped (.gz extension)
            std::string_view path(sam_path);
            if (path.size() >= 3 && path.substr(path.size()-3) == ".gz") {
                // BGZF support: write to a buffer, then compress
                std::stringstream buffer;
                write_header(buffer);
                aligner_.align_stream(reader, [&](const AlignmentResult& result) {
                    write_alignment(buffer, result);
                });
                write_gzipped(sam_path, buffer.str());
                return;
            }
            sam_file.open(sam_path, std::ios::binary);
            if (!sam_file) {
                throw std::runtime_error("Cannot open SAM output file");
            }
            out = &sam_file;
        }

        write_header(*out);

        aligner_.align_stream(reader, [&](const AlignmentResult& result) {
            write_alignment(*out, result);
        });

        if (sam_file.is_open()) sam_file.close();
    }

    // Write string to a gzip file
    void write_gzipped(const char* path, const std::string& data) const {
        gzFile gz = gzopen(path, "wb");
        if (!gz) throw std::runtime_error("Cannot open gzip output file");
        gzwrite(gz, data.data(), static_cast<unsigned>(data.size()));
        gzclose(gz);
    }

    // Align FASTQ file to BAM output (BGZF BAM + coordinate-sorted BAI).
    // Records are buffered, coordinate-sorted, then written so the BAI is
    // valid; memory use grows with the number of alignments.
    void align_to_bam(const char* fastq_path, const char* bam_path) const;

    // Paired-end variant of align_to_bam.
    void align_pair_to_bam(const char* fastq1, const char* fastq2,
                           const char* bam_path) const;

    // Parallel alignment using a bounded chunk pipeline: reads are consumed
    // in chunks of batch_size, each chunk is aligned by a work-stealing pool
    // (per-thread Aligner + TLS arena), and per-read SAM strings are
    // reassembled in read order. Output order is therefore identical to the
    // serial path. (.gz SAM falls back to the serial path.)
    void align_file_parallel(const char* fastq_path, const char* sam_path, int num_threads) const {
        io::SeqReader reader(fastq_path);
        std::ofstream sam_file;
        std::ostream* out = &std::cout;

        if (std::string_view(sam_path) != "-") {
            std::string_view path(sam_path);
            if (path.size() >= 3 && path.substr(path.size()-3) == ".gz") {
                // For BGZF, buffer all output then compress
                std::stringstream buffer;
                write_header(buffer);
                aligner_.align_stream(reader, [&](const AlignmentResult& result) {
                    write_alignment(buffer, result);
                });
                write_gzipped(sam_path, buffer.str());
                return;
            }
            sam_file.open(sam_path, std::ios::binary);
            if (!sam_file) {
                throw std::runtime_error("Cannot open SAM output file");
            }
            out = &sam_file;
        }

        write_header(*out);

        const size_t chunk_size =
            std::max<size_t>(static_cast<size_t>(config_.batch_size), 64);
        std::vector<io::SeqRecord> reads(chunk_size);
        std::vector<std::string> pieces(chunk_size);
        const int T = std::max(1, num_threads);

        for (;;) {
            size_t n = 0;
            while (n < chunk_size && reader.read(reads[n])) ++n;
            if (n == 0) break;

            std::atomic<size_t> next{0};
            auto worker = [&]() {
                Aligner local_aligner(index_, config_);
                std::ostringstream os;
                for (;;) {
                    size_t i = next.fetch_add(1, std::memory_order_relaxed);
                    if (i >= n) break;
                    AlignmentResult result = local_aligner.align(reads[i]);
                    os.str(std::string());
                    os.clear();
                    write_alignment(os, result);
                    pieces[i] = os.str();
                    memory::reset_tls_arena();
                }
            };
            std::vector<std::thread> threads;
            const int nt = static_cast<int>(std::min<size_t>(T, n));
            for (int t = 0; t < nt; ++t) threads.emplace_back(worker);
            for (auto& t : threads) t.join();

            for (size_t i = 0; i < n; ++i) {
                *out << pieces[i];
                pieces[i].clear();
            }
            // Reads were parsed into this thread's arena; release them before
            // parsing the next chunk.
            memory::reset_tls_arena();
        }

        if (sam_file.is_open()) sam_file.close();
    }

    // Parallel paired-end alignment (same chunked scheme as above).
    void align_pair_parallel(const char* fastq1, const char* fastq2,
                             const char* sam_path, int num_threads) const {
        io::SeqReader r1(fastq1), r2(fastq2);
        std::ofstream sam_file;
        std::ostream* out = &std::cout;

        if (std::string_view(sam_path) != "-") {
            std::string_view path(sam_path);
            if (path.size() >= 3 && path.substr(path.size()-3) == ".gz") {
                std::stringstream buffer;
                write_header(buffer);
                io::SeqRecord read1, read2;
                while (r1.read(read1) && r2.read(read2)) {
                    auto [res1, res2] = aligner_.align_pair(read1, read2);
                    write_pair(buffer, res1, res2);
                    memory::reset_tls_arena();
                }
                write_gzipped(sam_path, buffer.str());
                return;
            }
            sam_file.open(sam_path, std::ios::binary);
            if (!sam_file) {
                throw std::runtime_error("Cannot open SAM output file");
            }
            out = &sam_file;
        }

        write_header(*out);

        const size_t chunk_size =
            std::max<size_t>(static_cast<size_t>(config_.batch_size), 64);
        std::vector<io::SeqRecord> reads1(chunk_size), reads2(chunk_size);
        std::vector<std::string> pieces(chunk_size);
        const int T = std::max(1, num_threads);

        for (;;) {
            size_t n = 0;
            while (n < chunk_size && r1.read(reads1[n]) && r2.read(reads2[n])) ++n;
            if (n == 0) break;

            std::atomic<size_t> next{0};
            auto worker = [&]() {
                Aligner local_aligner(index_, config_);
                std::ostringstream os;
                for (;;) {
                    size_t i = next.fetch_add(1, std::memory_order_relaxed);
                    if (i >= n) break;
                    auto [res1, res2] =
                        local_aligner.align_pair(reads1[i], reads2[i]);
                    os.str(std::string());
                    os.clear();
                    write_pair(os, res1, res2);
                    pieces[i] = os.str();
                    memory::reset_tls_arena();
                }
            };
            std::vector<std::thread> threads;
            const int nt = static_cast<int>(std::min<size_t>(T, n));
            for (int t = 0; t < nt; ++t) threads.emplace_back(worker);
            for (auto& t : threads) t.join();

            for (size_t i = 0; i < n; ++i) {
                *out << pieces[i];
                pieces[i].clear();
            }
            // Reads were parsed into this thread's arena; release them before
            // parsing the next chunk.
            memory::reset_tls_arena();
        }

        if (sam_file.is_open()) sam_file.close();
    }

    // Align paired FASTQ files
    void align_pair(const char* fastq1, const char* fastq2, const char* sam_path = "-") const {
        ensure_insert_size(fastq1, fastq2);
        if (config_.num_threads > 1) {
            align_pair_parallel(fastq1, fastq2, sam_path, config_.num_threads);
            return;
        }
        io::SeqReader r1(fastq1), r2(fastq2);
        std::ofstream sam_file;
        std::ostream* out = &std::cout;

        if (std::string_view(sam_path) != "-") {
            std::string_view path(sam_path);
            if (path.size() >= 3 && path.substr(path.size()-3) == ".gz") {
                // For BGZF, buffer all output then compress
                std::stringstream buffer;
                write_header(buffer);
                io::SeqRecord read1, read2;
                while (r1.read(read1) && r2.read(read2)) {
                    auto [res1, res2] = aligner_.align_pair(read1, read2);
                    write_pair(buffer, res1, res2);
                    memory::reset_tls_arena();
                }
                write_gzipped(sam_path, buffer.str());
                return;
            }
            sam_file.open(sam_path, std::ios::binary);
            if (!sam_file) {
                throw std::runtime_error("Cannot open SAM output file");
            }
            out = &sam_file;
        }

        write_header(*out);

        io::SeqRecord read1, read2;
        while (r1.read(read1) && r2.read(read2)) {
            auto [res1, res2] = aligner_.align_pair(read1, read2);
            write_pair(*out, res1, res2);
            memory::reset_tls_arena();
        }

        if (sam_file.is_open()) sam_file.close();
    }

    // Align interleaved paired FASTQ (-p): consecutive records form pairs.
    void align_pair_interleaved(const char* fastq, const char* sam_path = "-") const {
        ensure_insert_size(fastq, nullptr);
        io::SeqReader reader(fastq);
        std::ofstream sam_file;
        std::ostream* out = &std::cout;

        if (std::string_view(sam_path) != "-") {
            sam_file.open(sam_path, std::ios::binary);
            if (!sam_file) {
                throw std::runtime_error("Cannot open SAM output file");
            }
            out = &sam_file;
        }

        write_header(*out);

        io::SeqRecord read1, read2;
        while (reader.read(read1) && reader.read(read2)) {
            auto [res1, res2] = aligner_.align_pair(read1, read2);
            write_pair(*out, res1, res2);
            memory::reset_tls_arena();
        }

        if (sam_file.is_open()) sam_file.close();
    }

    // Interleaved paired-end variant of align_pair_to_bam.
    void align_pair_interleaved_bam(const char* fastq, const char* bam_path) const;

private:
    // Estimate the insert-size distribution from a sample of read pairs
    // (single-threaded pre-pass; no-op for stdin or when already valid).
    // A null fastq2 means interleaved input (consecutive pairs in fastq1).
    void ensure_insert_size(const char* fastq1, const char* fastq2) const;

    void write_header(std::ostream& out) const;
    void write_alignment(std::ostream& out, const AlignmentResult& result) const;
    void write_pair(std::ostream& out, const AlignmentResult& res1,
                    const AlignmentResult& res2) const;
    void write_sam_record(std::ostream& out, const AlnRecord& aln) const;
};

} // namespace bwa