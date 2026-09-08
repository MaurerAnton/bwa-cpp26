// BAI (BAM Index) writer
// BAI format: http://samtools.github.io/hts-specs/BAIv1.pdf
#pragma once

#include <cstdint>

namespace bwa::io {

// BAI (BAM Index) writer
// BAI format: http://samtools.github.io/hts-specs/BAIv1.pdf
class BaiWriter {
    struct Impl;
    Impl* impl_;

public:
    BaiWriter();
    ~BaiWriter();
    BaiWriter(BaiWriter&&) noexcept;
    BaiWriter& operator=(BaiWriter&&) noexcept;
    BaiWriter(const BaiWriter&) = delete;
    BaiWriter& operator=(const BaiWriter&) = delete;

    [[nodiscard]] bool open(const char* path);
    void close();

    // Initialize with number of reference sequences
    void init(int32_t num_refs);

    // Record an alignment for BAI indexing
    void record_alignment(int32_t ref_idx, int32_t pos, int64_t file_offset);

    // Write the BAI index file
    void write_index();
};

} // namespace bwa::io