// BAI (BAM Index) writer
// BAI format: http://samtools.github.io/hts-specs/BAIv1.pdf
#pragma once

#include <cstdint>
#include <vector>

namespace bwa::io {

// BAI (BAM Index) writer. Records are expected to be added in coordinate
// order (same order as the BAM records they index).
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

    // Initialize with reference lengths (from the BAM header @SQ lines).
    void init(const std::vector<int64_t>& ref_lengths);

    // Record an alignment for BAI indexing. pos/end_pos are 0-based
    // half-open reference coordinates; file_offset is the virtual offset of
    // the record start. Records with ref_idx < 0 or pos < 0 count as
    // unplaced unmapped reads.
    void record_alignment(int32_t ref_idx, int32_t pos, int32_t end_pos,
                          int64_t file_offset);

    // Write the index. eof_offset is the virtual offset just past the last
    // BAM record (used as the end of the final chunk).
    void write_index(int64_t eof_offset);
};

} // namespace bwa::io
