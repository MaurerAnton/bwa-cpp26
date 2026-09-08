// BAM (Binary Alignment/Map) format writer
// BAM is a binary representation of SAM, compressed with BGZF
// Reference: https://samtools.github.io/hts-specs/SAMv1.pdf
#pragma once

#include <zlib.h>
#include <cstdint>
#include <cstring>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace bwa::io {

// Forward declarations to avoid including pipeline.hpp
class Pipeline;
class Config;
class Index;
class Aligner;
struct AlignmentResult;
struct AlnRecord;

// Helper to extract fields from AlnRecord without including pipeline.hpp
struct AlnRecordView {
    std::string qname;
    uint32_t flag;
    std::string rname;
    int32_t pos;
    uint8_t mapq;
    std::vector<uint32_t> cigar;
    std::string rnext;
    int32_t pnext;
    int32_t tlen;
    std::string seq;
    std::string qual;
    std::vector<std::pair<std::string, std::string>> tags;
    int32_t score;
};

// BGZF block constants
constexpr int BGZF_BLOCK_SIZE = 65280;  // 64KB minus headers
constexpr int BGZF_MAX_BLOCK_SIZE = 65536;
constexpr int BGZF_HEADER_SIZE = 18;
constexpr int BGZF_FOOTER_SIZE = 8;

// BAM magic bytes
constexpr uint8_t BAM_MAGIC[4] = {'B', 'A', 'M', 1};

// Variable-length integer encoding (BGZF/BAM)
// Writes an unsigned int in variable-length format
inline void write_itf8(std::vector<uint8_t>& buf, int32_t val) {
    if (val >= 0 && val < 0x80) {
        buf.push_back(static_cast<uint8_t>(val));
    } else if (val >= 0 && val < 0x4000) {
        buf.push_back(static_cast<uint8_t>((val >> 8) | 0x80));
        buf.push_back(static_cast<uint8_t>(val & 0xFF));
    } else if (val >= 0 && val < 0x200000) {
        buf.push_back(static_cast<uint8_t>((val >> 16) | 0xC0));
        buf.push_back(static_cast<uint8_t>((val >> 8) & 0xFF));
        buf.push_back(static_cast<uint8_t>(val & 0xFF));
    } else if (val >= 0 && val < 0x10000000) {
        buf.push_back(static_cast<uint8_t>((val >> 24) | 0xE0));
        buf.push_back(static_cast<uint8_t>((val >> 16) & 0xFF));
        buf.push_back(static_cast<uint8_t>((val >> 8) & 0xFF));
        buf.push_back(static_cast<uint8_t>(val & 0xFF));
    } else {
        buf.push_back(0xF0 | ((val >> 28) & 0x0F));
        buf.push_back(static_cast<uint8_t>((val >> 20) & 0xFF));
        buf.push_back(static_cast<uint8_t>((val >> 12) & 0xFF));
        buf.push_back(static_cast<uint8_t>((val >> 4) & 0xFF));
        buf.push_back(static_cast<uint8_t>(val & 0x0F));
    }
}

// Write a string in BAM format (length-prefixed)
inline void write_bam_string(std::vector<uint8_t>& buf, const std::string& s) {
    write_itf8(buf, static_cast<int32_t>(s.size()));
    buf.insert(buf.end(), s.begin(), s.end());
}

// Write little-endian int32
inline void write_le32(std::vector<uint8_t>& buf, int32_t val) {
    buf.push_back(val & 0xFF);
    buf.push_back((val >> 8) & 0xFF);
    buf.push_back((val >> 16) & 0xFF);
    buf.push_back((val >> 24) & 0xFF);
}

// Write little-endian uint32
inline void write_le_u32(std::vector<uint8_t>& buf, uint32_t val) {
    buf.push_back(val & 0xFF);
    buf.push_back((val >> 8) & 0xFF);
    buf.push_back((val >> 16) & 0xFF);
    buf.push_back((val >> 24) & 0xFF);
}

// Write little-endian int64
inline void write_le64(std::vector<uint8_t>& buf, int64_t val) {
    for (int i = 0; i < 8; ++i) {
        buf.push_back((val >> (i*8)) & 0xFF);
    }
}

// Write little-endian uint8
inline void write_u8(std::vector<uint8_t>& buf, uint8_t val) {
    buf.push_back(val);
}

// BGZF writer
class BgzfWriter {
    gzFile gz_ = nullptr;
    std::vector<uint8_t> block_;  // Current uncompressed block
    size_t block_pos_ = 0;        // Position in current block
    int64_t block_offset_ = 0;    // Virtual offset of current block
    int64_t total_offset_ = 0;    // Total bytes written (compressed)
    int64_t uncompressed_offset_ = 0; // Uncompressed byte offset

public:
    BgzfWriter() = default;
    ~BgzfWriter() { close(); }

    BgzfWriter(const BgzfWriter&) = delete;
    BgzfWriter& operator=(const BgzfWriter&) = delete;

    [[nodiscard]] bool open(const char* path) {
        gz_ = gzopen(path, "wb");
        if (!gz_) return false;
        block_.resize(BGZF_MAX_BLOCK_SIZE);
        return true;
    }

    void close() {
        if (gz_) {
            flush_block();
            gzclose(gz_);
            gz_ = nullptr;
        }
    }

    // Write raw bytes
    void write(const uint8_t* data, size_t len) {
        while (len > 0) {
            if (block_pos_ >= BGZF_BLOCK_SIZE) {
                flush_block();
            }
            size_t to_write = (len < BGZF_BLOCK_SIZE - block_pos_) ? len : (BGZF_BLOCK_SIZE - block_pos_);
            std::memcpy(block_.data() + block_pos_, data, to_write);
            block_pos_ += to_write;
            data += to_write;
            len -= to_write;
            uncompressed_offset_ += to_write;
        }
    }

    // Write a vector
    void write(const std::vector<uint8_t>& data) {
        write(data.data(), data.size());
    }

    // Flush current block
    void flush_block() {
        if (block_pos_ == 0) return;

        // Compress with raw deflate
        uLongf compressed_size = BGZF_MAX_BLOCK_SIZE;
        std::vector<uint8_t> compressed(BGZF_MAX_BLOCK_SIZE);
        compress2(compressed.data(), &compressed_size,
                  block_.data(), block_pos_, Z_DEFAULT_COMPRESSION);

        // Write BGZF header (18 bytes)
        std::vector<uint8_t> header;
        header.push_back(0x1F);  // gzip magic 1
        header.push_back(0x8B);  // gzip magic 2
        header.push_back(0x08);  // compression method (deflate)
        header.push_back(0x04);  // FLG (FEXTRA set)
        header.push_back(0x00);  // MTIME
        header.push_back(0x00);
        header.push_back(0x00);
        header.push_back(0x00);
        header.push_back(0xFF);  // XFL (unknown)
        header.push_back(0x06);  // OS (unknown)
        header.push_back(0x00);  // XLEN
        header.push_back(0x00);
        // Extra field: BC subfield
        header.push_back(0x42);  // SI1 = 'B'
        header.push_back(0x43);  // SI2 = 'C'
        header.push_back(0x02);  // SLEN
        header.push_back(0x00);
        // BSIZE = total block size - 1
        uint16_t bsize = static_cast<uint16_t>(header.size() + 4 + compressed_size - 1);
        header.push_back(bsize & 0xFF);
        header.push_back((bsize >> 8) & 0xFF);

        // Write header + compressed data
        gzwrite(gz_, header.data(), static_cast<unsigned>(header.size()));
        gzwrite(gz_, compressed.data(), static_cast<unsigned>(compressed_size));

        // Footer (CRC32 + ISIZE)
        uint32_t crc = crc32(0L, block_.data(), static_cast<uInt>(block_pos_));
        std::vector<uint8_t> footer;
        footer.push_back(crc & 0xFF);
        footer.push_back((crc >> 8) & 0xFF);
        footer.push_back((crc >> 16) & 0xFF);
        footer.push_back((crc >> 24) & 0xFF);
        uint32_t isize = static_cast<uint32_t>(block_pos_);
        footer.push_back(isize & 0xFF);
        footer.push_back((isize >> 8) & 0xFF);
        footer.push_back((isize >> 16) & 0xFF);
        footer.push_back((isize >> 24) & 0xFF);
        gzwrite(gz_, footer.data(), static_cast<unsigned>(footer.size()));

        // Update offsets
        total_offset_ += header.size() + compressed_size + 8;
        block_offset_ = total_offset_;
        block_pos_ = 0;
    }

    [[nodiscard]] int64_t virtual_offset() const {
        return (block_offset_ << 16) | block_pos_;
    }

    [[nodiscard]] int64_t total_offset() const { return total_offset_; }
};

// BAM writer
class BamWriter {
    BgzfWriter bgzf_;
    std::vector<int32_t> ref_seqs_;  // Index of each reference

public:
    BamWriter() = default;
    ~BamWriter() { close(); }

    [[nodiscard]] bool open(const char* path) {
        return bgzf_.open(path);
    }

    void close() {
        bgzf_.close();
    }

    [[nodiscard]] int64_t virtual_offset() const { return bgzf_.virtual_offset(); }
    [[nodiscard]] int64_t total_offset() const { return bgzf_.total_offset(); }

    // Write BAM header
    void write_header(const std::string& sam_header, int32_t num_refs) {
        std::vector<uint8_t> buf;

        // Magic
        for (int i = 0; i < 4; ++i) buf.push_back(BAM_MAGIC[i]);

        // l_text (length of header text)
        write_itf8(buf, static_cast<int32_t>(sam_header.size()));
        buf.insert(buf.end(), sam_header.begin(), sam_header.end());

        // n_ref (number of reference sequences)
        write_le32(buf, num_refs);

        // For each reference: l_name, name, l_ref
        // This is written separately when we have the reference info
        // For now, just write the header structure

        bgzf_.write(buf);
    }

    // Write a reference sequence entry
    void write_reference(int32_t ref_len, const std::string& ref_name) {
        std::vector<uint8_t> buf;
        write_bam_string(buf, ref_name);
        write_le32(buf, ref_len);
        bgzf_.write(buf);
    }

    // Encode a CIGAR operation
    static uint32_t encode_cigar_op(int32_t len, char op) {
        return (static_cast<uint32_t>(len) << 4) | (op & 0xF);
    }

    // Write an alignment record
    void write_alignment(const AlnRecordView& aln, int32_t ref_idx) {
        std::vector<uint8_t> buf;

        // block_size (placeholder, will be filled)
        size_t block_size_pos = buf.size();
        write_le32(buf, 0);

        int32_t start = buf.size();

        // refID
        write_le32(buf, ref_idx);

        // pos (0-based)
        write_le32(buf, aln.pos > 0 ? aln.pos - 1 : -1);

        // l_read_name
        write_u8(buf, static_cast<uint8_t>(aln.qname.size() + 1));

        // mapq
        write_u8(buf, aln.mapq);

        // bin
        // Calculate bin from alignment position
        int32_t pos = aln.pos > 0 ? aln.pos - 1 : 0;
        uint16_t bin = 4680; // Maximum bin
        if (pos < 0) {
            bin = 0;
        } else {
            // Use a simple binning scheme
            if (pos < 65536) bin = pos >> 14;
            else if (pos < 262144) bin = 9 + (pos >> 17);
            else if (pos < 1048576) bin = 15 + (pos >> 20);
            else if (pos < 4194304) bin = 21 + (pos >> 23);
            else if (pos < 16777216) bin = 27 + (pos >> 26);
            else bin = 4680;
        }
        buf.push_back(bin & 0xFF);
        buf.push_back((bin >> 8) & 0xFF);

        // n_cigar_op
        uint16_t n_cigar = static_cast<uint16_t>(aln.cigar.size());
        buf.push_back(n_cigar & 0xFF);
        buf.push_back((n_cigar >> 8) & 0xFF);

        // flag
        write_le_u32(buf, aln.flag);

        // l_seq
        int32_t l_seq = static_cast<int32_t>(aln.seq.size());
        write_le32(buf, l_seq);

        // next_refID
        write_le32(buf, aln.rnext == "*" ? -1 : ref_idx);

        // next_pos
        write_le32(buf, aln.pnext > 0 ? aln.pnext - 1 : -1);

        // tlen
        write_le32(buf, aln.tlen);

        // read_name (null-terminated)
        buf.insert(buf.end(), aln.qname.begin(), aln.qname.end());
        buf.push_back(0);

        // cigar
        for (uint32_t c : aln.cigar) {
            write_le_u32(buf, c);
        }

        // seq (4-bit packed, 0=ACGTN map 0,1,2,3,4)
        for (size_t i = 0; i < aln.seq.size(); i += 2) {
            uint8_t b1 = 0, b2 = 0;
            char c1 = aln.seq[i];
            char c2 = (i + 1 < aln.seq.size()) ? aln.seq[i + 1] : 0;
            // Map ACGT -> 0,1,2,3; N -> 4
            switch (c1) {
                case 'A': case 'a': b1 = 0; break;
                case 'C': case 'c': b1 = 1; break;
                case 'G': case 'g': b1 = 2; break;
                case 'T': case 't': b1 = 3; break;
                default: b1 = 4; break;
            }
            switch (c2) {
                case 'A': case 'a': b2 = 0; break;
                case 'C': case 'c': b2 = 1; break;
                case 'G': case 'g': b2 = 2; break;
                case 'T': case 't': b2 = 3; break;
                default: b2 = 0; break;
            }
            buf.push_back(b1 << 4 | b2);
        }

        // qual
        for (size_t i = 0; i < aln.qual.size(); ++i) {
            buf.push_back(static_cast<uint8_t>(aln.qual[i] - 33));
        }

        // aux data (tags)
        for (const auto& tag : aln.tags) {
            // Tag: 2 char tag name + type + value
            if (tag.first.size() < 2) continue;
            buf.push_back(tag.first[0]);
            buf.push_back(tag.first[1]);
            char type = tag.second.size() == 1 ? 'A' :
                       tag.second[0] == '-' ? 'i' : 'Z';
            buf.push_back(type);
            if (type == 'i') {
                int32_t val = std::stoi(tag.second);
                write_le32(buf, val);
            } else {
                buf.insert(buf.end(), tag.second.begin(), tag.second.end());
                buf.push_back(0);
            }
        }

        // Fill in block_size
        int32_t block_size = static_cast<int32_t>(buf.size() - start);
        buf[block_size_pos] = block_size & 0xFF;
        buf[block_size_pos + 1] = (block_size >> 8) & 0xFF;
        buf[block_size_pos + 2] = (block_size >> 16) & 0xFF;
        buf[block_size_pos + 3] = (block_size >> 24) & 0xFF;

        bgzf_.write(buf);
    }

}; // class BamWriter

} // namespace bwa::io
