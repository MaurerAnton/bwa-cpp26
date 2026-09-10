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
#include <algorithm>

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
    int32_t next_ref_idx = -1;  // resolved mate reference for BAM next_refID
};

// BAM/BAI binning (SAM spec reg2bin). beg is 0-based inclusive, end 0-based
// exclusive. Returns the smallest bin fully containing [beg, end).
[[nodiscard]] inline uint16_t reg2bin(int32_t beg, int32_t end) noexcept {
    if (beg < 0) beg = 0;
    if (end <= beg) end = beg + 1;
    --end;
    if ((beg >> 14) == (end >> 14)) return static_cast<uint16_t>(4681 + (beg >> 14));
    if ((beg >> 17) == (end >> 17)) return static_cast<uint16_t>(585 + (beg >> 17));
    if ((beg >> 20) == (end >> 20)) return static_cast<uint16_t>(73 + (beg >> 20));
    if ((beg >> 23) == (end >> 23)) return static_cast<uint16_t>(9 + (beg >> 23));
    if ((beg >> 26) == (end >> 26)) return static_cast<uint16_t>(1 + (beg >> 26));
    return 0;
}

// BGZF block constants
constexpr int BGZF_BLOCK_SIZE = 65280;  // 64KB minus headers, keeps worst case < 65536
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

// BGZF writer. Writes raw BGZF blocks (gzip members with a BC extra field)
// to a plain binary stream. Virtual offsets are (compressed block offset <<
// 16) | uncompressed offset within the block, per the BAM spec.
class BgzfWriter {
    std::ofstream out_;
    std::vector<uint8_t> block_;  // Current uncompressed block
    size_t block_pos_ = 0;        // Bytes used in current block
    int64_t block_offset_ = 0;    // Compressed offset of current block
    int64_t total_offset_ = 0;    // Compressed bytes written so far
    bool open_ = false;

public:
    BgzfWriter() = default;
    ~BgzfWriter() { close(); }

    BgzfWriter(const BgzfWriter&) = delete;
    BgzfWriter& operator=(const BgzfWriter&) = delete;

    [[nodiscard]] bool open(const char* path) {
        out_.open(path, std::ios::binary);
        if (!out_) return false;
        block_.resize(BGZF_BLOCK_SIZE);
        block_pos_ = 0;
        block_offset_ = 0;
        total_offset_ = 0;
        open_ = true;
        return true;
    }

    void close() {
        if (!open_) return;
        flush_block();
        write_eof_block();
        out_.close();
        open_ = false;
    }

    // Write raw bytes
    void write(const uint8_t* data, size_t len) {
        while (len > 0) {
            if (block_pos_ >= static_cast<size_t>(BGZF_BLOCK_SIZE)) flush_block();
            size_t to_write = std::min(len, static_cast<size_t>(BGZF_BLOCK_SIZE) - block_pos_);
            std::memcpy(block_.data() + block_pos_, data, to_write);
            block_pos_ += to_write;
            data += to_write;
            len -= to_write;
        }
    }

    // Write a vector
    void write(const std::vector<uint8_t>& data) {
        write(data.data(), data.size());
    }

    // Flush current block as one BGZF member.
    void flush_block() {
        if (block_pos_ == 0) return;

        // BGZF payload is RAW deflate (windowBits = -15), not zlib-wrapped.
        std::vector<uint8_t> compressed(compressBound(static_cast<uLong>(block_pos_)));
        z_stream strm{};
        deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8,
                     Z_DEFAULT_STRATEGY);
        strm.next_in = block_.data();
        strm.avail_in = static_cast<uInt>(block_pos_);
        strm.next_out = compressed.data();
        strm.avail_out = static_cast<uInt>(compressed.size());
        deflate(&strm, Z_FINISH);
        uLongf compressed_size = strm.total_out;
        deflateEnd(&strm);

        // BSIZE = total member size - 1 = 18 (header) + compressed + 8 (footer) - 1
        uint16_t bsize = static_cast<uint16_t>(BGZF_HEADER_SIZE + compressed_size +
                                               BGZF_FOOTER_SIZE - 1);
        uint8_t header[BGZF_HEADER_SIZE] = {
            0x1F, 0x8B, 0x08, 0x04,             // gzip magic, CM=deflate, FLG=FEXTRA
            0x00, 0x00, 0x00, 0x00,             // MTIME
            0xFF, 0x06,                         // XFL, OS=unknown
            0x06, 0x00,                         // XLEN=6
            0x42, 0x43, 0x02, 0x00,             // BC subfield, SLEN=2
            static_cast<uint8_t>(bsize & 0xFF),
            static_cast<uint8_t>((bsize >> 8) & 0xFF)
        };
        uint32_t crc = crc32(0L, block_.data(), static_cast<uInt>(block_pos_));
        uint32_t isize = static_cast<uint32_t>(block_pos_);
        uint8_t footer[BGZF_FOOTER_SIZE] = {
            static_cast<uint8_t>(crc & 0xFF),
            static_cast<uint8_t>((crc >> 8) & 0xFF),
            static_cast<uint8_t>((crc >> 16) & 0xFF),
            static_cast<uint8_t>((crc >> 24) & 0xFF),
            static_cast<uint8_t>(isize & 0xFF),
            static_cast<uint8_t>((isize >> 8) & 0xFF),
            static_cast<uint8_t>((isize >> 16) & 0xFF),
            static_cast<uint8_t>((isize >> 24) & 0xFF)
        };

        out_.write(reinterpret_cast<const char*>(header), BGZF_HEADER_SIZE);
        out_.write(reinterpret_cast<const char*>(compressed.data()),
                   static_cast<std::streamsize>(compressed_size));
        out_.write(reinterpret_cast<const char*>(footer), BGZF_FOOTER_SIZE);

        total_offset_ += BGZF_HEADER_SIZE + static_cast<int64_t>(compressed_size) +
                         BGZF_FOOTER_SIZE;
        block_offset_ = total_offset_;
        block_pos_ = 0;
    }

    // Standard 28-byte BGZF EOF marker.
    void write_eof_block() {
        static const uint8_t eof[28] = {
            0x1F, 0x8B, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0xFF, 0x06,
            0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1B, 0x00, 0x03, 0x00,
            0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
        };
        out_.write(reinterpret_cast<const char*>(eof), 28);
        total_offset_ += 28;
    }

    [[nodiscard]] int64_t virtual_offset() const {
        return (block_offset_ << 16) | static_cast<int64_t>(block_pos_);
    }

    [[nodiscard]] int64_t total_offset() const { return total_offset_; }
    [[nodiscard]] bool is_open() const noexcept { return open_; }
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

        // l_text is a plain little-endian int32 (not ITF8) per the BAM spec.
        write_le32(buf, static_cast<int32_t>(sam_header.size()));
        buf.insert(buf.end(), sam_header.begin(), sam_header.end());

        // n_ref (number of reference sequences)
        write_le32(buf, num_refs);

        // Reference entries follow via write_reference().
        bgzf_.write(buf);
    }

    // Write a reference sequence entry: l_name (int32 LE, includes NUL),
    // NUL-terminated name, l_ref (int32 LE).
    void write_reference(int32_t ref_len, const std::string& ref_name) {
        std::vector<uint8_t> buf;
        write_le32(buf, static_cast<int32_t>(ref_name.size() + 1));
        buf.insert(buf.end(), ref_name.begin(), ref_name.end());
        buf.push_back(0);
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

        // bin: smallest bin containing [pos, end) per SAM spec. Unmapped
        // (pos < 0) records use bin 4680.
        int32_t pos = aln.pos > 0 ? aln.pos - 1 : -1;
        uint16_t bin = 4680;
        if (pos >= 0) {
            int32_t ref_span = 0;
            for (uint32_t c : aln.cigar) {
                char op = static_cast<char>(c & 0xF);
                int len = static_cast<int>(c >> 4);
                // M(0) =/X(7,8) D(2) N(3) consume reference
                if (op == 0 || op == 2 || op == 3 || op == 7 || op == 8) {
                    ref_span += len;
                }
            }
            bin = reg2bin(pos, pos + (ref_span > 0 ? ref_span : 1));
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

        // next_refID (resolved by the caller; "=" handled there)
        write_le32(buf, aln.next_ref_idx);

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
        auto base_code = [](char c) -> uint8_t {
            switch (c) {
                case 'A': case 'a': return 0;
                case 'C': case 'c': return 1;
                case 'G': case 'g': return 2;
                case 'T': case 't': return 3;
                case '=': return 0;
                default: return 4;
            }
        };
        for (size_t i = 0; i < aln.seq.size(); i += 2) {
            uint8_t b1 = base_code(aln.seq[i]);
            uint8_t b2 = (i + 1 < aln.seq.size()) ? base_code(aln.seq[i + 1]) : 0;
            buf.push_back(static_cast<uint8_t>((b1 << 4) | b2));
        }

        // qual (SAM ASCII-33; '*' means unavailable -> 0xFF)
        if (aln.qual == "*") {
            for (int32_t i = 0; i < l_seq; ++i) buf.push_back(0xFF);
        } else {
            for (size_t i = 0; i < aln.qual.size(); ++i) {
                buf.push_back(static_cast<uint8_t>(aln.qual[i] - 33));
            }
        }

        // aux data (tags). Keys carry their SAM type suffix ("NM:i", "MD:Z",
        // "SA:Z"); legacy bare keys default to Z.
        for (const auto& tag : aln.tags) {
            if (tag.first.size() < 2) continue;
            char type = 'Z';
            std::string val = tag.second;
            if (tag.first.size() >= 4 && tag.first[2] == ':') {
                type = tag.first[3];
            }
            buf.push_back(static_cast<uint8_t>(tag.first[0]));
            buf.push_back(static_cast<uint8_t>(tag.first[1]));
            buf.push_back(static_cast<uint8_t>(type));
            if (type == 'i') {
                write_le32(buf, static_cast<int32_t>(std::stol(val)));
            } else if (type == 'A') {
                buf.push_back(val.empty() ? static_cast<uint8_t>('?')
                                          : static_cast<uint8_t>(val[0]));
            } else if (type == 'f') {
                float f = std::stof(val);
                uint32_t bits;
                std::memcpy(&bits, &f, sizeof(bits));
                write_le_u32(buf, bits);
            } else {  // Z/H/B treated as string for our tag set
                buf.insert(buf.end(), val.begin(), val.end());
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
