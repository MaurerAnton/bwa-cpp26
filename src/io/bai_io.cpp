#include <bwa/io/bai_io.hpp>
#include <bwa/io/bam_io.hpp>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <cstdint>
#include <memory>

namespace bwa::io {

struct BaiWriter::Impl {
    struct Rec {
        int32_t ref = -1;
        int32_t pos = -1;
        int32_t end = -1;
        int64_t off = 0;
    };
    std::ofstream out_;
    std::vector<int64_t> ref_lengths_;
    std::vector<Rec> recs_;       // in file (coordinate) order
    int64_t n_no_coor_ = 0;       // unplaced unmapped reads
    bool open_ = false;
};

BaiWriter::BaiWriter() : impl_(new Impl()) {}
BaiWriter::~BaiWriter() { delete impl_; }
BaiWriter::BaiWriter(BaiWriter&& other) noexcept : impl_(other.impl_) { other.impl_ = nullptr; }
BaiWriter& BaiWriter::operator=(BaiWriter&& other) noexcept {
    if (this != &other) {
        delete impl_;
        impl_ = other.impl_;
        other.impl_ = nullptr;
    }
    return *this;
}

[[nodiscard]] bool BaiWriter::open(const char* path) {
    // BAI is an uncompressed binary file (not BGZF).
    impl_->out_.open(path, std::ios::binary);
    impl_->open_ = static_cast<bool>(impl_->out_);
    return impl_->open_;
}

void BaiWriter::close() {
    if (impl_->out_.is_open()) impl_->out_.close();
    impl_->open_ = false;
}

void BaiWriter::init(const std::vector<int64_t>& ref_lengths) {
    impl_->ref_lengths_ = ref_lengths;
    impl_->recs_.clear();
    impl_->n_no_coor_ = 0;
}

void BaiWriter::record_alignment(int32_t ref_idx, int32_t pos, int32_t end_pos,
                                 int64_t file_offset) {
    if (ref_idx < 0 || pos < 0) {
        impl_->n_no_coor_++;
        return;
    }
    if (ref_idx >= static_cast<int32_t>(impl_->ref_lengths_.size())) return;
    if (end_pos <= pos) end_pos = pos + 1;
    impl_->recs_.push_back({ref_idx, pos, end_pos, file_offset});
}

void BaiWriter::write_index(int64_t eof_offset) {
    std::vector<uint8_t> buf;
    buf.push_back('B');
    buf.push_back('A');
    buf.push_back('I');
    buf.push_back(1);

    const int32_t num_refs = static_cast<int32_t>(impl_->ref_lengths_.size());
    write_le32(buf, num_refs);

    // Group records per reference and bin. Records are coordinate-sorted, so
    // a bin's records form contiguous file-order runs.
    struct Chunk {
        int64_t beg = 0;
        int64_t end = 0;
    };
    std::vector<std::map<uint32_t, std::vector<Chunk>>> ref_bins(
        static_cast<size_t>(num_refs));
    std::vector<std::vector<int64_t>> linear(static_cast<size_t>(num_refs));

    const int64_t total = static_cast<int64_t>(impl_->recs_.size());
    for (int64_t i = 0; i < total; ++i) {
        const auto& r = impl_->recs_[static_cast<size_t>(i)];
        const int64_t rec_end = (i + 1 < total)
                                    ? impl_->recs_[static_cast<size_t>(i + 1)].off
                                    : eof_offset;
        uint32_t bin = reg2bin(r.pos, r.end);
        auto& chunks = ref_bins[static_cast<size_t>(r.ref)][bin];
        if (!chunks.empty() && chunks.back().end == r.off) {
            chunks.back().end = rec_end;
        } else {
            chunks.push_back({r.off, rec_end});
        }

        // Linear index: a record contributes its offset to every 16kb window
        // it overlaps. Keep the minimum offset per window.
        auto& lin = linear[static_cast<size_t>(r.ref)];
        int32_t first_w = r.pos >> 14;
        int32_t last_w = (r.end - 1) >> 14;
        for (int32_t w = first_w; w <= last_w; ++w) {
            size_t wi = static_cast<size_t>(w);
            if (lin.size() <= wi) lin.resize(wi + 1, -1);
            if (lin[wi] < 0 || r.off < lin[wi]) lin[wi] = r.off;
        }
    }

    for (int32_t ri = 0; ri < num_refs; ++ri) {
        const auto& bins = ref_bins[static_cast<size_t>(ri)];

        // n_bin and bin entries
        write_le32(buf, static_cast<int32_t>(bins.size()));
        for (const auto& [bin, chunks] : bins) {
            write_le32(buf, static_cast<int32_t>(bin));
            write_le32(buf, static_cast<int32_t>(chunks.size()));
            for (const auto& c : chunks) {
                write_le64(buf, static_cast<uint64_t>(c.beg));
                write_le64(buf, static_cast<uint64_t>(c.end));
            }
        }

        // Linear index. htslib expects at least ceil(ref_len/16384) entries
        // plus one sentinel; unknown windows inherit the next known offset.
        auto& lin = linear[static_cast<size_t>(ri)];
        size_t want = 0;
        if (impl_->ref_lengths_[static_cast<size_t>(ri)] > 0) {
            want = static_cast<size_t>(
                (impl_->ref_lengths_[static_cast<size_t>(ri)] + 0x3FFF) >> 14);
        }
        if (lin.size() < want) lin.resize(want, -1);
        int64_t next = eof_offset;
        for (size_t i = lin.size(); i-- > 0;) {
            if (lin[i] < 0) lin[i] = next;
            else next = lin[i];
        }
        // Sentinel entry: offset just past the last record in this reference.
        lin.push_back(eof_offset);

        write_le32(buf, static_cast<int32_t>(lin.size()));
        for (int64_t off : lin) {
            write_le64(buf, static_cast<uint64_t>(off));
        }
    }

    // Global count of unplaced unmapped reads.
    write_le64(buf, static_cast<uint64_t>(impl_->n_no_coor_));

    impl_->out_.write(reinterpret_cast<const char*>(buf.data()),
                      static_cast<std::streamsize>(buf.size()));
}

} // namespace bwa::io
