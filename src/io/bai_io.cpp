#include <bwa/io/bai_io.hpp>
#include <bwa/io/bam_io.hpp>
#include <vector>
#include <string>
#include <cstdint>
#include <memory>

namespace bwa::io {

struct BaiWriter::Impl {
    BgzfWriter bgzf_;
    struct RefStats {
        int32_t n_mapped = 0;
        int32_t n_unmapped = 0;
        std::vector<std::vector<int64_t>> bin_offsets;
        std::vector<std::pair<int64_t, int64_t>> linear_offsets;

        RefStats() : bin_offsets(4681) {}
    };
    std::vector<RefStats> ref_stats_;
    int32_t num_refs_ = 0;
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
    return impl_->bgzf_.open(path);
}

void BaiWriter::close() {
    impl_->bgzf_.close();
}

void BaiWriter::init(int32_t num_refs) {
    impl_->num_refs_ = num_refs;
    impl_->ref_stats_.assign(num_refs, Impl::RefStats{});
}

void BaiWriter::record_alignment(int32_t ref_idx, int32_t pos, int64_t file_offset) {
    if (ref_idx < 0 || ref_idx >= impl_->num_refs_) return;

    auto& stats = impl_->ref_stats_[ref_idx];
    if (pos < 0) {
        stats.n_unmapped++;
        return;
    }
    stats.n_mapped++;

    uint16_t bin = 4680;
    if (pos < 65536) bin = pos >> 14;
    else if (pos < 262144) bin = 9 + (pos >> 17);
    else if (pos < 1048576) bin = 15 + (pos >> 20);
    else if (pos < 4194304) bin = 21 + (pos >> 23);
    else if (pos < 16777216) bin = 27 + (pos >> 26);

    impl_->ref_stats_[ref_idx].bin_offsets[bin].push_back(file_offset);

    int32_t linear_idx = pos >> 14;
    if (linear_idx >= 0) {
        auto& stats = impl_->ref_stats_[ref_idx];
        if (static_cast<int32_t>(stats.linear_offsets.size()) <= linear_idx) {
            stats.linear_offsets.resize(linear_idx + 1, {0, 0});
        }
        if (stats.linear_offsets[linear_idx].first == 0 ||
            file_offset < stats.linear_offsets[linear_idx].first) {
            stats.linear_offsets[linear_idx] = {file_offset, 0};
        }
    }
}

void BaiWriter::write_index() {
    std::vector<uint8_t> buf;

    buf.push_back('B');
    buf.push_back('A');
    buf.push_back('I');
    buf.push_back(1);

    write_le32(buf, impl_->num_refs_);

    for (int32_t i = 0; i < impl_->num_refs_; ++i) {
        auto& stats = impl_->ref_stats_[i];

        int32_t n_bin = 0;
        for (int b = 0; b < 16; ++b) {
            if (!stats.bin_offsets[b].empty()) n_bin++;
        }
        if (!stats.bin_offsets[4680].empty()) n_bin++;
        write_le32(buf, n_bin);

        for (int b = 0; b < 16; ++b) {
            if (!stats.bin_offsets[b].empty()) {
                write_le32(buf, b);
                int32_t n_chunk = static_cast<int32_t>(stats.bin_offsets[b].size());
                write_le32(buf, n_chunk);
                for (int64_t offset : stats.bin_offsets[b]) {
                    write_le64(buf, static_cast<uint64_t>(offset));
                }
            }
        }
        if (!stats.bin_offsets[4680].empty()) {
            write_le32(buf, 4680);
            int32_t n_chunk = static_cast<int32_t>(stats.bin_offsets[4680].size());
            write_le32(buf, n_chunk);
            for (int64_t offset : stats.bin_offsets[4680]) {
                write_le64(buf, static_cast<uint64_t>(offset));
            }
        }

        int32_t n_intv = static_cast<int32_t>(stats.linear_offsets.size());
        write_le32(buf, n_intv);
        for (int32_t j = 0; j < n_intv; ++j) {
            write_le64(buf, static_cast<uint64_t>(stats.linear_offsets[j].first));
        }

        write_le64(buf, static_cast<uint64_t>(stats.n_mapped));
        write_le64(buf, static_cast<uint64_t>(stats.n_unmapped));
    }

    impl_->bgzf_.write(buf);
}

} // namespace bwa::io