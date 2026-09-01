#pragma once

#include <bwa/core/string.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/arena.hpp>
#include <span>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <expected>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <memory>
#include <zlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

namespace bwa::io {

// Error types for I/O operations
enum class IoError {
    None = 0,
    FileNotFound,
    PermissionDenied,
    InvalidFormat,
    TruncatedRecord,
    BufferTooSmall,
    ZlibError,
    MmapFailed,
    Eof
};

inline const char* to_string(IoError e) noexcept {
    switch (e) {
        case IoError::None: return "No error";
        case IoError::FileNotFound: return "File not found";
        case IoError::PermissionDenied: return "Permission denied";
        case IoError::InvalidFormat: return "Invalid format";
        case IoError::TruncatedRecord: return "Truncated record";
        case IoError::BufferTooSmall: return "Buffer too small";
        case IoError::ZlibError: return "Zlib error";
        case IoError::MmapFailed: return "Memory mapping failed";
        case IoError::Eof: return "End of file";
        default: return "Unknown error";
    }
}

// Compressed file handle (bgzip/gzip)
class GzFile {
    gzFile gz_ = nullptr;
    bool own_ = true;

public:
    GzFile() = default;
    explicit GzFile(const char* path, const char* mode = "rb") { open(path, mode); }
    ~GzFile() { close(); }

    GzFile(const GzFile&) = delete;
    GzFile& operator=(const GzFile&) = delete;

    GzFile(GzFile&& other) noexcept : gz_(other.gz_), own_(other.own_) {
        other.gz_ = nullptr;
        other.own_ = false;
    }

    GzFile& operator=(GzFile&& other) noexcept {
        if (this != &other) {
            close();
            gz_ = other.gz_;
            own_ = other.own_;
            other.gz_ = nullptr;
            other.own_ = false;
        }
        return *this;
    }

    [[nodiscard]] bool open(const char* path, const char* mode = "rb") noexcept {
        close();
        gz_ = gzopen(path, mode);
        return gz_ != nullptr;
    }

    void close() noexcept {
        if (gz_ && own_) {
            gzclose(gz_);
            gz_ = nullptr;
        }
    }

    [[nodiscard]] bool is_open() const noexcept { return gz_ != nullptr; }
    [[nodiscard]] gzFile get() const noexcept { return gz_; }

    // Read raw bytes
    [[nodiscard]] std::expected<size_t, IoError> read(void* buf, size_t len) noexcept {
        if (!gz_) return std::unexpected(IoError::FileNotFound);
        int n = gzread(gz_, buf, static_cast<unsigned>(len));
        if (n < 0) {
            int err;
            const char* msg = gzerror(gz_, &err);
            return std::unexpected(IoError::ZlibError);
        }
        return static_cast<size_t>(n);
    }

    // Read line (up to newline or EOF)
    [[nodiscard]] std::expected<size_t, IoError> readline(core::PmrString& buf) noexcept {
        buf.clear();
        char ch;
        while (true) {
            auto result = read(&ch, 1);
            if (!result) return std::unexpected(result.error());
            if (*result == 0) {
                return buf.empty() ? std::unexpected(IoError::Eof) : std::expected<size_t, IoError>(buf.size());
            }
            if (ch == '\n') break;
            if (ch != '\r') buf.kputc(ch);
        }
        return buf.size();
    }

    // Seek (limited support for gzip)
    [[nodiscard]] bool seek(int64_t offset, int whence = SEEK_SET) noexcept {
        return gz_ && gzseek(gz_, offset, whence) != -1;
    }

    [[nodiscard]] int64_t tell() const noexcept {
        return gz_ ? gztell(gz_) : -1;
    }

    [[nodiscard]] bool eof() const noexcept {
        return gz_ && gzeof(gz_);
    }
};

// Memory-mapped file for fast random access
class MmapFile {
    void* data_ = nullptr;
    size_t size_ = 0;
    int fd_ = -1;

public:
    MmapFile() = default;
    explicit MmapFile(const char* path) { open(path); }
    ~MmapFile() { close(); }

    MmapFile(const MmapFile&) = delete;
    MmapFile& operator=(const MmapFile&) = delete;

    MmapFile(MmapFile&& other) noexcept
        : data_(other.data_), size_(other.size_), fd_(other.fd_) {
        other.data_ = nullptr;
        other.size_ = 0;
        other.fd_ = -1;
    }

    MmapFile& operator=(MmapFile&& other) noexcept {
        if (this != &other) {
            close();
            data_ = other.data_;
            size_ = other.size_;
            fd_ = other.fd_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.fd_ = -1;
        }
        return *this;
    }

    [[nodiscard]] bool open(const char* path) noexcept {
        close();
        fd_ = ::open(path, O_RDONLY);
        if (fd_ == -1) return false;

        off_t sz = lseek(fd_, 0, SEEK_END);
        if (sz == -1) { close(); return false; }
        size_ = static_cast<size_t>(sz);
        lseek(fd_, 0, SEEK_SET);

        data_ = mmap(nullptr, size_, PROT_READ, MAP_PRIVATE, fd_, 0);
        if (data_ == MAP_FAILED) { close(); return false; }
        return true;
    }

    void close() noexcept {
        if (data_ && data_ != MAP_FAILED) {
            munmap(data_, size_);
            data_ = nullptr;
        }
        if (fd_ != -1) {
            ::close(fd_);
            fd_ = -1;
        }
        size_ = 0;
    }

    [[nodiscard]] bool is_open() const noexcept { return data_ && data_ != MAP_FAILED; }
    [[nodiscard]] const void* data() const noexcept { return data_; }
    [[nodiscard]] void* data() noexcept { return data_; }
    [[nodiscard]] size_t size() const noexcept { return size_; }
    [[nodiscard]] std::span<const std::byte> span() const noexcept {
        return {static_cast<const std::byte*>(data_), size_};
    }
    [[nodiscard]] std::span<std::byte> span() noexcept {
        return {static_cast<std::byte*>(data_), size_};
    }
};

// FASTQ/FASTA record
struct SeqRecord {
    core::PmrString name;      // Sequence name (without @/>)
    core::PmrString comment;   // Comment after space in header (optional)
    core::PmrString seq;       // Sequence
    core::PmrString qual;      // Quality (empty for FASTA)

    [[nodiscard]] bool is_fasta() const noexcept { return qual.empty(); }
    [[nodiscard]] bool is_fastq() const noexcept { return !qual.empty(); }
    [[nodiscard]] size_t length() const noexcept { return seq.size(); }

    void clear() noexcept {
        name.clear();
        comment.clear();
        seq.clear();
        qual.clear();
    }
};

// Streaming FASTQ/FASTA parser (kseq.h replacement)
class SeqReader {
public:
    SeqReader() = default;
    explicit SeqReader(const char* path) { open(path); }

    [[nodiscard]] bool open(const char* path) noexcept {
        last_error_ = IoError::None;
        if (!gz_.open(path, "rb")) {
            last_error_ = IoError::FileNotFound;
            return false;
        }
        // Peek first char to detect format
        char first;
        if (gzread(gz_.get(), &first, 1) != 1) {
            last_error_ = IoError::InvalidFormat;
            return false;
        }
        gzseek(gz_.get(), 0, SEEK_SET);
        last_line_read_ = false;
        return true;
    }

    void close() noexcept { gz_.close(); }
    [[nodiscard]] bool is_open() const noexcept { return gz_.is_open(); }
    [[nodiscard]] IoError last_error() const noexcept { return last_error_; }

    // Read next record
    // Returns: true = record read, false = EOF or error (check last_error())
    [[nodiscard]] std::expected<bool, IoError> read(SeqRecord& record) noexcept {
        record.clear();
        last_error_ = IoError::None;

        // Read header line (@ or >)
        while (true) {
            auto result = gz_.readline(line_buf_);
            if (!result) {
                last_error_ = result.error();
                return std::unexpected(last_error_);
            }
            if (*result == 0) {
                last_error_ = IoError::Eof;
                return false; // EOF
            }
            if (!line_buf_.empty() && (line_buf_[0] == '>' || line_buf_[0] == '@')) {
                break;
            }
            // Skip empty lines
        }

        // Parse header: @name comment or >name comment
        std::string_view header = line_buf_.view().substr(1); // Skip @ or >
        size_t space_pos = header.find(' ');
        if (space_pos == std::string_view::npos) {
            record.name.assign(header);
        } else {
            record.name.assign(header.substr(0, space_pos));
            record.comment.assign(header.substr(space_pos + 1));
        }

        // Read sequence lines
        seq_buf_.clear();
        while (true) {
            auto result = gz_.readline(line_buf_);
            if (!result) {
                last_error_ = result.error();
                return std::unexpected(last_error_);
            }
            if (*result == 0) {
                last_error_ = IoError::TruncatedRecord;
                return std::unexpected(last_error_); // Unexpected EOF
            }
            if (!line_buf_.empty() && (line_buf_[0] == '+' || line_buf_[0] == '>')) {
                // Next record header or quality header
                last_line_read_ = true;
                break;
            }
            seq_buf_.append(line_buf_.view());
        }
        record.seq = std::move(seq_buf_);

        // Check if FASTA (next record started with >)
        if (last_line_read_ && !line_buf_.empty() && line_buf_[0] == '>') {
            // This is actually the next FASTA record - put it back by seeking
            // For simplicity, we'll store it for next read
            record.qual.clear(); // FASTA
            return true;
        }

        // FASTQ: read quality lines
        // line_buf_ now contains '+' line (may have comment after +)
        qual_buf_.clear();
        size_t qual_needed = record.seq.size();
        while (qual_buf_.size() < qual_needed) {
            auto result = gz_.readline(line_buf_);
            if (!result) {
                last_error_ = result.error();
                return std::unexpected(last_error_);
            }
            if (*result == 0) {
                last_error_ = IoError::TruncatedRecord;
                return std::unexpected(last_error_);
            }
            qual_buf_.append(line_buf_.view());
        }
        // Trim to sequence length
        if (qual_buf_.size() > qual_needed) {
            qual_buf_.resize(qual_needed);
        }
        record.qual = std::move(qual_buf_);

        last_line_read_ = false;
        return true;
    }

    // Iterator support
    class iterator {
        SeqReader* reader_;
        SeqRecord record_;
        bool valid_ = false;

        void read_next() {
            auto result = reader_->read(record_);
            valid_ = result && *result;
        }

    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = SeqRecord;
        using difference_type = std::ptrdiff_t;
        using pointer = SeqRecord*;
        using reference = SeqRecord&;

        iterator() = default;
        explicit iterator(SeqReader* r) : reader_(r) { read_next(); }

        reference operator*() noexcept { return record_; }
        pointer operator->() noexcept { return &record_; }

        iterator& operator++() { read_next(); return *this; }
        void operator++(int) { read_next(); }

        friend bool operator==(const iterator& a, const iterator& b) noexcept {
            return (!a.valid_ && !b.valid_) || (a.reader_ == b.reader_ && a.valid_ == b.valid_);
        }
        friend bool operator!=(const iterator& a, const iterator& b) noexcept {
            return !(a == b);
        }
    };

    [[nodiscard]] iterator begin() noexcept { return iterator(this); }
    [[nodiscard]] iterator end() noexcept { return iterator(); }

private:
    GzFile gz_;
    core::PmrString line_buf_;
    core::PmrString seq_buf_;
    core::PmrString qual_buf_;
    bool last_line_read_ = false;
    IoError last_error_ = IoError::None;
};

// Batch reader for parallel processing
class BatchSeqReader {
    SeqReader reader_;
    core::Vector<SeqRecord> buffer_;
    size_t batch_size_;
    size_t current_ = 0;
    bool exhausted_ = false;

public:
    BatchSeqReader(const char* path, size_t batch_size = 1000)
        : batch_size_(batch_size) {
        reader_.open(path);
        buffer_.reserve(batch_size);
    }

    [[nodiscard]] bool is_open() const noexcept { return reader_.is_open(); }
    [[nodiscard]] IoError last_error() const noexcept { return reader_.last_error(); }

    // Fill buffer, return number of records read
    size_t fill() {
        buffer_.clear();
        current_ = 0;
        SeqRecord rec;
        while (buffer_.size() < batch_size_) {
            auto result = reader_.read(rec);
            if (!result || !*result) {
                exhausted_ = true;
                break;
            }
            buffer_.push_back(std::move(rec));
        }
        return buffer_.size();
    }

    [[nodiscard]] std::span<SeqRecord> current_batch() noexcept {
        return buffer_.span();
    }

    [[nodiscard]] bool next_batch() {
        if (exhausted_) return false;
        return fill() > 0;
    }

    [[nodiscard]] bool exhausted() const noexcept { return exhausted_; }
};

// Writer for FASTQ/FASTA output
class SeqWriter {
    GzFile gz_;
    core::PmrString out_buf_;

public:
    SeqWriter() = default;
    explicit SeqWriter(const char* path, const char* mode = "wb") { open(path, mode); }

    [[nodiscard]] bool open(const char* path, const char* mode = "wb") noexcept {
        return gz_.open(path, mode);
    }

    void close() noexcept { gz_.close(); }
    [[nodiscard]] bool is_open() const noexcept { return gz_.is_open(); }

    // Write FASTA record
    [[nodiscard]] bool write_fasta(std::string_view name, std::string_view comment,
                                   std::string_view seq, int line_width = 80) noexcept {
        out_buf_.clear();
        out_buf_.kputc('>');
        out_buf_.kputsn(name.data(), name.size());
        if (!comment.empty()) {
            out_buf_.kputc(' ');
            out_buf_.kputsn(comment.data(), comment.size());
        }
        out_buf_.kputc('\n');

        for (size_t i = 0; i < seq.size(); i += line_width) {
            size_t len = std::min<size_t>(line_width, seq.size() - i);
            out_buf_.kputsn(seq.data() + i, len);
            out_buf_.kputc('\n');
        }

        auto result = gz_.read(out_buf_.data(), out_buf_.size());
        return result && *result == out_buf_.size();
    }

    // Write FASTQ record
    [[nodiscard]] bool write_fastq(std::string_view name, std::string_view comment,
                                   std::string_view seq, std::string_view qual,
                                   int line_width = 80) noexcept {
        out_buf_.clear();
        out_buf_.kputc('@');
        out_buf_.kputsn(name.data(), name.size());
        if (!comment.empty()) {
            out_buf_.kputc(' ');
            out_buf_.kputsn(comment.data(), comment.size());
        }
        out_buf_.kputc('\n');

        for (size_t i = 0; i < seq.size(); i += line_width) {
            size_t len = std::min<size_t>(line_width, seq.size() - i);
            out_buf_.kputsn(seq.data() + i, len);
            out_buf_.kputc('\n');
        }

        out_buf_.kputc('+');
        if (!comment.empty()) {
            out_buf_.kputc(' ');
            out_buf_.kputsn(comment.data(), comment.size());
        }
        out_buf_.kputc('\n');

        for (size_t i = 0; i < qual.size(); i += line_width) {
            size_t len = std::min<size_t>(line_width, qual.size() - i);
            out_buf_.kputsn(qual.data() + i, len);
            out_buf_.kputc('\n');
        }

        auto result = gz_.read(out_buf_.data(), out_buf_.size());
        return result && *result == out_buf_.size();
    }

    [[nodiscard]] bool write(const SeqRecord& rec, int line_width = 80) noexcept {
        if (rec.is_fasta()) {
            return write_fasta(rec.name.view(), rec.comment.view(), rec.seq.view(), line_width);
        } else {
            return write_fastq(rec.name.view(), rec.comment.view(), rec.seq.view(), rec.qual.view(), line_width);
        }
    }

    void flush() noexcept {
        if (gz_.is_open()) {
            gzflush(gz_.get(), Z_SYNC_FLUSH);
        }
    }
};

} // namespace bwa::io