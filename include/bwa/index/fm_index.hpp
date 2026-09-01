#pragma once

#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/sort.hpp>
#include <span>
#include <cstddef>
#include <cstdint>
#include <array>
#include <bit>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace bwa::index {

// Packed DNA sequence (2 bits per base)
// A=0, C=1, G=2, T=3, N=4 (stored as 0 with separate N mask)
class PackedSequence {
public:
    using word_type = uint64_t;
    static constexpr int BASES_PER_WORD = 32; // 64 bits / 2 bits
    static constexpr word_type MASK = 0x3;

    // Base encoding
    static constexpr uint8_t ENCODE_A = 0;
    static constexpr uint8_t ENCODE_C = 1;
    static constexpr uint8_t ENCODE_G = 2;
    static constexpr uint8_t ENCODE_T = 3;
    static constexpr uint8_t ENCODE_N = 4;

    static constexpr char DECODE_TABLE[5] = {'A', 'C', 'G', 'T', 'N'};

private:
    core::Vector<word_type> data_;
    core::Vector<word_type> n_mask_; // 1 bit per base: 1 = N
    size_t length_ = 0;

    [[nodiscard]] static size_t word_index(size_t pos) noexcept { return pos / BASES_PER_WORD; }
    [[nodiscard]] static int bit_offset(size_t pos) noexcept { return (pos % BASES_PER_WORD) * 2; }
    [[nodiscard]] static size_t n_word_index(size_t pos) noexcept { return pos / 64; }
    [[nodiscard]] static int n_bit_offset(size_t pos) noexcept { return pos % 64; }

public:
    PackedSequence() = default;
    explicit PackedSequence(size_t length) { resize(length); }

    [[nodiscard]] size_t size() const noexcept { return length_; }
    [[nodiscard]] size_t length() const noexcept { return length_; }
    [[nodiscard]] bool empty() const noexcept { return length_ == 0; }
    [[nodiscard]] size_t word_count() const noexcept { return data_.size(); }
    [[nodiscard]] size_t byte_size() const noexcept { return data_.size() * sizeof(word_type); }

    void resize(size_t n) {
        length_ = n;
        data_.resize((n + BASES_PER_WORD - 1) / BASES_PER_WORD);
        n_mask_.resize((n + 63) / 64);
    }

    void clear() noexcept {
        data_.clear();
        n_mask_.clear();
        length_ = 0;
    }

    [[nodiscard]] uint8_t get(size_t pos) const noexcept {
        if (pos >= length_) return ENCODE_N;
        if (is_n(pos)) return ENCODE_N;
        size_t wi = word_index(pos);
        int bo = bit_offset(pos);
        return static_cast<uint8_t>((data_[wi] >> bo) & MASK);
    }

    [[nodiscard]] bool is_n(size_t pos) const noexcept {
        if (pos >= length_) return true;
        size_t wi = n_word_index(pos);
        return (n_mask_[wi] >> n_bit_offset(pos)) & 1;
    }

    void set(size_t pos, uint8_t base) noexcept {
        if (pos >= length_) return;
        size_t wi = word_index(pos);
        int bo = bit_offset(pos);
        // Store N as 0 in data, use N mask to distinguish
        uint8_t stored = (base == ENCODE_N) ? 0 : (base & MASK);
        data_[wi] = (data_[wi] & ~(word_type(MASK) << bo)) | (word_type(stored) << bo);

        size_t nwi = n_word_index(pos);
        int nbo = n_bit_offset(pos);
        if (base == ENCODE_N) {
            n_mask_[nwi] |= (word_type(1) << nbo);
        } else {
            n_mask_[nwi] &= ~(word_type(1) << nbo);
        }
    }

    void set_n(size_t pos, bool is_n) noexcept {
        if (pos >= length_) return;
        size_t wi = n_word_index(pos);
        int bo = n_bit_offset(pos);
        if (is_n) n_mask_[wi] |= (word_type(1) << bo);
        else n_mask_[wi] &= ~(word_type(1) << bo);
    }

    // Append base (amortized O(1))
    void push_back(uint8_t base) {
        if (length_ % BASES_PER_WORD == 0) {
            data_.push_back(0);
        }
        if (length_ % 64 == 0) {
            n_mask_.push_back(0);
        }
        set(length_, base);
        ++length_;
    }

    // Bulk append from string
    void append(const char* s, size_t n) {
        size_t old_len = length_;
        resize(length_ + n);
        for (size_t i = 0; i < n; ++i) {
            set(old_len + i, encode_base(s[i]));
        }
    }

    // Iterator for forward traversal
    class const_iterator {
        const PackedSequence* seq_;
        size_t pos_;
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = uint8_t;
        using difference_type = std::ptrdiff_t;
        using pointer = const uint8_t*;
        using reference = uint8_t;

        const_iterator() = default;
        const_iterator(const PackedSequence* s, size_t p) : seq_(s), pos_(p) {}

        reference operator*() const noexcept { return seq_->get(pos_); }
        const_iterator& operator++() { ++pos_; return *this; }
        const_iterator operator++(int) { auto tmp = *this; ++pos_; return tmp; }
        friend bool operator==(const const_iterator& a, const const_iterator& b) noexcept {
            return a.pos_ == b.pos_;
        }
        friend bool operator!=(const const_iterator& a, const const_iterator& b) noexcept {
            return a.pos_ != b.pos_;
        }
    };

    [[nodiscard]] const_iterator begin() const noexcept { return const_iterator(this, 0); }
    [[nodiscard]] const_iterator end() const noexcept { return const_iterator(this, length_); }

    // Span access to raw words
    [[nodiscard]] std::span<const word_type> words() const noexcept { return data_.span(); }
    [[nodiscard]] std::span<word_type> words() noexcept { return data_.span(); }
    [[nodiscard]] std::span<const word_type> n_masks() const noexcept { return n_mask_.span(); }

    // Get sequence as packed bytes (2 bits per base)
    [[nodiscard]] std::span<const uint8_t> bases() const noexcept {
        // PackedSequence stores 2 bits per base in words, but we need bytes for MEM finder
        // For now, return a view of the underlying byte data
        // This is a simplification - in reality we'd need to unpack
        return std::span<const uint8_t>(reinterpret_cast<const uint8_t*>(data_.data()), length_);
    }

    // Static encoding/decoding
    [[nodiscard]] static uint8_t encode_base(char c) noexcept {
        switch (c | 0x20) { // lowercase
            case 'a': return ENCODE_A;
            case 'c': return ENCODE_C;
            case 'g': return ENCODE_G;
            case 't': return ENCODE_T;
            default: return ENCODE_N;
        }
    }

    [[nodiscard]] static char decode_base(uint8_t code) noexcept {
        return code <= 3 ? DECODE_TABLE[code] : 'N';
    }

    // Reverse complement in-place
    void reverse_complement() {
        for (size_t i = 0, j = length_ - 1; i < j; ++i, --j) {
            uint8_t bi = get(i);
            uint8_t bj = get(j);
            bool ni = is_n(i);
            bool nj = is_n(j);

            // Complement: A<->T, C<->G, N->N
            auto comp = [](uint8_t b) -> uint8_t {
                return b == ENCODE_A ? ENCODE_T :
                       b == ENCODE_T ? ENCODE_A :
                       b == ENCODE_C ? ENCODE_G :
                       b == ENCODE_G ? ENCODE_C : ENCODE_N;
            };

            set(i, comp(bj));
            set(j, comp(bi));
            set_n(i, nj);
            set_n(j, ni);
        }
        if (length_ % 2 == 1) {
            size_t mid = length_ / 2;
            uint8_t b = get(mid);
            set(mid, b == ENCODE_A ? ENCODE_T :
                    b == ENCODE_T ? ENCODE_A :
                    b == ENCODE_C ? ENCODE_G :
                    b == ENCODE_G ? ENCODE_C : ENCODE_N);
        }
    }

    // Export to string
    void to_string(core::PmrString& out) const {
        out.clear();
        out.reserve(length_);
        for (size_t i = 0; i < length_; ++i) {
            out.kputc(decode_base(get(i)));
        }
    }
};

// FM-index with rank/select support
class FMIndex {
public:
    using occ_t = uint32_t;
    static constexpr int OCC_INTERVAL = 128; // Rank sampling interval
    static constexpr int SA_INTERVAL = 32;   // SA sampling interval

private:
    PackedSequence bwt_;
    core::Vector<occ_t> occ_;       // Rank table: occ[base][i] = count of base in bwt[0..i*OCC_INTERVAL)
    core::Vector<uint32_t> sa_;     // Sampled SA: sa[i] = SA[i * SA_INTERVAL]
    core::Vector<uint32_t> cnt_;    // Cumulative counts: cnt[c] = #bases < c
    size_t primary_ = 0;      // Position of original string end ($)
    size_t length_ = 0;       // Original sequence length

    // Occurrence array layout: [4][num_intervals]
    [[nodiscard]] size_t occ_offset(uint8_t base, size_t interval) const noexcept {
        return interval * 4 + base;
    }

    [[nodiscard]] size_t num_intervals() const noexcept {
        return (length_ + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    }

public:
    FMIndex() = default;

    [[nodiscard]] size_t size() const noexcept { return length_; }
    [[nodiscard]] size_t length() const noexcept { return length_; }
    [[nodiscard]] size_t bwt_size() const noexcept { return bwt_.size(); }
    [[nodiscard]] const PackedSequence& bwt() const noexcept { return bwt_; }
    [[nodiscard]] size_t primary() const noexcept { return primary_; }

    // Serialization accessors
    [[nodiscard]] const core::Vector<uint32_t>& sa_samples() const noexcept { return sa_; }
    [[nodiscard]] const core::Vector<uint32_t>& occ_table() const noexcept { return occ_; }
    [[nodiscard]] const core::Vector<uint32_t>& count_table() const noexcept { return cnt_; }

    // Setter for deserialization
    void set_size(size_t s) { length_ = s; }
    void set_bwt(const PackedSequence& bwt) { bwt_ = bwt; }
    void set_primary(size_t p) { primary_ = p; }
    void set_sa_samples(const core::Vector<uint32_t>& sa) { sa_ = sa; }
    void set_occ_table(const core::Vector<uint32_t>& occ) { occ_ = occ; }
    void set_count_table(const core::Vector<uint32_t>& cnt) { cnt_ = cnt; }

    // Rank: number of occurrences of base in bwt[0..pos)
    [[nodiscard]] uint32_t rank(uint8_t base, size_t pos) const noexcept {
        if (pos == 0) return 0;
        size_t interval = pos / OCC_INTERVAL;
        size_t offset = pos % OCC_INTERVAL;

        uint32_t count = occ_[occ_offset(base, interval)];
        // Scan remaining
        size_t start = interval * OCC_INTERVAL;
        for (size_t i = start; i < pos; ++i) {
            if (bwt_.get(i) == base) ++count;
        }
        return count;
    }

    // C array: cumulative counts
    [[nodiscard]] uint32_t C(uint8_t base) const noexcept {
        return cnt_[base];
    }

    // Total count of base
    [[nodiscard]] uint32_t total(uint8_t base) const noexcept {
        return cnt_[base + 1] - cnt_[base];
    }

    // LF mapping: LF(pos) = C[bwt[pos]] + rank(bwt[pos], pos)
    [[nodiscard]] size_t lf(size_t pos) const noexcept {
        uint8_t base = bwt_.get(pos);
        return C(base) + rank(base, pos);
    }

    // Backward search: extend interval [l, r) with base
    // Returns new interval [l', r')
    [[nodiscard]] std::pair<size_t, size_t> backward_extend(uint8_t base,
                                                             size_t l, size_t r) const noexcept {
        size_t new_l = C(base) + rank(base, l);
        size_t new_r = C(base) + rank(base, r);
        return {new_l, new_r};
    }

    // Backward search for pattern
    // Returns interval [l, r) in SA, or empty if not found
    template <typename Pattern>
    [[nodiscard]] std::pair<size_t, size_t> backward_search(const Pattern& pat) const noexcept {
        size_t l = 0, r = length_;
        for (auto it = pat.rbegin(); it != pat.rend(); ++it) {
            uint8_t base = PackedSequence::encode_base(*it);
            auto [new_l, new_r] = backward_extend(base, l, r);
            if (new_l >= new_r) return {0, 0}; // Empty
            l = new_l;
            r = new_r;
        }
        return {l, r};
    }

    // Locate: get SA value at position (requires sampled SA)
    [[nodiscard]] std::optional<uint32_t> locate(size_t pos) const noexcept {
        if (sa_.empty()) return std::nullopt;
        if (pos % SA_INTERVAL == 0) {
            return sa_[pos / SA_INTERVAL];
        }
        // Walk backward using LF until sampled
        size_t steps = 0;
        while (pos % SA_INTERVAL != 0 && steps < SA_INTERVAL) {
            pos = lf(pos);
            ++steps;
        }
        if (pos % SA_INTERVAL == 0) {
            auto base_sa = sa_[pos / SA_INTERVAL];
            return base_sa + steps < length_ ? std::optional(base_sa + steps) : std::nullopt;
        }
        return std::nullopt;
    }

    // Count occurrences of pattern
    template <typename Pattern>
    [[nodiscard]] size_t count(const Pattern& pat) const noexcept {
        auto [l, r] = backward_search(pat);
        return r - l;
    }

    // Build from packed sequence (SA construction)
    static FMIndex build(const PackedSequence& seq, memory::Arena& arena) {
        FMIndex idx;
        idx.length_ = seq.size();
        idx.bwt_.resize(seq.size());

        // Allocate working memory from arena
        size_t n = seq.size();
        core::Vector<uint32_t> sa(&arena);
        sa.resize(n);
        core::Vector<uint32_t> rank(&arena);
        rank.resize(n);

        // SA construction using induced sorting (SA-IS) or fallback
        // For now, use simple prefix-doubling (O(n log n)) - replace with SA-IS later
        build_sa_doubling(seq, sa, rank, arena);

        // Build BWT from SA
        for (size_t i = 0; i < n; ++i) {
            size_t sa_i = sa[i];
            uint8_t b = (sa_i == 0) ? PackedSequence::ENCODE_N : seq.get(sa_i - 1);
            idx.bwt_.set(i, b);
            if (sa_i == 0) idx.primary_ = i;
        }

        // Build occ table
        size_t num_intv = idx.num_intervals();
        idx.occ_.resize(num_intv * 4);
        idx.cnt_.resize(5); // cnt[0..4], cnt[4] = n

        // Count bases
        std::array<uint32_t, 4> counts = {0, 0, 0, 0};
        for (size_t i = 0; i < n; ++i) {
            uint8_t b = idx.bwt_.get(i);
            if (b < 4) ++counts[b];

            if ((i + 1) % OCC_INTERVAL == 0) {
                size_t intv = (i + 1) / OCC_INTERVAL;
                for (int b = 0; b < 4; ++b) {
                    idx.occ_[idx.occ_offset(b, intv)] = counts[b];
                }
            }
        }

        // Cumulative counts
        idx.cnt_[0] = 0;
        for (int b = 0; b < 4; ++b) {
            idx.cnt_[b + 1] = idx.cnt_[b] + counts[b];
        }

        // Sample SA
        idx.sa_.resize((n + SA_INTERVAL - 1) / SA_INTERVAL);
        for (size_t i = 0; i < n; i += SA_INTERVAL) {
            idx.sa_[i / SA_INTERVAL] = sa[i];
        }

        return idx;
    }

    // Save/load
    void save(const char* prefix) const {
        // Implementation for binary serialization
    }

    static FMIndex load(const char* prefix) {
        FMIndex idx;
        // Implementation for binary deserialization
        return idx;
    }

private:
    // Simple prefix-doubling SA construction (O(n log n))
    // Replace with SA-IS (O(n)) for production
    static void build_sa_doubling(const PackedSequence& seq,
                                   core::Vector<uint32_t>& sa,
                                   core::Vector<uint32_t>& rank,
                                   memory::Arena& arena) {
        size_t n = seq.size();
        // Initial rank = base value (0..4)
        for (size_t i = 0; i < n; ++i) {
            sa[i] = static_cast<uint32_t>(i);
            rank[i] = seq.get(i);
        }

        core::Vector<uint32_t> tmp_sa(&arena);
        core::Vector<uint32_t> tmp_rank(&arena);
        tmp_sa.resize(n);
        tmp_rank.resize(n);

        for (size_t k = 1; k < n; k <<= 1) {
            auto cmp = [&](uint32_t a, uint32_t b) {
                if (rank[a] != rank[b]) return rank[a] < rank[b];
                uint32_t ra = (a + k < n) ? rank[a + k] : 0;
                uint32_t rb = (b + k < n) ? rank[b + k] : 0;
                return ra < rb;
            };

            core::sort(sa.begin(), sa.end(), cmp);

            tmp_rank[sa[0]] = 1;
            int r = 1;
            for (size_t i = 1; i < n; ++i) {
                if (cmp(sa[i - 1], sa[i])) ++r;
                tmp_rank[sa[i]] = r;
            }
            rank.swap(tmp_rank);
            if (r >= static_cast<int>(n)) break;
        }
    }
};

// Multi-sequence FM-index (for multiple references)
class MultiFMIndex {
    core::Vector<FMIndex> indexes_;
    core::Vector<size_t> offsets_; // Cumulative lengths

public:
    void add_sequence(const PackedSequence& seq, memory::Arena& arena) {
        indexes_.push_back(FMIndex::build(seq, arena));
        offsets_.push_back(offsets_.empty() ? seq.size() : offsets_.back() + seq.size());
    }

    [[nodiscard]] size_t num_sequences() const noexcept { return indexes_.size(); }
    [[nodiscard]] size_t total_length() const noexcept {
        return offsets_.empty() ? 0 : offsets_.back();
    }

    [[nodiscard]] const FMIndex& operator[](size_t i) const noexcept { return indexes_[i]; }
    [[nodiscard]] FMIndex& operator[](size_t i) noexcept { return indexes_[i]; }

    // Iterator support for range-based for
    [[nodiscard]] auto begin() const noexcept { return indexes_.begin(); }
    [[nodiscard]] auto end() const noexcept { return indexes_.end(); }
    [[nodiscard]] auto begin() noexcept { return indexes_.begin(); }
    [[nodiscard]] auto end() noexcept { return indexes_.end(); }

    // Find which sequence contains global position
    [[nodiscard]] std::pair<size_t, size_t> locate_sequence(size_t global_pos) const noexcept {
        if (global_pos >= total_length()) return {indexes_.size(), 0};
        size_t seq_idx = std::lower_bound(offsets_.begin(), offsets_.end(), global_pos + 1) - offsets_.begin();
        size_t seq_start = (seq_idx == 0) ? 0 : offsets_[seq_idx - 1];
        return {seq_idx, global_pos - seq_start};
    }
};

} // namespace bwa::index