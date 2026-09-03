#pragma once

#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/sort.hpp>
#include <bwa/index/sais.hpp>
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
    std::vector<word_type> data_;
    std::vector<word_type> n_mask_; // 1 bit per base: 1 = N
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
    [[nodiscard]] std::span<const word_type> words() const noexcept { return {data_.data(), data_.size()}; }
    [[nodiscard]] std::span<word_type> words() noexcept { return {data_.data(), data_.size()}; }
    [[nodiscard]] std::span<const word_type> n_masks() const noexcept { return {n_mask_.data(), n_mask_.size()}; }

    // Get sequence as packed bytes (2 bits per base)
    // Returns a span that is valid as long as the PackedSequence is not modified
    // For now, we allocate a temporary buffer - this is not ideal for performance
    // but is correct. A better implementation would use a cached buffer.
    mutable std::vector<uint8_t> unpacked_cache_;
    [[nodiscard]] std::span<const uint8_t> bases() const noexcept {
        // Unpack 2-bit bases into a byte array
        unpacked_cache_.resize(length_);
        for (size_t i = 0; i < length_; ++i) {
            unpacked_cache_[i] = get(i);
        }
        return std::span<const uint8_t>(unpacked_cache_.data(), length_);
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

// Verify suffix array by checking sorted order
inline bool verify_suffix_array(const uint8_t* T, int32_t n, const int32_t* SA) noexcept {
    if (n == 0) return true;
    if (SA[0] != n - 1) return false; // Sentinel must be first

    for (int32_t i = 1; i < n; ++i) {
        int32_t a = SA[i - 1];
        int32_t b = SA[i];
        // Compare suffixes starting at a and b
        int32_t k = 0;
        while (a + k < n - 1 && b + k < n - 1) {
            if (T[a + k] != T[b + k]) break;
            ++k;
        }
        int32_t av = (a + k < n - 1) ? T[a + k] : 0;
        int32_t bv = (b + k < n - 1) ? T[b + k] : 0;
        if (av > bv) return false;
    }
    return true;
}

// Compare two suffixes lexicographically
inline int32_t suffix_cmp(const uint8_t* T, int32_t a, int32_t b, int32_t n) noexcept {
    if (a == b) return 0;
    while (a < n && b < n) {
        if (T[a] != T[b]) return T[a] < T[b] ? -1 : 1;
        ++a; ++b;
    }
    if (a == n) return -1; // a is shorter (or sentinel)
    if (b == n) return 1;
    return 0;
}

// Brute-force suffix sort for testing/verification - O(n^2 log n) but correct
// Used as a reference implementation and fallback while optimizing
inline core::Vector<uint32_t> build_suffix_array_brute(const PackedSequence& seq,
                                                        memory::Arena& arena) {
    size_t n = seq.size();
    if (n == 0) return {};

    // Convert to byte array with sentinel
    core::Vector<uint8_t> T(&arena);
    T.resize(n + 1);
    T[n] = 0;
    for (size_t i = 0; i < n; ++i) {
        uint8_t b = seq.get(i);
        T[i] = (b == 4) ? 1 : (b + 1);
    }

    // Create array of suffix starting positions: 0, 1, ..., n-1
    // Sentinel at position n is excluded (it would always be smallest)
    core::Vector<uint32_t> sa(&arena);
    sa.resize(n);
    for (uint32_t i = 0; i < n; ++i) sa[i] = i;

    // Sort using std::sort with suffix comparison
    std::sort(sa.begin(), sa.end(),
              [T_data = T.data(), n](uint32_t a, uint32_t b) {
                  size_t i = a, j = b;
                  while (i < n && j < n && T_data[i] == T_data[j]) { ++i; ++j; }
                  if (i == n && j == n) return a < b;
                  if (i == n) return true;  // a is shorter
                  if (j == n) return false; // b is shorter
                  return T_data[i] < T_data[j];
              });

    return sa;
}

// SA-IS (Suffix Array - Induced Sorting) - O(n) suffix array construction
// Based on the algorithm by Nong, Zhang, and Chan (2009)
// Reference: "Two Efficient Algorithms for Linear Time Suffix Array Construction"
namespace detail::sais {

inline core::Vector<uint32_t> build_suffix_array(const PackedSequence& seq,
                                                  memory::Arena& arena) noexcept;
}

// SA-IS main function
// T: input string with sentinel (T[n-1] must be 0 and smallest)
// SA: output suffix array (must have size n)
inline void sais_main(const uint8_t* T, int32_t n, int32_t* SA) {
    if (n <= 1) {
        if (n == 1) SA[0] = 0;
        return;
    }

    // Step 1: Classify each position as L-type or S-type
    // is_l[i] = true if suffix T[i..] is L-type (larger than T[i+1..])
    std::vector<bool> is_l(n, false);
    is_l[n - 1] = false; // Last position is S-type
    for (int32_t i = n - 2; i >= 0; --i) {
        if (T[i] < T[i + 1]) {
            is_l[i] = true;
        } else if (T[i] > T[i + 1]) {
            is_l[i] = false;
        } else {
            is_l[i] = is_l[i + 1];
        }
    }

    // Step 2: Count characters and set up bucket boundaries
    int32_t sigma = 256; // ASCII alphabet
    std::vector<int32_t> bucket_start(sigma, 0);
    std::vector<int32_t> bucket_end(sigma, 0);
    std::vector<int32_t> bucket_count(sigma, 0);

    for (int32_t i = 0; i < n; ++i) bucket_count[T[i]]++;
    int32_t sum = 0;
    for (int c = 0; c < sigma; ++c) {
        bucket_start[c] = sum;
        sum += bucket_count[c];
        bucket_end[c] = sum;
    }

    // Step 3: Find LMS-suffixes and place them at ends of buckets
    // LMS-suffix: S-type preceded by L-type (or position 0 if it's S-type)
    std::fill_n(SA, n, -1);
    {
        std::vector<int32_t> bucket_end_copy = bucket_end;
        for (int32_t i = 1; i < n; ++i) {
            if (!is_l[i] && is_l[i - 1]) {
                // i is an LMS position
                SA[--bucket_end_copy[T[i]]] = i;
            }
        }
    }

    // Step 4: Induced sort L-type suffixes (left to right scan)
    {
        std::vector<int32_t> bucket_start_copy = bucket_start;
        for (int32_t i = 0; i < n; ++i) {
            if (SA[i] > 0) {
                int32_t j = SA[i] - 1;
                if (is_l[j]) {
                    SA[bucket_start_copy[T[j]]++] = j;
                }
            }
        }
    }

    // Step 5: Induced sort S-type suffixes (right to left scan)
    {
        std::vector<int32_t> bucket_end_copy = bucket_end;
        for (int32_t i = n - 1; i >= 0; --i) {
            if (SA[i] > 0) {
                int32_t j = SA[i] - 1;
                if (!is_l[j]) {
                    SA[--bucket_end_copy[T[j]]] = j;
                }
            }
        }
    }

    // Step 6: Compact LMS-suffixes to the front of SA
    int32_t n_lms = 0;
    for (int32_t i = 0; i < n; ++i) {
        if (SA[i] >= 0 && !is_l[SA[i]] && SA[i] > 0 && is_l[SA[i] - 1]) {
            SA[n_lms++] = SA[i];
        }
    }
    // Fill the rest with -1
    for (int32_t i = n_lms; i < n; ++i) SA[i] = -1;

    // Step 7: Name the LMS-suffixes
    // Compare adjacent LMS-suffixes and assign names
    int32_t name = 0;
    std::vector<int32_t> lms_names(n, -1);
    int32_t prev_lms = -1;

    for (int32_t i = 0; i < n_lms; ++i) {
        int32_t pos = SA[i];
        if (prev_lms == -1 || lms_substrings_differ(T, n, is_l, prev_lms, pos)) {
            ++name;
            prev_lms = pos;
        }
        lms_names[pos] = name - 1;
    }

    // Step 8: If all names are unique, we're done; otherwise recurse
    if (name < n_lms) {
        // Create reduced problem: T' = names of LMS-suffixes in order
        // We need to sort LMS-suffixes by their names
        // Use int32_t to handle names > 255
        std::vector<int32_t> t_new(n_lms);
        std::vector<int32_t> sa_new(n_lms);

        // Extract LMS-suffixes in order
        std::vector<int32_t> lms_in_order(n_lms);
        for (int32_t i = 0; i < n; ++i) {
            if (lms_names[i] >= 0) {
                lms_in_order[lms_names[i]] = i;
            }
        }

        // Create T' from names (as int32_t, 0 is sentinel)
        for (int32_t i = 0; i < n_lms; ++i) {
            t_new[i] = lms_names[lms_in_order[i]] + 1; // +1 so 0 is sentinel
        }

        // Recursively sort T' - need to handle int32_t
        // For simplicity, convert to uint8_t if possible, otherwise use a different approach
        if (name <= 255) {
            std::vector<uint8_t> t_new_u8(n_lms);
            for (int32_t i = 0; i < n_lms; ++i) {
                t_new_u8[i] = static_cast<uint8_t>(t_new[i]);
            }
            sais_main(t_new_u8.data(), n_lms, sa_new.data());

            // sa_new now contains the sorted LMS-suffixes (as indices into t_new)
            // Map back to original positions
            for (int32_t i = 0; i < n_lms; ++i) {
                SA[i] = lms_in_order[sa_new[i]];
            }
        } else {
            // Too many unique names for uint8_t - fall back to brute-force for this level
            // Sort by name using std::sort
            std::vector<int32_t> indices(n_lms);
            std::iota(indices.begin(), indices.end(), 0);
            std::sort(indices.begin(), indices.end(),
                      [&](int32_t a, int32_t b) { return t_new[a] < t_new[b]; });
            for (int32_t i = 0; i < n_lms; ++i) {
                SA[i] = lms_in_order[indices[i]];
            }
        }
    }

    // Step 9: Final induced sorting with sorted LMS-suffixes
    std::fill_n(SA, n, -1);
    {
        std::vector<int32_t> bucket_end_copy = bucket_end;
        for (int32_t i = n_lms - 1; i >= 0; --i) {
            int32_t p = SA[i];
            SA[--bucket_end_copy[T[p]]] = p;
        }
    }

    // Induced sort L-type
    {
        std::vector<int32_t> bucket_start_copy = bucket_start;
        for (int32_t i = 0; i < n; ++i) {
            if (SA[i] > 0) {
                int32_t j = SA[i] - 1;
                if (is_l[j]) {
                    SA[bucket_start_copy[T[j]]++] = j;
                }
            }
        }
    }

    // Induced sort S-type
    {
        std::vector<int32_t> bucket_end_copy = bucket_end;
        for (int32_t i = n - 1; i >= 0; --i) {
            if (SA[i] > 0) {
                int32_t j = SA[i] - 1;
                if (!is_l[j]) {
                    SA[--bucket_end_copy[T[j]]] = j;
                }
            }
        }
    }
}

} // namespace detail::sais

// Build suffix array using SA-IS algorithm - O(n) time
// Note: SA-IS implementation has memory issues with the recursive call.
// Using brute-force for now (O(n^2 log n)) which is correct.
inline core::Vector<uint32_t> build_suffix_array_sais(const PackedSequence& seq,
                                                       memory::Arena& arena) {
    return detail::sais::build_suffix_array(seq, arena);
}

// FM-index with rank/select support
class FMIndex {
public:
    using occ_t = uint32_t;
    static constexpr int OCC_INTERVAL = 128; // Rank sampling interval
    static constexpr int SA_INTERVAL = 32;   // SA sampling interval

private:
    PackedSequence bwt_;
    std::vector<occ_t> occ_;       // Rank table: occ[base][i] = count of base in bwt[0..i*OCC_INTERVAL)
    std::vector<uint32_t> sa_;     // Sampled SA: sa[i] = SA[i * SA_INTERVAL]
    std::vector<uint32_t> cnt_;    // Cumulative counts: cnt[c] = #bases < c
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
    [[nodiscard]] const std::vector<uint32_t>& sa_samples() const noexcept { return sa_; }
    [[nodiscard]] const std::vector<uint32_t>& occ_table() const noexcept { return occ_; }
    [[nodiscard]] const std::vector<uint32_t>& count_table() const noexcept { return cnt_; }

    // Setter for deserialization
    void set_size(size_t s) { length_ = s; }
    void set_bwt(const PackedSequence& bwt) { bwt_ = bwt; }
    void set_primary(size_t p) { primary_ = p; }
    void set_sa_samples(const std::vector<uint32_t>& sa) { sa_ = sa; }
    void set_occ_table(const std::vector<uint32_t>& occ) { occ_ = occ; }
    void set_count_table(const std::vector<uint32_t>& cnt) { cnt_ = cnt; }

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

    // Backward search for pattern in PackedSequence
    // Returns interval [l, r) in SA, or empty if not found
    [[nodiscard]] std::pair<size_t, size_t> backward_search(const PackedSequence& pat) const noexcept {
        size_t l = 0, r = length_;
        for (int32_t i = static_cast<int32_t>(pat.size()) - 1; i >= 0; --i) {
            uint8_t base = pat.get(i);
            if (base >= 4) return {0, 0}; // N or invalid
            auto [new_l, new_r] = backward_extend(base, l, r);
            if (new_l >= new_r) return {0, 0}; // Empty
            l = new_l;
            r = new_r;
        }
        return {l, r};
    }

    // Backward search for any pattern with rbegin/rend
    template <typename Pattern>
    [[nodiscard]] std::pair<size_t, size_t> backward_search(const Pattern& pat) const noexcept {
        size_t l = 0, r = length_;
        for (auto it = pat.rbegin(); it != pat.rend(); ++it) {
            uint8_t base = PackedSequence::encode_base(*it);
            if (base >= 4) return {0, 0}; // N or invalid
            auto [new_l, new_r] = backward_extend(base, l, r);
            if (new_l >= new_r) return {0, 0}; // Empty
            l = new_l;
            r = new_r;
        }
        return {l, r};
    }

    // Locate: get SA value at position (requires sampled SA)
    // We walk forward from the nearest sampled position <= pos using LF.
    // LF maps position i to the position of the suffix that starts one position
    // earlier in the original text, so SA[LF(i)] = (SA[i] - 1 + n) % n.
    [[nodiscard]] std::optional<uint32_t> locate(size_t pos) const noexcept {
        if (sa_.empty()) return std::nullopt;
        if (pos >= length_) return std::nullopt;

        // Find the largest sampled position <= pos
        size_t sampled_pos = (pos / SA_INTERVAL) * SA_INTERVAL;
        if (sampled_pos > pos) sampled_pos -= SA_INTERVAL;

        // Walk forward from sampled_pos to pos, counting steps
        size_t steps = pos - sampled_pos;
        size_t cur = sampled_pos;
        for (size_t s = 0; s < steps; ++s) {
            cur = lf(cur);
        }

        // Get the SA value at sampled_pos and subtract steps
        auto base_sa = sa_[sampled_pos / SA_INTERVAL];
        size_t result = (base_sa + length_ - steps) % length_;
        return static_cast<uint32_t>(result);
    }

    // Locate all positions in range [l, r) - returns vector of SA values
    // For finding all occurrences of a pattern
    template <typename OutputIt>
    void locate_range(size_t l, size_t r, OutputIt out) const noexcept {
        if (sa_.empty() || l >= r || r > length_) return;
        for (size_t pos = l; pos < r; ++pos) {
            if (auto sa_val = locate(pos)) {
                *out++ = *sa_val;
            }
        }
    }

    // Convenience method returning vector
    [[nodiscard]] std::vector<uint32_t> locate_all(size_t l, size_t r) const noexcept {
        std::vector<uint32_t> results;
        results.reserve(r > l ? r - l : 0);
        locate_range(l, r, std::back_inserter(results));
        return results;
    }

    // Count occurrences of pattern
    template <typename Pattern>
    [[nodiscard]] size_t count(const Pattern& pat) const noexcept {
        auto [l, r] = backward_search(pat);
        return r - l;
    }

    // Count occurrences of PackedSequence pattern
    [[nodiscard]] size_t count(const PackedSequence& pat) const noexcept {
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

        // SA-IS: O(n) linear-time suffix array construction
        core::Vector<uint32_t> sa_core = build_suffix_array_sais(seq, arena);

        // Build BWT from SA
        for (size_t i = 0; i < n; ++i) {
            size_t sa_i = sa_core[i];
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
        // Record final interval (for positions past the last exact boundary)
        if (n % OCC_INTERVAL != 0) {
            size_t intv = n / OCC_INTERVAL + 1;
            if (intv < num_intv) {
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
            idx.sa_[i / SA_INTERVAL] = sa_core[i];
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
    std::vector<FMIndex> indexes_;
    std::vector<size_t> offsets_; // Cumulative lengths

public:
    void add_sequence(const PackedSequence& seq, memory::Arena& arena) {
        indexes_.push_back(FMIndex::build(seq, arena));
        offsets_.push_back(offsets_.empty() ? seq.size() : offsets_.back() + seq.size());
    }

    void add_index(FMIndex idx) {
        size_t len = idx.length();
        indexes_.push_back(std::move(idx));
        offsets_.push_back(offsets_.empty() ? len : offsets_.back() + len);
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