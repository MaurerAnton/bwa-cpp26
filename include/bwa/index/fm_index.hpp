#pragma once

#include <span>
#include <cstddef>
#include <cstdint>
#include <array>
#include <bit>
#include <algorithm>
#include <numeric>
#include <stdexcept>

#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/sort.hpp>
#include <bwa/index/sais.hpp>

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
// Sorts all n+1 suffixes INCLUDING the empty suffix (value n), which sorts
// first. Any replacement (e.g. SA-IS) must preserve this: row 0 has SA
// value n, and codes must follow sentinel=0 < N=1 < A=2 < C=3 < G=4 < T=5.
inline core::Vector<uint32_t> build_suffix_array_brute(const PackedSequence& seq,
                                                        memory::Arena& arena) {
    size_t n = seq.size();

    // Convert to byte array with sentinel.
    // Order: sentinel=0 < N=1 < A=2 < C=3 < G=4 < T=5 (N sorts first so its
    // rows form one contiguous block; see FMIndex docs).
    core::Vector<uint8_t> T(&arena);
    T.resize(n + 1);
    T[n] = 0;
    for (size_t i = 0; i < n; ++i) {
        uint8_t b = seq.get(i);
        T[i] = (b == 4) ? 1 : (b + 2);
    }

    // Suffix starting positions: 0, 1, ..., n (n = empty suffix, smallest)
    core::Vector<uint32_t> sa(&arena);
    sa.resize(n + 1);
    for (uint32_t i = 0; i <= n; ++i) sa[i] = i;

    // Sort by direct comparison over T[0..n]
    std::sort(sa.begin(), sa.end(),
              [T_data = T.data(), n](uint32_t a, uint32_t b) {
                  size_t i = a, j = b;
                  while (i <= n && j <= n && T_data[i] == T_data[j]) { ++i; ++j; }
                  if (i > n && j > n) return a < b;
                  if (i > n) return true;   // a exhausted: smaller
                  if (j > n) return false;  // b exhausted: smaller
                  return T_data[i] < T_data[j];
              });

    return sa;
}

// Prefix-doubling suffix sort - O(n log n) via radix-sort passes.
// Production SA builder: exactly as correct as brute force, far faster
// (brute force is O(n^2 log n) with O(n) comparisons). Sorts all n+1
// suffixes INCLUDING the empty suffix (value n, sorts first), using the
// same codes as brute force: sentinel=0 < N=1 < A=2 < C=3 < G=4 < T=5.
// Each pass sorts 64-bit combined keys (rank[i]+1)<<32 | (rank[i+k]+1)
// (0 = past-end, smallest) with a stable radix sort, then re-ranks.
// Terminates once all ranks are distinct (within ceil(log2(n+1)) passes).
inline core::Vector<uint32_t> build_suffix_array_doubling(const PackedSequence& seq,
                                                          memory::Arena& arena) {
    size_t n = seq.size();
    size_t m = n + 1;

    core::Vector<uint32_t> rank(&arena), newrank(&arena);
    core::Vector<uint64_t> keys(&arena), key_buf(&arena);
    core::Vector<uint32_t> vals(&arena), val_buf(&arena);
    rank.resize(m);
    newrank.resize(m);
    keys.resize(m);
    key_buf.resize(m);
    vals.resize(m);
    val_buf.resize(m);

    for (size_t i = 0; i < n; ++i) {
        uint8_t b = seq.get(i);
        rank[i] = (b == 4) ? 1 : (b + 2);
        vals[i] = static_cast<uint32_t>(i);
    }
    rank[n] = 0;
    vals[n] = static_cast<uint32_t>(n);

    for (size_t k = 1; k < m; k <<= 1) {
        for (size_t i = 0; i < m; ++i) {
            uint32_t r1 = rank[i] + 1;
            uint32_t r2 = (i + k <= n) ? rank[i + k] + 1 : 0;
            keys[i] = (static_cast<uint64_t>(r1) << 32) | r2;
            vals[i] = static_cast<uint32_t>(i);
        }
        core::detail::radix_sort_pairs(keys.data(), vals.data(), m, key_buf, val_buf);
        // keys[] sorted, vals[] carries the SA order; re-rank by key runs
        newrank[vals[0]] = 0;
        uint32_t r = 0;
        for (size_t t = 1; t < m; ++t) {
            if (keys[t] != keys[t - 1]) ++r;
            newrank[vals[t]] = r;
        }
        rank.swap(newrank);
        if (r + 1 >= m) break;  // all ranks distinct
    }
    return vals;
}

// Permutation check for a candidate suffix array: nrows entries holding
// each value in [0, nrows) exactly once. Guards the SA-IS fast path.
inline bool is_valid_sa(const core::Vector<uint32_t>& sa, size_t nrows) {
    if (sa.size() != nrows || nrows == 0) return false;
    std::vector<char> seen(nrows, 0);
    for (size_t i = 0; i < nrows; ++i) {
        if (sa[i] >= nrows || seen[sa[i]]) return false;
        seen[sa[i]] = 1;
    }
    return true;
}

// FM-index with rank/select support (standard construction with explicit
// sentinel). Rows = text length + 1; SA includes the empty suffix (value n,
// always row 0). BWT codes match suffix-sort order: $=0 < N=1 < A=2 < C=3 <
// G=4 < T=5. Query/text bases use A=0..N=4 and are mapped at the search
// boundary; query Ns never match (rejected in backward_search), so reference
// Ns act as natural separators. LF is exact on every row (single cycle),
// hence locate() walks need no special cases.
class FMIndex {
public:
    using occ_t = uint32_t;
    static constexpr int OCC_INTERVAL = 128; // Rank sampling interval
    static constexpr int SA_INTERVAL = 32;   // SA sampling interval
    static constexpr int ALPHABET = 6;       // $, A, C, G, T, N
    static constexpr uint8_t CODE_SENTINEL = 0;

    // Query/text base (A=0..N=4) -> BWT code ($=0, N=1, A=2, C=3, G=4, T=5)
    [[nodiscard]] static constexpr uint8_t to_code(uint8_t base) noexcept {
        // base: A=0, C=1, G=2, T=3, N=4
        // BWT codes: $=0, N=1, A=2, C=3, G=4, T=5
        return base == 4 ? 1 : static_cast<uint8_t>(base + 2);
    }

private:
    std::vector<uint8_t> bwt_;   // BWT sort codes, length_ + 1 cells
    std::vector<occ_t> occ_;       // Rank table: occ[code][i] = count of code in bwt[0..i*OCC_INTERVAL)
    std::vector<uint32_t> sa_;     // Sampled SA: sa[i] = SA[i * SA_INTERVAL]
    std::vector<uint32_t> cnt_;    // Cumulative counts: cnt[c] = #codes < c
    size_t primary_ = 0;      // Row holding the sentinel cell (SA[i] == 0)
    size_t length_ = 0;       // Text length (rows = length_ + 1)

    // Occurrence array layout: [6][num_intervals]
    [[nodiscard]] size_t occ_offset(uint8_t code, size_t interval) const noexcept {
        return interval * ALPHABET + code;
    }

    [[nodiscard]] size_t num_intervals() const noexcept {
        return (rows() + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
    }

public:
    // Row count (text length + 1 sentinel row). Public: search code must
    // initialize intervals to [0, rows()), not [0, length()).
    [[nodiscard]] size_t rows() const noexcept { return length_ + 1; }

    FMIndex() = default;

    [[nodiscard]] size_t size() const noexcept { return length_; }
    [[nodiscard]] size_t length() const noexcept { return length_; }
    [[nodiscard]] size_t bwt_size() const noexcept { return bwt_.size(); }
    [[nodiscard]] const std::vector<uint8_t>& bwt() const noexcept { return bwt_; }
    [[nodiscard]] size_t primary() const noexcept { return primary_; }

    // Serialization accessors
    [[nodiscard]] const std::vector<uint32_t>& sa_samples() const noexcept { return sa_; }
    [[nodiscard]] const std::vector<uint32_t>& occ_table() const noexcept { return occ_; }
    [[nodiscard]] const std::vector<uint32_t>& count_table() const noexcept { return cnt_; }

    // Setter for deserialization
    void set_size(size_t s) { length_ = s; }
    void set_bwt(const std::vector<uint8_t>& bwt) { bwt_ = bwt; }
    void set_primary(size_t p) { primary_ = p; }
    void set_sa_samples(const std::vector<uint32_t>& sa) { sa_ = sa; }
    void set_occ_table(const std::vector<uint32_t>& occ) { occ_ = occ; }
    void set_count_table(const std::vector<uint32_t>& cnt) { cnt_ = cnt; }

    // Rank: number of occurrences of code in bwt[0..pos)
    [[nodiscard]] uint32_t rank(uint8_t code, size_t pos) const noexcept {
        if (pos == 0) return 0;
        size_t interval = pos / OCC_INTERVAL;
        size_t offset = pos % OCC_INTERVAL;
        (void)offset;

        uint32_t count = occ_[occ_offset(code, interval)];
        // Scan remaining
        size_t start = interval * OCC_INTERVAL;
        for (size_t i = start; i < pos; ++i) {
            if (bwt_[i] == code) ++count;
        }
        return count;
    }

    // C array: cumulative counts (code order $ < A < C < G < T < N)
    [[nodiscard]] uint32_t C(uint8_t code) const noexcept {
        return cnt_[code];
    }

    // Total count of code
    [[nodiscard]] uint32_t total(uint8_t code) const noexcept {
        return cnt_[code + 1] - cnt_[code];
    }

    // LF mapping: LF(pos) = C[bwt[pos]] + rank(bwt[pos], pos)
    [[nodiscard]] size_t lf(size_t pos) const noexcept {
        uint8_t code = bwt_[pos];
        return C(code) + rank(code, pos);
    }

    // Backward search: extend interval [l, r) with a query base (A=0..T=3)
    // Returns new interval [l', r')
    [[nodiscard]] std::pair<size_t, size_t> backward_extend(uint8_t base,
                                                             size_t l, size_t r) const noexcept {
        uint8_t code = to_code(base);
        size_t new_l = C(code) + rank(code, l);
        size_t new_r = C(code) + rank(code, r);
        return {new_l, new_r};
    }

    // Backward search for pattern in PackedSequence
    // Returns interval [l, r) in SA, or empty if not found
    [[nodiscard]] std::pair<size_t, size_t> backward_search(const PackedSequence& pat) const noexcept {
        size_t l = 0, r = rows();
        for (int32_t i = static_cast<int32_t>(pat.size()) - 1; i >= 0; --i) {
            uint8_t base = pat.get(i);
            if (base >= 4) return {0, 0}; // N or invalid never matches
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
        size_t l = 0, r = rows();
        for (auto it = pat.rbegin(); it != pat.rend(); ++it) {
            uint8_t base = PackedSequence::encode_base(*it);
            if (base >= 4) return {0, 0}; // N or invalid never matches
            auto [new_l, new_r] = backward_extend(base, l, r);
            if (new_l >= new_r) return {0, 0}; // Empty
            l = new_l;
            r = new_r;
        }
        return {l, r};
    }

    // Locate: get SA value at row pos via the sampled SA.
    // LF maps row i to the row of the suffix starting one position earlier:
    // SA[LF(i)] = (SA[i]-1+n+1)%(n+1), hence SA[i] = (SA[LF(i)]+1)%(n+1).
    // Walk LF from pos until a sampled row, counting steps k, then
    // Locate: get SA value at row pos via the sampled SA.
    // LF maps row i to the row of the suffix starting one position earlier:
    // SA[LF(i)] = (SA[i]-1+n+1)%(n+1), hence SA[i] = (SA[LF(i)]+1)%(n+1).
    // Walk LF from pos until a row with known SA value: an explicitly sampled
    // text row (SA indices 1, 1+SA_INTERVAL, ...) or row 0, whose SA value is
    // implicitly n (empty suffix). SA[pos] = (known + steps) % (n+1).
    // Row 0 itself is rejected (not a text position); search intervals never
    // contain it since every extension starts at C[code] >= 1.
    [[nodiscard]] std::optional<uint32_t> locate(size_t pos) const noexcept {
        if (sa_.empty() || length_ == 0) return std::nullopt;
        // Text rows are 1..length_ (row 0 is the empty suffix)
        if (pos == 0 || pos >= rows()) return std::nullopt;

        size_t cur = pos;
        size_t steps = 0;
        uint32_t base = 0;
        bool have_base = false;
        // lf() permutes [0, n+1) in a single cycle: a known row is reached
        // within n+1 steps; the cap guards against corrupt tables.
        while (steps <= rows()) {
            if (cur == 0) {
                base = static_cast<uint32_t>(length_);  // SA[row 0] = n
                have_base = true;
                break;
            }
            if ((cur - 1) % SA_INTERVAL == 0) {
                size_t samp = (cur - 1) / SA_INTERVAL;
                if (samp >= sa_.size()) return std::nullopt;
                base = sa_[samp];
                if (base > length_) return std::nullopt;
                have_base = true;
                break;
            }
            cur = lf(cur);
            ++steps;
        }
        if (!have_base) return std::nullopt;
        uint32_t result = (base + static_cast<uint32_t>(steps)) %
                          static_cast<uint32_t>(rows());
        // Text positions are 0..length_-1
        if (result >= length_) return std::nullopt;
        return result;
    }

    // Locate all positions in range [l, r) - returns vector of SA values
    // For finding all occurrences of a pattern (rows span [0, n+1))
    template <typename OutputIt>
    void locate_range(size_t l, size_t r, OutputIt out) const noexcept {
        if (sa_.empty() || l >= r || r > rows()) return;
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

    // Build from packed sequence (SA construction, rows = n + 1)
    static FMIndex build(const PackedSequence& seq, memory::Arena& arena) {
        FMIndex idx;
        idx.length_ = seq.size();
        idx.bwt_.resize(seq.size() + 1);

        // Allocate working memory from arena
        size_t n = seq.size();
        size_t nrows = n + 1;

        // SA construction: SA-IS O(n) is the production path (rewritten and
        // differentially tested against brute force). Validation keeps the
        // prefix-doubling O(n log n) and brute-force fallbacks for safety.
        core::Vector<uint32_t> sa_core = detail::sais::build_suffix_array(seq, arena);
        if (!is_valid_sa(sa_core, nrows)) {
            sa_core = build_suffix_array_doubling(seq, arena);
        }
        if (!is_valid_sa(sa_core, nrows)) {
            sa_core = build_suffix_array_brute(seq, arena);
        }

        // Build BWT from SA: BWT[i] = text[(SA[i]-1) mod (n+1)] in sort codes
        // ($=0, N=1, A=2, C=3, G=4, T=5). The row of suffix 0 holds the
        // sentinel (its predecessor is the empty suffix); the row of the
        // empty suffix (SA[i]==n) holds text[n-1].
        for (size_t i = 0; i < nrows; ++i) {
            size_t sa_i = sa_core[i];
            uint8_t code;
            if (sa_i == 0) {
                code = CODE_SENTINEL;
            } else {
                uint8_t b = seq.get(sa_i - 1);
                code = to_code(b);
            }
            idx.bwt_[i] = code;
            if (sa_i == 0) idx.primary_ = i;
        }

        // Build occ table
        size_t num_intv = idx.num_intervals();
        idx.occ_.resize(num_intv * ALPHABET);
        idx.cnt_.resize(ALPHABET + 1); // cnt[0..6], cnt[6] = n + 1

        // Count codes (all 6 classes)
        std::array<uint32_t, ALPHABET> counts = {0, 0, 0, 0, 0, 0};
        for (size_t i = 0; i < nrows; ++i) {
            uint8_t code = idx.bwt_[i];
            if (code < ALPHABET) ++counts[code];

            if ((i + 1) % OCC_INTERVAL == 0) {
                size_t intv = (i + 1) / OCC_INTERVAL;
                for (int b = 0; b < ALPHABET; ++b) {
                    idx.occ_[idx.occ_offset(b, intv)] = counts[b];
                }
            }
        }
        // Record final interval (for positions past the last exact boundary)
        if (nrows % OCC_INTERVAL != 0) {
            size_t intv = nrows / OCC_INTERVAL + 1;
            if (intv < num_intv) {
                for (int b = 0; b < ALPHABET; ++b) {
                    idx.occ_[idx.occ_offset(b, intv)] = counts[b];
                }
            }
        }

        // Cumulative counts
        idx.cnt_[0] = 0;
        for (int b = 0; b < ALPHABET; ++b) {
            idx.cnt_[b + 1] = idx.cnt_[b] + counts[b];
        }

        // Sample SA at text positions (SA array indices 1, 1+SA_INTERVAL, 1+2*SA_INTERVAL...)
// These correspond to text positions in suffix array order, excluding sentinel row 0.
        size_t n_text = nrows - 1;
        idx.sa_.resize((n_text + SA_INTERVAL - 1) / SA_INTERVAL);
        for (size_t i = 1; i < nrows; i += SA_INTERVAL) {
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