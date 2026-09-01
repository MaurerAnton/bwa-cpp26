#pragma once

#include <bwa/core/arena.hpp>
#include <string_view>
#include <span>
#include <charconv>
#include <format>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <expected>

namespace bwa::core {

// PMR-compatible string with SSO and kstring_t-compatible API
// Inline buffer: 224 bytes (3.5 cache lines) - covers most BWA strings (read names, CIGAR, etc.)
class PmrString {
public:
    using value_type = char;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using allocator_type = memory::PolymorphicAllocator<char>;
    using reference = char&;
    using const_reference = const char&;
    using pointer = char*;
    using const_pointer = const char*;
    using iterator = char*;
    using const_iterator = const char*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    static constexpr size_type INLINE_CAPACITY = 224;
    static constexpr size_type NPOS = static_cast<size_type>(-1);

private:
    union Storage {
        alignas(16) char inline_buffer[INLINE_CAPACITY];
        char* heap_ptr;
    };

    Storage storage_;
    size_type size_ = 0;
    size_type capacity_ = INLINE_CAPACITY;
    [[no_unique_address]] allocator_type alloc_;

    [[nodiscard]] bool is_inline() const noexcept {
        return capacity_ == INLINE_CAPACITY;
    }

    [[nodiscard]] pointer data_ptr() noexcept {
        return is_inline() ? storage_.inline_buffer : storage_.heap_ptr;
    }

    [[nodiscard]] const_pointer data_ptr() const noexcept {
        return is_inline() ? storage_.inline_buffer : storage_.heap_ptr;
    }

    void allocate_heap(size_type new_cap) {
        if (new_cap == 0) {
            storage_.heap_ptr = nullptr;
            capacity_ = 0;
            return;
        }
        storage_.heap_ptr = std::allocator_traits<allocator_type>::allocate(alloc_, new_cap);
        capacity_ = new_cap;
    }

    void deallocate_heap() noexcept {
        if (!is_inline() && storage_.heap_ptr) {
            std::allocator_traits<allocator_type>::deallocate(alloc_, storage_.heap_ptr, capacity_);
            storage_.heap_ptr = nullptr;
            capacity_ = 0;
        }
    }

    void grow(size_type min_capacity) {
        size_type new_cap = std::max(capacity_ == 0 ? 1 : capacity_ * 2, min_capacity);
        if (new_cap <= INLINE_CAPACITY) new_cap = INLINE_CAPACITY + 1;

        pointer new_data = std::allocator_traits<allocator_type>::allocate(alloc_, new_cap);
        if (size_ > 0) {
            std::memcpy(new_data, data_ptr(), size_);
        }
        deallocate_heap();
        storage_.heap_ptr = new_data;
        capacity_ = new_cap;
    }

    void ensure_capacity(size_type additional) {
        if (size_ + additional > capacity_) {
            grow(size_ + additional);
        }
    }

public:
    // Constructors
    PmrString() noexcept = default;

    explicit PmrString(const allocator_type& alloc) noexcept
        : alloc_(alloc) {}

    PmrString(const PmrString& other)
        : alloc_(std::allocator_traits<allocator_type>::
                 select_on_container_copy_construction(other.alloc_)) {
        if (other.size_ > 0) {
            ensure_capacity(other.size_);
            std::memcpy(data_ptr(), other.data_ptr(), other.size_);
            size_ = other.size_;
        }
    }

    PmrString(PmrString&& other) noexcept
        : alloc_(std::move(other.alloc_))
        , size_(other.size_)
        , capacity_(other.capacity_) {
        if (other.is_inline()) {
            std::memcpy(storage_.inline_buffer, other.storage_.inline_buffer, size_);
        } else {
            storage_.heap_ptr = other.storage_.heap_ptr;
        }
        other.size_ = 0;
        other.capacity_ = INLINE_CAPACITY;
    }

    PmrString(std::string_view sv, const allocator_type& alloc = allocator_type{})
        : alloc_(alloc) {
        assign(sv);
    }

    PmrString(const char* s, const allocator_type& alloc = allocator_type{})
        : alloc_(alloc) {
        assign(s);
    }

    PmrString(const char* s, size_type n, const allocator_type& alloc = allocator_type{})
        : alloc_(alloc) {
        assign(s, n);
    }

    template <std::input_iterator It>
    PmrString(It first, It last, const allocator_type& alloc = allocator_type{})
        : alloc_(alloc) {
        append(first, last);
    }

    ~PmrString() {
        deallocate_heap();
    }

    // Assignment
    PmrString& operator=(const PmrString& other) {
        if (this != &other) {
            clear();
            if (std::allocator_traits<allocator_type>::propagate_on_container_copy_assignment::value) {
                alloc_ = other.alloc_;
            }
            assign(other.data(), other.size_);
        }
        return *this;
    }

    PmrString& operator=(PmrString&& other) noexcept(
        std::allocator_traits<allocator_type>::propagate_on_container_move_assignment::value ||
        std::allocator_traits<allocator_type>::is_always_equal::value) {
        if (this != &other) {
            clear();
            deallocate_heap();

            if (std::allocator_traits<allocator_type>::propagate_on_container_move_assignment::value) {
                alloc_ = std::move(other.alloc_);
            }

            size_ = other.size_;
            capacity_ = other.capacity_;

            if (other.is_inline()) {
                std::memcpy(storage_.inline_buffer, other.storage_.inline_buffer, size_);
            } else {
                storage_.heap_ptr = other.storage_.heap_ptr;
            }

            other.size_ = 0;
            other.capacity_ = INLINE_CAPACITY;
        }
        return *this;
    }

    PmrString& operator=(std::string_view sv) {
        assign(sv);
        return *this;
    }

    PmrString& operator=(const char* s) {
        assign(s);
        return *this;
    }

    // Iterators
    [[nodiscard]] iterator begin() noexcept { return data(); }
    [[nodiscard]] const_iterator begin() const noexcept { return data(); }
    [[nodiscard]] const_iterator cbegin() const noexcept { return data(); }
    [[nodiscard]] iterator end() noexcept { return data() + size_; }
    [[nodiscard]] const_iterator end() const noexcept { return data() + size_; }
    [[nodiscard]] const_iterator cend() const noexcept { return data() + size_; }
    [[nodiscard]] reverse_iterator rbegin() noexcept { return reverse_iterator(end()); }
    [[nodiscard]] const_reverse_iterator rbegin() const noexcept { return const_reverse_iterator(end()); }
    [[nodiscard]] const_reverse_iterator crbegin() const noexcept { return const_reverse_iterator(end()); }
    [[nodiscard]] reverse_iterator rend() noexcept { return reverse_iterator(begin()); }
    [[nodiscard]] const_reverse_iterator rend() const noexcept { return const_reverse_iterator(begin()); }
    [[nodiscard]] const_reverse_iterator crend() const noexcept { return const_reverse_iterator(begin()); }

    // Capacity
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
    [[nodiscard]] size_type size() const noexcept { return size_; }
    [[nodiscard]] size_type length() const noexcept { return size_; }
    [[nodiscard]] size_type capacity() const noexcept { return capacity_; }
    [[nodiscard]] size_type max_size() const noexcept { return std::allocator_traits<allocator_type>::max_size(alloc_); }

    void reserve(size_type new_cap) {
        if (new_cap > capacity_) {
            if (new_cap <= INLINE_CAPACITY) return;
            grow(new_cap);
        }
    }

    void shrink_to_fit() {
        if (size_ == capacity_) return;
        if (size_ <= INLINE_CAPACITY) {
            if (!is_inline()) {
                std::memcpy(storage_.inline_buffer, data(), size_);
                deallocate_heap();
                capacity_ = INLINE_CAPACITY;
            }
        } else {
            pointer new_data = std::allocator_traits<allocator_type>::allocate(alloc_, size_);
            std::memcpy(new_data, data_ptr(), size_);
            deallocate_heap();
            storage_.heap_ptr = new_data;
            capacity_ = size_;
        }
    }

    void clear() noexcept { size_ = 0; }

    // Element access
    [[nodiscard]] reference operator[](size_type pos) noexcept { return data_ptr()[pos]; }
    [[nodiscard]] const_reference operator[](size_type pos) const noexcept { return data_ptr()[pos]; }
    [[nodiscard]] reference at(size_type pos) {
        if (pos >= size_) throw std::out_of_range("PmrString::at");
        return data_ptr()[pos];
    }
    [[nodiscard]] const_reference at(size_type pos) const {
        if (pos >= size_) throw std::out_of_range("PmrString::at");
        return data_ptr()[pos];
    }
    [[nodiscard]] reference front() noexcept { return data_ptr()[0]; }
    [[nodiscard]] const_reference front() const noexcept { return data_ptr()[0]; }
    [[nodiscard]] reference back() noexcept { return data_ptr()[size_ - 1]; }
    [[nodiscard]] const_reference back() const noexcept { return data_ptr()[size_ - 1]; }
    [[nodiscard]] pointer data() noexcept { return data_ptr(); }
    [[nodiscard]] const_pointer data() const noexcept { return data_ptr(); }
    [[nodiscard]] pointer c_str() noexcept {
        // Ensure null-termination (not counted in size)
        if (size_ == capacity_) {
            ensure_capacity(1);
        }
        data_ptr()[size_] = '\0';
        return data_ptr();
    }

    [[nodiscard]] const_pointer c_str() const noexcept {
        // For const version, just return data (may not be null-terminated)
        return data_ptr();
    }

    // String view interop
    [[nodiscard]] std::string_view view() const noexcept { return {data_ptr(), size_}; }
    [[nodiscard]] operator std::string_view() const noexcept { return view(); }

    // Modifiers - kstring_t compatible API
    void kputc(char c) {
        ensure_capacity(1);
        data_ptr()[size_++] = c;
    }

    void kputs(const char* s) {
        if (!s) return;
        size_type n = std::strlen(s);
        ensure_capacity(n);
        std::memcpy(data_ptr() + size_, s, n);
        size_ += n;
    }

    void kputsn(const char* s, size_type n) {
        if (!s || n == 0) return;
        ensure_capacity(n);
        std::memcpy(data_ptr() + size_, s, n);
        size_ += n;
    }

    void kputw(int32_t x) {
        // Reserve max digits for int32_t: -2147483648 (11 chars)
        ensure_capacity(12);
        auto [ptr, ec] = std::to_chars(data_ptr() + size_, data_ptr() + capacity_, x);
        if (ec == std::errc{}) {
            size_ = ptr - data_ptr();
        }
    }

    void kputuw(uint32_t x) {
        ensure_capacity(11); // 4294967295 (10 digits)
        auto [ptr, ec] = std::to_chars(data_ptr() + size_, data_ptr() + capacity_, x);
        if (ec == std::errc{}) {
            size_ = ptr - data_ptr();
        }
    }

    void kputl(int64_t x) {
        ensure_capacity(21); // -9223372036854775808 (20 digits)
        auto [ptr, ec] = std::to_chars(data_ptr() + size_, data_ptr() + capacity_, x);
        if (ec == std::errc{}) {
            size_ = ptr - data_ptr();
        }
    }

    void kputul(uint64_t x) {
        ensure_capacity(20); // 18446744073709551615 (20 digits)
        auto [ptr, ec] = std::to_chars(data_ptr() + size_, data_ptr() + capacity_, x);
        if (ec == std::errc{}) {
            size_ = ptr - data_ptr();
        }
    }

    void kputd(double x) {
        // %g style - shortest representation
        ensure_capacity(32);
        auto [ptr, ec] = std::to_chars(data_ptr() + size_, data_ptr() + capacity_, x,
                                        std::chars_format::general);
        if (ec == std::errc{}) {
            size_ = ptr - data();
        }
    }

    void kputf(double x, int precision) {
        // %f style - fixed precision
        ensure_capacity(32 + precision);
        auto [ptr, ec] = std::to_chars(data() + size_, data() + capacity_, x,
                                        std::chars_format::fixed, precision);
        if (ec == std::errc{}) {
            size_ = ptr - data();
        }
    }

    // ksprintf - fmt-style formatting
    template <typename... Args>
    int ksprintf(std::format_string<Args...> fmt, Args&&... args) {
        try {
            size_type old_size = size_;
            // Reserve heuristic: format string length + 50% for args
            // format_string doesn't have .str() in C++26, use get()
            std::string_view fmt_view = fmt.get();
            ensure_capacity(fmt_view.size() + fmt_view.size() / 2 + 32);
            
            auto out = std::format_to(data() + size_, fmt, std::forward<Args>(args)...);
            size_ = out - data();
            return static_cast<int>(size_ - old_size);
        } catch (...) {
            return -1;
        }
    }

    // Modern C++ API
    void push_back(char c) { kputc(c); }
    void pop_back() noexcept { if (size_ > 0) --size_; }

    PmrString& append(std::string_view sv) {
        ensure_capacity(sv.size());
        std::memcpy(data() + size_, sv.data(), sv.size());
        size_ += sv.size();
        return *this;
    }

    PmrString& append(const char* s) { return append(std::string_view{s}); }
    PmrString& append(const char* s, size_type n) { return append(std::string_view{s, n}); }
    PmrString& append(const PmrString& other) { return append(other.view()); }

    template <std::input_iterator It>
    PmrString& append(It first, It last) {
        for (; first != last; ++first) {
            push_back(*first);
        }
        return *this;
    }

    PmrString& assign(std::string_view sv) {
        clear();
        return append(sv);
    }

    PmrString& assign(const char* s) { return assign(std::string_view{s}); }
    PmrString& assign(const char* s, size_type n) { return assign(std::string_view{s, n}); }

    void resize(size_type n, char c = '\0') {
        if (n < size_) {
            size_ = n;
        } else if (n > size_) {
            ensure_capacity(n - size_);
            std::fill_n(data() + size_, n - size_, c);
            size_ = n;
        }
    }

    void swap(PmrString& other) noexcept(
        std::allocator_traits<allocator_type>::propagate_on_container_swap::value ||
        std::allocator_traits<allocator_type>::is_always_equal::value) {
        if (this != &other) {
            if (std::allocator_traits<allocator_type>::propagate_on_container_swap::value) {
                std::swap(alloc_, other.alloc_);
            }
            std::swap(size_, other.size_);
            std::swap(capacity_, other.capacity_);

            if (is_inline() && other.is_inline()) {
                char temp[INLINE_CAPACITY];
                std::memcpy(temp, storage_.inline_buffer, std::max(size_, other.size_));
                std::memcpy(storage_.inline_buffer, other.storage_.inline_buffer, other.size_);
                std::memcpy(other.storage_.inline_buffer, temp, size_);
            } else if (!is_inline() && !other.is_inline()) {
                std::swap(storage_.heap_ptr, other.storage_.heap_ptr);
            } else if (is_inline()) {
                char temp[INLINE_CAPACITY];
                std::memcpy(temp, storage_.inline_buffer, size_);
                storage_.heap_ptr = other.storage_.heap_ptr;
                other.capacity_ = INLINE_CAPACITY;
                std::memcpy(other.storage_.inline_buffer, temp, size_);
            } else {
                char temp[INLINE_CAPACITY];
                std::memcpy(temp, other.storage_.inline_buffer, other.size_);
                other.storage_.heap_ptr = storage_.heap_ptr;
                capacity_ = INLINE_CAPACITY;
                std::memcpy(storage_.inline_buffer, temp, other.size_);
            }
        }
    }

    // Search/find (subset)
    [[nodiscard]] size_type find(char c, size_type pos = 0) const noexcept {
        if (pos >= size_) return NPOS;
        const_pointer p = std::find(data() + pos, data() + size_, c);
        return p == data() + size_ ? NPOS : p - data();
    }

    [[nodiscard]] size_type find(std::string_view sv, size_type pos = 0) const noexcept {
        if (pos >= size_) return NPOS;
        const_pointer p = std::search(data() + pos, data() + size_, sv.begin(), sv.end());
        return p == data() + size_ ? NPOS : p - data();
    }

    [[nodiscard]] size_type rfind(char c, size_type pos = NPOS) const noexcept {
        if (size_ == 0) return NPOS;
        if (pos >= size_) pos = size_ - 1;
        for (size_type i = pos; i != NPOS; --i) {
            if (data()[i] == c) return i;
        }
        return NPOS;
    }

    [[nodiscard]] PmrString substr(size_type pos = 0, size_type len = NPOS) const {
        if (pos > size_) throw std::out_of_range("PmrString::substr");
        size_type rlen = std::min(len, size_ - pos);
        return PmrString{data() + pos, rlen, alloc_};
    }

    // Comparison
    friend bool operator==(const PmrString& lhs, const PmrString& rhs) noexcept {
        return lhs.view() == rhs.view();
    }
    friend bool operator==(const PmrString& lhs, std::string_view rhs) noexcept {
        return lhs.view() == rhs;
    }
    friend bool operator==(std::string_view lhs, const PmrString& rhs) noexcept {
        return lhs == rhs.view();
    }

    friend auto operator<=>(const PmrString& lhs, const PmrString& rhs) noexcept {
        return lhs.view() <=> rhs.view();
    }

    [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc_; }

    // kstring_t compatibility: direct access to internals
    [[nodiscard]] size_type& l() noexcept { return size_; }
    [[nodiscard]] const size_type& l() const noexcept { return size_; }
    [[nodiscard]] size_type& m() noexcept { return capacity_; }
    [[nodiscard]] const size_type& m() const noexcept { return capacity_; }
    [[nodiscard]] pointer& s() noexcept { return storage_.heap_ptr; }
    [[nodiscard]] const_pointer s() const noexcept { return data_ptr(); }
};

// Non-member swap
inline void swap(PmrString& a, PmrString& b) noexcept { a.swap(b); }

// kstring_t compatibility functions (for gradual translation)
inline void ks_init(PmrString& s) { s.clear(); }
inline void ks_free(PmrString& s) { s.clear(); s.shrink_to_fit(); }
inline int ks_resize(PmrString& s, std::size_t n) {
    if (n > s.capacity()) s.reserve(n);
    return 0;
}
inline int kputc(char c, PmrString& s) { s.kputc(c); return 0; }
inline int kputs(const char* str, PmrString& s) { s.kputs(str); return 0; }
inline int kputsn(const char* str, std::size_t n, PmrString& s) { s.kputsn(str, n); return 0; }
inline int kputw(int x, PmrString& s) { s.kputw(x); return 0; }
inline int kputuw(unsigned x, PmrString& s) { s.kputuw(x); return 0; }
inline int kputl(long x, PmrString& s) { s.kputl(x); return 0; }
inline int kputul(unsigned long x, PmrString& s) { s.kputul(x); return 0; }
inline int kputd(double x, PmrString& s) { s.kputd(x); return 0; }
inline int kputf(double x, int prec, PmrString& s) { s.kputf(x, prec); return 0; }

// Type aliases
using KString = PmrString; // Direct kstring_t replacement

} // namespace bwa::core

// Hash support
namespace std {
template <>
struct hash<bwa::core::PmrString> {
    size_t operator()(const bwa::core::PmrString& s) const noexcept {
        return std::hash<std::string_view>{}(s.view());
    }
};
}