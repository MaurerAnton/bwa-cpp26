#pragma once

#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <algorithm>
#include <functional>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <concepts>
#include <execution>

namespace bwa::core {

// Sorting utilities (ksort.h replacement)
// Provides:
// - Radix sort for integers (O(n), cache-friendly)
// - Introsort for general types (std::sort wrapper)
// - Parallel sort (C++17 execution policies)
// - Key-value pair sorting

namespace detail {

// Radix sort for unsigned integers (LSD, 11 bits per pass)
template <typename T>
concept UnsignedIntegral = std::unsigned_integral<T>;

template <UnsignedIntegral T>
void radix_sort_unsigned(T* data, size_t n, Vector<T>& buffer) {
    if (n <= 1) return;
    buffer.resize(n);

    constexpr int BITS = 11;
    constexpr int MASK = (1 << BITS) - 1;
    constexpr int PASSES = (sizeof(T) * 8 + BITS - 1) / BITS;

    T* src = data;
    T* dst = buffer.data();

    for (int pass = 0, shift = 0; pass < PASSES; ++pass, shift += BITS) {
        // Counting sort
        size_t counts[1 << BITS] = {};
        for (size_t i = 0; i < n; ++i) {
            ++counts[(src[i] >> shift) & MASK];
        }

        // Prefix sum
        for (int i = 1; i < (1 << BITS); ++i) {
            counts[i] += counts[i - 1];
        }

        // Distribute (backward for stability)
        for (size_t i = n; i-- > 0;) {
            size_t idx = (src[i] >> shift) & MASK;
            dst[--counts[idx]] = src[i];
        }

        std::swap(src, dst);
    }

    // Copy back if odd number of passes
    if (src != data) {
        std::copy_n(src, n, data);
    }
}

// Radix sort for signed integers (flip sign bit)
template <std::signed_integral T>
void radix_sort_signed(T* data, size_t n, Vector<T>& buffer) {
    using UT = std::make_unsigned_t<T>;
    constexpr UT SIGN_BIT = UT(1) << (sizeof(T) * 8 - 1);

    // Flip sign bit to make negative < positive
    UT* udata = reinterpret_cast<UT*>(data);
    for (size_t i = 0; i < n; ++i) {
        udata[i] ^= SIGN_BIT;
    }

    Vector<UT> ubuf;
    ubuf.resize(n);
    radix_sort_unsigned(udata, n, ubuf);

    // Flip back
    for (size_t i = 0; i < n; ++i) {
        udata[i] ^= SIGN_BIT;
    }
}

// Radix sort for floats (IEEE 754 compliant)
template <std::floating_point T>
void radix_sort_float(T* data, size_t n, Vector<T>& buffer) {
    using UT = std::conditional_t<sizeof(T) == 4, uint32_t, uint64_t>;
    constexpr UT SIGN_BIT = UT(1) << (sizeof(T) * 8 - 1);

    UT* udata = reinterpret_cast<UT*>(data);
    for (size_t i = 0; i < n; ++i) {
        // Flip sign bit, and for negative numbers flip all bits
        UT u = udata[i];
        udata[i] = (u & SIGN_BIT) ? ~u : (u ^ SIGN_BIT);
    }

    Vector<UT> ubuf;
    ubuf.resize(n);
    radix_sort_unsigned(udata, n, ubuf);

    // Reverse transform
    for (size_t i = 0; i < n; ++i) {
        UT u = udata[i];
        udata[i] = (u & SIGN_BIT) ? (u ^ SIGN_BIT) : ~u;
    }
}

// Key-value radix sort
template <UnsignedIntegral Key, typename Value>
void radix_sort_pairs(Key* keys, Value* values, size_t n,
                       Vector<Key>& key_buf, Vector<Value>& val_buf) {
    if (n <= 1) return;
    key_buf.resize(n);
    val_buf.resize(n);

    constexpr int BITS = 11;
    constexpr int MASK = (1 << BITS) - 1;
    constexpr int PASSES = (sizeof(Key) * 8 + BITS - 1) / BITS;

    Key* ksrc = keys;
    Key* kdst = key_buf.data();
    Value* vsrc = values;
    Value* vdst = val_buf.data();

    for (int pass = 0, shift = 0; pass < PASSES; ++pass, shift += BITS) {
        size_t counts[1 << BITS] = {};
        for (size_t i = 0; i < n; ++i) {
            ++counts[(ksrc[i] >> shift) & MASK];
        }
        for (int i = 1; i < (1 << BITS); ++i) {
            counts[i] += counts[i - 1];
        }
        for (size_t i = n; i-- > 0;) {
            size_t idx = (ksrc[i] >> shift) & MASK;
            kdst[--counts[idx]] = ksrc[i];
            vdst[counts[idx]] = vsrc[i];
        }
        std::swap(ksrc, kdst);
        std::swap(vsrc, vdst);
    }

    if (ksrc != keys) {
        std::copy_n(ksrc, n, keys);
        std::copy_n(vsrc, n, values);
    }
}

} // namespace detail

// Public API

// Sort unsigned integers
template <detail::UnsignedIntegral T>
void radix_sort(T* data, size_t n) {
    Vector<T> buffer;
    detail::radix_sort_unsigned(data, n, buffer);
}

template <detail::UnsignedIntegral T, typename Alloc>
void radix_sort(Vector<T, 0, Alloc>& vec) {
    radix_sort(vec.data(), vec.size());
}

template <detail::UnsignedIntegral T, size_t N, typename Alloc>
void radix_sort(Vector<T, N, Alloc>& vec) {
    radix_sort(vec.data(), vec.size());
}

// Sort signed integers
template <std::signed_integral T>
void radix_sort(T* data, size_t n) {
    Vector<T> buffer;
    detail::radix_sort_signed(data, n, buffer);
}

template <std::signed_integral T, typename Alloc>
void radix_sort(Vector<T, 0, Alloc>& vec) {
    radix_sort(vec.data(), vec.size());
}

template <std::signed_integral T, size_t N, typename Alloc>
void radix_sort(Vector<T, N, Alloc>& vec) {
    radix_sort(vec.data(), vec.size());
}

// Sort floats
template <std::floating_point T>
void radix_sort(T* data, size_t n) {
    Vector<T> buffer;
    detail::radix_sort_float(data, n, buffer);
}

template <std::floating_point T, typename Alloc>
void radix_sort(Vector<T, 0, Alloc>& vec) {
    radix_sort(vec.data(), vec.size());
}

// Sort key-value pairs (keys are unsigned integers)
template <detail::UnsignedIntegral Key, typename Value>
void radix_sort_pairs(Key* keys, Value* values, size_t n) {
    Vector<Key> key_buf;
    Vector<Value> val_buf;
    detail::radix_sort_pairs(keys, values, n, key_buf, val_buf);
}

template <detail::UnsignedIntegral Key, typename Value, typename Alloc>
void radix_sort_pairs(Vector<Key, 0, Alloc>& keys, Vector<Value, 0, Alloc>& values) {
    if (keys.size() != values.size()) return;
    radix_sort_pairs(keys.data(), values.data(), keys.size());
}

// Parallel introsort (wrapper around std::sort with execution policy)
template <std::random_access_iterator It, typename Comp = std::less<>>
void parallel_sort(It first, It last, Comp comp = {}) {
    std::sort(std::execution::par, first, last, comp);
}

template <typename Range, typename Comp = std::less<>>
void parallel_sort(Range&& r, Comp comp = {}) {
    parallel_sort(std::begin(r), std::end(r), comp);
}

// Stable parallel sort
template <std::random_access_iterator It, typename Comp = std::less<>>
void parallel_stable_sort(It first, It last, Comp comp = {}) {
    std::stable_sort(std::execution::par, first, last, comp);
}

// Small sort optimization (insertion sort for n <= 32)
template <std::random_access_iterator It, typename Comp = std::less<>>
void small_sort(It first, It last, Comp comp = {}) {
    size_t n = last - first;
    if (n <= 1) return;
    for (size_t i = 1; i < n; ++i) {
        auto key = std::move(first[i]);
        size_t j = i;
        while (j > 0 && comp(key, first[j - 1])) {
            first[j] = std::move(first[j - 1]);
            --j;
        }
        first[j] = std::move(key);
    }
}

// Hybrid sort: small -> introsort
template <std::random_access_iterator It, typename Comp = std::less<>>
void sort(It first, It last, Comp comp = {}) {
    size_t n = last - first;
    if (n <= 32) {
        small_sort(first, last, comp);
    } else {
        std::sort(first, last, comp);
    }
}

template <typename Range, typename Comp = std::less<>>
void sort(Range&& r, Comp comp = {}) {
    sort(std::begin(r), std::end(r), comp);
}

// Partial sort (top-k)
template <std::random_access_iterator It, typename Comp = std::less<>>
void partial_sort(It first, It middle, It last, Comp comp = {}) {
    std::partial_sort(first, middle, last, comp);
}

// nth_element
template <std::random_access_iterator It, typename Comp = std::less<>>
void nth_element(It first, It nth, It last, Comp comp = {}) {
    std::nth_element(first, nth, last, comp);
}

// Sort indices (argsort) - returns permutation vector
template <typename T, typename Comp = std::less<>>
Vector<size_t> argsort(const T* data, size_t n, Comp comp = {}) {
    Vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) { return comp(data[a], data[b]); });
    return idx;
}

template <typename T, typename Alloc, typename Comp = std::less<>>
Vector<size_t> argsort(const Vector<T, 0, Alloc>& vec, Comp comp = {}) {
    return argsort(vec.data(), vec.size(), comp);
}

// Apply permutation in-place
template <typename T>
void apply_permutation(T* data, const size_t* perm, size_t n) {
    Vector<char> visited(n, 0);
    for (size_t i = 0; i < n; ++i) {
        if (visited[i]) continue;
        size_t cur = i;
        T tmp = std::move(data[cur]);
        while (true) {
            size_t next = perm[cur];
            visited[cur] = 1;
            if (next == i) break;
            data[cur] = std::move(data[next]);
            cur = next;
        }
        data[cur] = std::move(tmp);
    }
}

template <typename T, typename Alloc>
void apply_permutation(Vector<T, 0, Alloc>& vec, const Vector<size_t>& perm) {
    if (vec.size() != perm.size()) return;
    apply_permutation(vec.data(), perm.data(), vec.size());
}

// Check if sorted
template <std::forward_iterator It, typename Comp = std::less<>>
bool is_sorted(It first, It last, Comp comp = {}) {
    return std::is_sorted(first, last, comp);
}

template <typename Range, typename Comp = std::less<>>
bool is_sorted(const Range& r, Comp comp = {}) {
    return is_sorted(std::begin(r), std::end(r), comp);
}

// Unique (remove consecutive duplicates)
template <std::forward_iterator It, typename BinaryPred = std::equal_to<>>
It unique(It first, It last, BinaryPred pred = {}) {
    return std::unique(first, last, pred);
}

template <typename Range, typename BinaryPred = std::equal_to<>>
auto unique(Range&& r, BinaryPred pred = {}) {
    return unique(std::begin(r), std::end(r), pred);
}

// Lower/upper bound
template <std::forward_iterator It, typename T, typename Comp = std::less<>>
It lower_bound(It first, It last, const T& value, Comp comp = {}) {
    return std::lower_bound(first, last, value, comp);
}

template <std::forward_iterator It, typename T, typename Comp = std::less<>>
It upper_bound(It first, It last, const T& value, Comp comp = {}) {
    return std::upper_bound(first, last, value, comp);
}

template <typename Range, typename T, typename Comp = std::less<>>
auto lower_bound(const Range& r, const T& value, Comp comp = {}) {
    return lower_bound(std::begin(r), std::end(r), value, comp);
}

template <typename Range, typename T, typename Comp = std::less<>>
auto upper_bound(const Range& r, const T& value, Comp comp = {}) {
    return upper_bound(std::begin(r), std::end(r), value, comp);
}

} // namespace bwa::core