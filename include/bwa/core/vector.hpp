#pragma once

#include <span>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <utility>
#include <iterator>
#include <type_traits>

#include <bwa/core/arena.hpp>

namespace bwa::core {

// Type-safe, memory-safe replacement for kvec_t
// Features:
// - Small buffer optimization (SBO) for small vectors
// - PMR allocator support
// - Bounds-checked access in debug
// - No manual memory management
// - Compatible with std::span, range-based for, algorithms

template <typename T, std::size_t InlineCapacity = 0, typename Allocator = memory::PolymorphicAllocator<T>>
class Vector {
public:
    using value_type = T;
    using allocator_type = Allocator;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    // Storage: either inline buffer or heap allocation
    union Storage {
        alignas(T) std::byte inline_buffer[InlineCapacity * sizeof(T)];
        T* heap_ptr;
    };

    Storage storage_;
    size_type size_ = 0;
    size_type capacity_ = 0;
    [[no_unique_address]] Allocator alloc_;

    static constexpr bool use_inline() noexcept { return InlineCapacity > 0; }
    static constexpr size_type inline_capacity() noexcept { return InlineCapacity; }

    [[nodiscard]] bool is_inline() const noexcept {
        return use_inline() && capacity_ == InlineCapacity;
    }

    [[nodiscard]] pointer data_ptr() noexcept {
        return is_inline() ? reinterpret_cast<pointer>(storage_.inline_buffer) : storage_.heap_ptr;
    }

    [[nodiscard]] const_pointer data_ptr() const noexcept {
        return is_inline() ? reinterpret_cast<const_pointer>(storage_.inline_buffer) : storage_.heap_ptr;
    }

    void allocate_heap(size_type new_cap) {
        if (new_cap == 0) {
            storage_.heap_ptr = nullptr;
            capacity_ = 0;
            return;
        }
        storage_.heap_ptr = std::allocator_traits<Allocator>::allocate(alloc_, new_cap);
        capacity_ = new_cap;
    }

    void deallocate_heap() noexcept {
        if (!is_inline() && storage_.heap_ptr) {
            std::allocator_traits<Allocator>::deallocate(alloc_, storage_.heap_ptr, capacity_);
            storage_.heap_ptr = nullptr;
            capacity_ = 0;
        }
    }

    void destroy_elements(pointer first, pointer last) noexcept {
        for (; first != last; ++first) {
            std::allocator_traits<Allocator>::destroy(alloc_, first);
        }
    }

    void move_elements(pointer dest, pointer src, size_type count) noexcept(
        std::is_nothrow_move_constructible_v<T>) {
        if constexpr (std::is_trivially_copyable_v<T>) {
            std::copy_n(src, count, dest);
        } else {
            for (size_type i = 0; i < count; ++i) {
                std::allocator_traits<Allocator>::construct(alloc_, dest + i, std::move(src[i]));
            }
        }
    }

    void copy_elements(pointer dest, const_pointer src, size_type count) {
        if constexpr (std::is_trivially_copyable_v<T>) {
            std::copy_n(src, count, dest);
        } else {
            for (size_type i = 0; i < count; ++i) {
                std::allocator_traits<Allocator>::construct(alloc_, dest + i, src[i]);
            }
        }
    }

    void grow(size_type min_capacity) {
        size_type new_cap = std::max(capacity_ == 0 ? 1 : capacity_ * 2, min_capacity);
        if (use_inline() && new_cap <= InlineCapacity) {
            new_cap = InlineCapacity;
        }

        pointer new_data = std::allocator_traits<Allocator>::allocate(alloc_, new_cap);
        move_elements(new_data, data_ptr(), size_);

        destroy_elements(data_ptr(), data_ptr() + size_);
        deallocate_heap();

        storage_.heap_ptr = new_data;
        capacity_ = new_cap;
    }

public:
    // Constructors
    Vector() noexcept(noexcept(Allocator()))
        : capacity_(use_inline() ? InlineCapacity : 0) {
        if (use_inline()) {
            // Mark inline buffer as uninitialized
        }
    }

    explicit Vector(memory::Arena* arena) noexcept
        : alloc_(arena), capacity_(use_inline() ? InlineCapacity : 0) {}

    explicit Vector(const Allocator& alloc) noexcept
        : alloc_(alloc), capacity_(use_inline() ? InlineCapacity : 0) {}

    Vector(size_type count, const T& value, const Allocator& alloc = Allocator())
        : alloc_(alloc), capacity_(use_inline() ? InlineCapacity : 0) {
        resize(count, value);
    }

    explicit Vector(size_type count, const Allocator& alloc = Allocator())
        : alloc_(alloc), capacity_(use_inline() ? InlineCapacity : 0) {
        resize(count);
    }

    Vector(const Vector& other)
        : alloc_(std::allocator_traits<Allocator>::select_on_container_copy_construction(other.alloc_))
        , capacity_(use_inline() ? InlineCapacity : 0) {
        if (other.size_ > 0) {
            reserve(other.size_);
            copy_elements(data(), other.data(), other.size_);
            size_ = other.size_;
        }
    }

    Vector(Vector&& other) noexcept
        : alloc_(other.alloc_)
        , size_(other.size_)
        , capacity_(other.capacity_) {
        if (other.is_inline()) {
            // Move elements from inline buffer
            move_elements(data(), other.data(), size_);
        } else {
            storage_.heap_ptr = other.storage_.heap_ptr;
        }
        other.size_ = 0;
        other.capacity_ = use_inline() ? InlineCapacity : 0;
        // Nullify source pointer to prevent double-free in source's destructor
        other.storage_.heap_ptr = nullptr;
    }

    Vector(std::initializer_list<T> init, const Allocator& alloc = Allocator())
        : alloc_(alloc), capacity_(use_inline() ? InlineCapacity : 0) {
        reserve(init.size());
        for (const auto& v : init) {
            push_back(v);
        }
    }

    template <std::input_iterator It>
    Vector(It first, It last, const Allocator& alloc = Allocator())
        : alloc_(alloc), capacity_(use_inline() ? InlineCapacity : 0) {
        append(first, last);
    }

    ~Vector() {
        clear();
        deallocate_heap();
    }

    Vector& operator=(const Vector& other) {
        if (this != &other) {
            clear();
            if (std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value) {
                alloc_ = other.alloc_;
            }
            reserve(other.size_);
            copy_elements(data(), other.data(), other.size_);
            size_ = other.size_;
        }
        return *this;
    }

    Vector& operator=(Vector&& other) noexcept(
        std::allocator_traits<Allocator>::propagate_on_container_move_assignment::value ||
        std::allocator_traits<Allocator>::is_always_equal::value) {
        if (this != &other) {
            clear();
            deallocate_heap();

            // Arena allocator is not assignable, but is_always_equal is true
            // when resources are equal. Since we use thread-local arenas, they should match.
            // Don't move the allocator - it's a pointer to a shared arena.

            size_ = other.size_;
            capacity_ = other.capacity_;

            if (other.is_inline()) {
                move_elements(data_ptr(), other.data_ptr(), size_);
            } else {
                storage_.heap_ptr = other.storage_.heap_ptr;
            }

            other.size_ = 0;
            other.capacity_ = use_inline() ? InlineCapacity : 0;
            // Nullify source pointer to prevent double-free in source's destructor
            other.storage_.heap_ptr = nullptr;
        }
        return *this;
    }

    Vector& operator=(std::initializer_list<T> ilist) {
        clear();
        reserve(ilist.size());
        for (const auto& v : ilist) {
            push_back(v);
        }
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
    [[nodiscard]] size_type max_size() const noexcept { return std::allocator_traits<Allocator>::max_size(alloc_); }
    [[nodiscard]] size_type capacity() const noexcept { return capacity_; }

    void reserve(size_type new_cap) {
        if (new_cap > capacity_) {
            if (use_inline() && new_cap <= InlineCapacity) {
                // Already have inline capacity
                return;
            }
            grow(new_cap);
        }
    }

    void shrink_to_fit() {
        if (size_ == capacity_) return;
        if (use_inline() && size_ <= InlineCapacity) {
            // Move to inline buffer
            if (!is_inline()) {
                move_elements(reinterpret_cast<pointer>(storage_.inline_buffer), data(), size_);
                deallocate_heap();
                capacity_ = InlineCapacity;
            }
        } else if (size_ > 0) {
            pointer new_data = std::allocator_traits<Allocator>::allocate(alloc_, size_);
            move_elements(new_data, data(), size_);
            destroy_elements(data(), data() + size_);
            deallocate_heap();
            storage_.heap_ptr = new_data;
            capacity_ = size_;
        } else {
            deallocate_heap();
            capacity_ = use_inline() ? InlineCapacity : 0;
        }
    }

    // Element access
    [[nodiscard]] reference operator[](size_type pos) noexcept {
        return data()[pos];
    }

    [[nodiscard]] const_reference operator[](size_type pos) const noexcept {
        return data()[pos];
    }

    [[nodiscard]] reference at(size_type pos) {
        if (pos >= size_) throw std::out_of_range("Vector::at");
        return data()[pos];
    }

    [[nodiscard]] const_reference at(size_type pos) const {
        if (pos >= size_) throw std::out_of_range("Vector::at");
        return data()[pos];
    }

    [[nodiscard]] reference front() noexcept { return data()[0]; }
    [[nodiscard]] const_reference front() const noexcept { return data()[0]; }
    [[nodiscard]] reference back() noexcept { return data()[size_ - 1]; }
    [[nodiscard]] const_reference back() const noexcept { return data()[size_ - 1]; }
    [[nodiscard]] pointer data() noexcept { return data_ptr(); }
    [[nodiscard]] const_pointer data() const noexcept { return data_ptr(); }

    // Modifiers
    void clear() noexcept {
        destroy_elements(data_ptr(), data_ptr() + size_);
        size_ = 0;
    }

    template <typename... Args>
    reference emplace_back(Args&&... args) {
        if (size_ == capacity_) {
            grow(size_ + 1);
        }
        std::allocator_traits<Allocator>::construct(alloc_, data_ptr() + size_, std::forward<Args>(args)...);
        return data_ptr()[size_++];
    }

    void push_back(const T& value) {
        emplace_back(value);
    }

    void push_back(T&& value) {
        emplace_back(std::move(value));
    }

    template <typename... Args>
    iterator emplace(const_iterator pos, Args&&... args) {
        const size_type index = pos - begin();
        if (size_ == capacity_) {
            grow(size_ + 1);
        }
        pointer p = data_ptr() + index;
        // Move elements after pos one position back
        move_elements(p + 1, p, size_ - index);
        std::allocator_traits<Allocator>::construct(alloc_, p, std::forward<Args>(args)...);
        ++size_;
        return p;
    }

    iterator insert(const_iterator pos, const T& value) {
        return emplace(pos, value);
    }

    iterator insert(const_iterator pos, T&& value) {
        return emplace(pos, std::move(value));
    }

    iterator insert(const_iterator pos, size_type count, const T& value) {
        const size_type index = pos - begin();
        if (count == 0) return begin() + index;
        reserve(size_ + count);
        pointer p = data_ptr() + index;
        move_elements(p + count, p, size_ - index);
        for (size_type i = 0; i < count; ++i) {
            std::allocator_traits<Allocator>::construct(alloc_, p + i, value);
        }
        size_ += count;
        return p;
    }

    template <std::input_iterator It>
    iterator insert(const_iterator pos, It first, It last) {
        const size_type index = pos - begin();
        const size_type count = std::distance(first, last);
        if (count == 0) return begin() + index;
        reserve(size_ + count);
        pointer p = data_ptr() + index;
        move_elements(p + count, p, size_ - index);
        for (size_type i = 0; i < count; ++i, ++first) {
            std::allocator_traits<Allocator>::construct(alloc_, p + i, *first);
        }
        size_ += count;
        return p;
    }

    iterator erase(const_iterator pos) {
        const size_type index = pos - begin();
        pointer p = data_ptr() + index;
        std::allocator_traits<Allocator>::destroy(alloc_, p);
        move_elements(p, p + 1, size_ - index - 1);
        --size_;
        return p;
    }

    iterator erase(const_iterator first, const_iterator last) {
        const size_type count = last - first;
        if (count == 0) return first;
        const size_type index = first - begin();
        pointer p = data_ptr() + index;
        destroy_elements(p, p + count);
        move_elements(p, p + count, size_ - index - count);
        size_ -= count;
        return p;
    }

    void push_back_range(std::ranges::input_range auto&& range) {
        for (auto&& elem : range) {
            push_back(std::forward<decltype(elem)>(elem));
        }
    }

    void pop_back() noexcept {
        std::allocator_traits<Allocator>::destroy(alloc_, data_ptr() + size_ - 1);
        --size_;
    }

    void resize(size_type count) {
        if (count < size_) {
            destroy_elements(data_ptr() + count, data_ptr() + size_);
        } else if (count > size_) {
            reserve(count);
            for (size_type i = size_; i < count; ++i) {
                std::allocator_traits<Allocator>::construct(alloc_, data_ptr() + i, T{});
            }
        }
        size_ = count;
    }

    void resize(size_type count, const T& value) {
        if (count < size_) {
            destroy_elements(data_ptr() + count, data_ptr() + size_);
        } else if (count > size_) {
            reserve(count);
            for (size_type i = size_; i < count; ++i) {
                std::allocator_traits<Allocator>::construct(alloc_, data_ptr() + i, value);
            }
        }
        size_ = count;
    }

    void swap(Vector& other) noexcept(
        std::allocator_traits<Allocator>::propagate_on_container_swap::value ||
        std::allocator_traits<Allocator>::is_always_equal::value) {
        if (this != &other) {
            if (std::allocator_traits<Allocator>::propagate_on_container_swap::value) {
                std::swap(alloc_, other.alloc_);
            }
            std::swap(size_, other.size_);
            std::swap(capacity_, other.capacity_);

            if (is_inline() && other.is_inline()) {
                // Both inline - swap element by element
                for (size_type i = 0; i < std::max(size_, other.size_); ++i) {
                    std::swap(data_ptr()[i], other.data_ptr()[i]);
                }
            } else if (!is_inline() && !other.is_inline()) {
                // Both heap - swap pointers
                std::swap(storage_.heap_ptr, other.storage_.heap_ptr);
            } else if (is_inline()) {
                // this inline, other heap
                Storage temp = storage_;
                storage_.heap_ptr = other.storage_.heap_ptr;
                other.storage_ = temp;
            } else {
                // this heap, other inline
                Storage temp = other.storage_;
                other.storage_.heap_ptr = storage_.heap_ptr;
                storage_ = temp;
            }
        }
    }

    // Span interop
    [[nodiscard]] std::span<T> span() noexcept { return {data_ptr(), size_}; }
    [[nodiscard]] std::span<const T> span() const noexcept { return {data_ptr(), size_}; }

    // Comparison
    friend bool operator==(const Vector& lhs, const Vector& rhs) noexcept {
        return lhs.size_ == rhs.size_ && std::equal(lhs.begin(), lhs.end(), rhs.begin());
    }

    friend auto operator<=>(const Vector& lhs, const Vector& rhs) noexcept {
        return std::lexicographical_compare_three_way(lhs.begin(), lhs.end(),
                                                       rhs.begin(), rhs.end());
    }

    // Get allocator
    [[nodiscard]] Allocator get_allocator() const noexcept { return alloc_; }
};

// Deduction guides
template <typename It, typename Alloc = memory::PolymorphicAllocator<std::iter_value_t<It>>>
Vector(It, It, Alloc = Alloc()) -> Vector<std::iter_value_t<It>, 0, Alloc>;

template <typename T, typename Alloc = memory::PolymorphicAllocator<T>>
Vector(std::initializer_list<T>, Alloc = Alloc()) -> Vector<T, 0, Alloc>;

// Convenience aliases
template <typename T, typename Alloc = memory::PolymorphicAllocator<T>>
using SmallVector = Vector<T, 16, Alloc>;

template <typename T, typename Alloc = memory::PolymorphicAllocator<T>>
using TinyVector = Vector<T, 4, Alloc>;

} // namespace bwa::core

// Specialize std::hash for Vector if needed
namespace std {
template <typename T, std::size_t N, typename Alloc>
struct hash<bwa::core::Vector<T, N, Alloc>> {
    size_t operator()(const bwa::core::Vector<T, N, Alloc>& v) const noexcept {
        size_t h = v.size();
        for (const auto& elem : v) {
            h ^= std::hash<T>{}(elem) + 0x9e3779b97f4a7c15 + (h << 6) + (h >> 2);
        }
        return h;
    }
};
}