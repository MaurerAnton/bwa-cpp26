#pragma once

#include <cstddef>
#include <cstdint>
#include <new>
#include <limits>
#include <algorithm>
#include <span>
#include <memory>
#include <vector>
#include <string>

namespace bwa::memory {

// Simple memory resource using plain new/delete (no std::pmr dependency)
class DefaultMemoryResource {
public:
    void* allocate(std::size_t bytes, std::size_t alignment) {
        // Align the allocation
        void* ptr = ::operator new(bytes, std::align_val_t(alignment));
        return ptr;
    }

    void deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
        ::operator delete(ptr, std::align_val_t(alignment));
    }

    bool is_equal(const DefaultMemoryResource&) const noexcept { return true; }
};

// Arena allocator using our own memory resource
// Provides fast bump-pointer allocation with automatic cleanup
class Arena {
    DefaultMemoryResource upstream_;
    std::byte* buffer_ = nullptr;
    std::size_t capacity_ = 0;
    std::size_t offset_ = 0;

public:
    explicit Arena(std::size_t initial_capacity = 64 * 1024)
        : buffer_(static_cast<std::byte*>(upstream_.allocate(initial_capacity, alignof(std::max_align_t))))
        , capacity_(initial_capacity)
        , offset_(0) {}

    ~Arena() {
        if (buffer_) {
            upstream_.deallocate(buffer_, capacity_, alignof(std::max_align_t));
        }
    }

    Arena(const Arena&) = delete;
    Arena& operator=(const Arena&) = delete;

    Arena(Arena&& other) noexcept
        : upstream_()
        , buffer_(other.buffer_)
        , capacity_(other.capacity_)
        , offset_(other.offset_) {
        other.buffer_ = nullptr;
        other.capacity_ = 0;
        other.offset_ = 0;
    }

    Arena& operator=(Arena&& other) noexcept {
        if (this != &other) {
            if (buffer_) {
                upstream_.deallocate(buffer_, capacity_, alignof(std::max_align_t));
            }
            buffer_ = other.buffer_;
            capacity_ = other.capacity_;
            offset_ = other.offset_;
            other.buffer_ = nullptr;
            other.capacity_ = 0;
            other.offset_ = 0;
        }
        return *this;
    }

    [[nodiscard]] std::size_t used() const noexcept { return offset_; }
    [[nodiscard]] std::size_t capacity() const noexcept { return capacity_; }
    [[nodiscard]] std::size_t remaining() const noexcept { return capacity_ - offset_; }

    void reset() noexcept { offset_ = 0; }

    void shrink_to_fit() {
        if (offset_ == 0) return;
        std::byte* new_buffer = static_cast<std::byte*>(upstream_.allocate(offset_, alignof(std::max_align_t)));
        std::copy_n(buffer_, offset_, new_buffer);
        upstream_.deallocate(buffer_, capacity_, alignof(std::max_align_t));
        buffer_ = new_buffer;
        capacity_ = offset_;
    }

    void* allocate(std::size_t bytes, std::size_t alignment) {
        const std::size_t aligned_offset = (offset_ + alignment - 1) & ~(alignment - 1);
        if (aligned_offset + bytes <= capacity_) {
            void* ptr = buffer_ + aligned_offset;
            offset_ = aligned_offset + bytes;
            return ptr;
        }
        // Fallback to upstream for large allocations or when arena is full
        return upstream_.allocate(bytes, alignment);
    }

    void deallocate(void* ptr, std::size_t bytes, std::size_t alignment) {
        // Arena doesn't support individual deallocation - use reset() or let arena die
        // Only deallocate if it was allocated from upstream (not from our buffer)
        std::byte* byte_ptr = static_cast<std::byte*>(ptr);
        if (byte_ptr < buffer_ || byte_ptr >= buffer_ + capacity_) {
            upstream_.deallocate(ptr, bytes, alignment);
        }
    }

    // Thread-local arena for per-thread allocations
    static Arena* get_tls_arena_ptr(std::size_t initial_capacity = 64 * 1024) noexcept {
        thread_local Arena arena{initial_capacity};
        return &arena;
    }

    static Arena& get_tls_arena(std::size_t initial_capacity = 64 * 1024) noexcept {
        return *get_tls_arena_ptr(initial_capacity);
    }

    static void reset_tls_arena() noexcept {
        if (auto* arena = get_tls_arena_ptr(0)) {
            arena->reset();
        }
    }
};

// Free functions for backward compatibility
inline Arena* get_tls_arena_ptr(std::size_t initial_capacity = 64 * 1024) noexcept {
    return Arena::get_tls_arena_ptr(initial_capacity);
}

inline Arena& get_tls_arena(std::size_t initial_capacity = 64 * 1024) noexcept {
    return Arena::get_tls_arena(initial_capacity);
}

inline void reset_tls_arena() noexcept {
    Arena::reset_tls_arena();
}

// Simple assignable allocator that wraps an Arena pointer
template <typename T>
class ArenaAllocator {
    Arena* arena_ = nullptr;
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using void_pointer = void*;
    using const_void_pointer = const void*;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template <typename U> struct rebind { using other = ArenaAllocator<U>; };

    ArenaAllocator() noexcept : arena_(Arena::get_tls_arena_ptr()) {}
    explicit ArenaAllocator(Arena* arena) noexcept : arena_(arena ? arena : Arena::get_tls_arena_ptr()) {}
    ArenaAllocator(const ArenaAllocator&) noexcept = default;
    ArenaAllocator& operator=(const ArenaAllocator&) noexcept = default;

    template <typename U>
    ArenaAllocator(const ArenaAllocator<U>& other) noexcept : arena_(other.arena_) {}

    [[nodiscard]] T* allocate(size_type n) {
        return static_cast<T*>(arena_->allocate(n * sizeof(T), alignof(T)));
    }

    void deallocate(T* p, size_type n) noexcept {
        arena_->deallocate(p, n * sizeof(T), alignof(T));
    }

    template <typename U, typename... Args>
    void construct(U* p, Args&&... args) {
        ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...);
    }

    template <typename U>
    void default_construct(U* p) {
        ::new (static_cast<void*>(p)) U();
    }

    template <typename U>
    void destroy(U* p) {
        p->~U();
    }

    [[nodiscard]] size_type max_size() const noexcept {
        return std::numeric_limits<size_type>::max() / sizeof(T);
    }

    [[nodiscard]] Arena* resource() const noexcept { return arena_; }

    friend bool operator==(const ArenaAllocator& a, const ArenaAllocator& b) noexcept {
        return a.arena_ == b.arena_;
    }
    friend bool operator!=(const ArenaAllocator& a, const ArenaAllocator& b) noexcept {
        return a.arena_ != b.arena_;
    }
};

// Type aliases using ArenaAllocator
template <typename T>
using PolymorphicAllocator = ArenaAllocator<T>;

template <typename T>
using PmrVector = std::vector<T, ArenaAllocator<T>>;

// std::pmr::deque not available until C++23 - use vector instead
template <typename T>
using PmrDeque = std::vector<T, ArenaAllocator<T>>;

using PmrString = std::basic_string<char, std::char_traits<char>, ArenaAllocator<char>>;
using PmrWString = std::basic_string<wchar_t, std::char_traits<wchar_t>, ArenaAllocator<wchar_t>>;

} // namespace bwa::memory