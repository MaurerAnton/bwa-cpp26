#pragma once

#include <bwa/core/arena.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <functional>
#include <memory>
#include <utility>
#include <iterator>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <bit>
#include <concepts>

namespace bwa::core {

// Hash map replacement for khash.h
// Features:
// - Robin hood hashing with backward shift deletion (fast, cache-friendly)
// - PMR allocator support
// - Heterogeneous lookup (C++20)
// - No exceptions in hot path
// - Compatible with khash_t API for gradual translation

namespace detail {

// Robin hood hash table implementation
template <typename Key, typename Value, typename Hash, typename KeyEqual, typename Allocator>
class RobinHoodHashMap {
public:
    using key_type = Key;
    using mapped_type = Value;
    using value_type = std::pair<const Key, Value>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using hasher = Hash;
    using key_equal = KeyEqual;
    using allocator_type = Allocator;
    using reference = value_type&;
    using const_reference = const value_type&;
    using pointer = value_type*;
    using const_pointer = const value_type*;

private:
    // Entry states
    enum class State : uint8_t { Empty, Occupied, Deleted };

    struct Entry {
        Key key;
        Value value;
        State state = State::Empty;
        uint32_t probe_dist = 0; // Distance from ideal position

        Entry() = default;
        Entry(const Key& k, const Value& v, uint32_t pd)
            : key(k), value(v), state(State::Occupied), probe_dist(pd) {}
        Entry(Key&& k, Value&& v, uint32_t pd)
            : key(std::move(k)), value(std::move(v)), state(State::Occupied), probe_dist(pd) {}
    };

    using EntryAllocator = typename std::allocator_traits<Allocator>::template rebind_alloc<Entry>;
    using EntryAllocTraits = std::allocator_traits<EntryAllocator>;

    Entry* entries_ = nullptr;
    size_type capacity_ = 0;
    size_type size_ = 0;
    size_type deleted_ = 0;
    [[no_unique_address]] Hash hasher_;
    [[no_unique_address]] KeyEqual key_equal_;
    [[no_unique_address]] EntryAllocator alloc_;

    static constexpr float MAX_LOAD_FACTOR = 0.9f;
    static constexpr size_type INITIAL_CAPACITY = 16;

    [[nodiscard]] size_type hash_to_bucket(size_type hash) const noexcept {
        return hash & (capacity_ - 1);
    }

    [[nodiscard]] size_type next_bucket(size_type bucket) const noexcept {
        return (bucket + 1) & (capacity_ - 1);
    }

    void destroy_entries() noexcept {
        if (entries_) {
            for (size_type i = 0; i < capacity_; ++i) {
                if (entries_[i].state == State::Occupied) {
                    EntryAllocTraits::destroy(alloc_, &entries_[i]);
                }
            }
            EntryAllocTraits::deallocate(alloc_, entries_, capacity_);
            entries_ = nullptr;
            capacity_ = 0;
        }
    }

    void allocate_entries(size_type new_cap) {
        if (new_cap == 0) new_cap = INITIAL_CAPACITY;
        // Round up to power of 2
        new_cap = std::bit_ceil(new_cap);
        entries_ = EntryAllocTraits::allocate(alloc_, new_cap);
        for (size_type i = 0; i < new_cap; ++i) {
            EntryAllocTraits::construct(alloc_, &entries_[i]);
        }
        capacity_ = new_cap;
    }

    void rehash(size_type new_cap) {
        Entry* old_entries = entries_;
        size_type old_cap = capacity_;
        size_type old_size = size_;

        allocate_entries(new_cap);
        size_ = 0;
        deleted_ = 0;

        for (size_type i = 0; i < old_cap; ++i) {
            if (old_entries[i].state == State::Occupied) {
                insert_impl(std::move(old_entries[i].key), std::move(old_entries[i].value));
                EntryAllocTraits::destroy(alloc_, &old_entries[i]);
            }
        }
        EntryAllocTraits::deallocate(alloc_, old_entries, old_cap);
    }

    template <typename K, typename V>
    std::pair<Entry*, bool> insert_impl(K&& k, V&& v) {
        if (size_ + deleted_ >= capacity_ * MAX_LOAD_FACTOR) {
            rehash(capacity_ * 2);
        }

        size_type hash = hasher_(k);
        size_type bucket = hash_to_bucket(hash);
        uint32_t probe_dist = 0;

        while (true) {
            Entry& entry = entries_[bucket];
            if (entry.state == State::Empty) {
                // Found empty slot - insert here
                EntryAllocTraits::construct(alloc_, &entry, std::forward<K>(k), std::forward<V>(v), probe_dist);
                entry.state = State::Occupied;
                ++size_;
                return {&entry, true};
            }

            if (entry.state == State::Occupied) {
                size_type existing_hash = hasher_(entry.key);
                size_type ideal_bucket = hash_to_bucket(existing_hash);
                uint32_t existing_probe_dist = (bucket >= ideal_bucket) ?
                    bucket - ideal_bucket : bucket + capacity_ - ideal_bucket;

                if (key_equal_(entry.key, k)) {
                    // Key exists - update value
                    entry.value = std::forward<V>(v);
                    return {&entry, false};
                }

                // Robin hood: if we've probed further than existing entry, swap
                if (probe_dist > existing_probe_dist) {
                    std::swap(entry.key, k);
                    std::swap(entry.value, v);
                    std::swap(entry.probe_dist, probe_dist);
                    // Continue with displaced entry
                }
            }
            // State::Deleted - continue probing

            bucket = next_bucket(bucket);
            ++probe_dist;
        }
    }

    template <typename K>
    pointer find_impl(const K& k) const noexcept {
        if (capacity_ == 0) return nullptr;

        size_type hash = hasher_(k);
        size_type bucket = hash_to_bucket(hash);
        uint32_t probe_dist = 0;

        while (true) {
            const Entry& entry = entries_[bucket];
            if (entry.state == State::Empty) {
                return nullptr; // Not found
            }
            if (entry.state == State::Occupied && key_equal_(entry.key, k)) {
                return reinterpret_cast<pointer>(const_cast<Entry*>(&entry));
            }

            // Early exit: if we've probed further than any occupied entry could be
            if (entry.state == State::Occupied) {
                size_type existing_hash = hasher_(entry.key);
                size_type ideal_bucket = hash_to_bucket(existing_hash);
                uint32_t existing_probe_dist = (bucket >= ideal_bucket) ?
                    bucket - ideal_bucket : bucket + capacity_ - ideal_bucket;
                if (probe_dist > existing_probe_dist) {
                    return nullptr;
                }
            }

            bucket = next_bucket(bucket);
            ++probe_dist;
            if (probe_dist >= capacity_) return nullptr; // Safety
        }
    }

    template <typename K>
    bool erase_impl(const K& k) {
        if (capacity_ == 0) return false;

        size_type hash = hasher_(k);
        size_type bucket = hash_to_bucket(hash);
        uint32_t probe_dist = 0;

        while (true) {
            Entry& entry = entries_[bucket];
            if (entry.state == State::Empty) {
                return false;
            }
            if (entry.state == State::Occupied && key_equal_(entry.key, k)) {
                EntryAllocTraits::destroy(alloc_, &entry);
                entry.state = State::Deleted;
                --size_;
                ++deleted_;

                // Backward shift deletion to maintain Robin hood property
                size_type next = next_bucket(bucket);
                while (entries_[next].state == State::Occupied) {
                    size_type next_hash = hasher_(entries_[next].key);
                    size_type next_ideal = hash_to_bucket(next_hash);
                    // Check if next entry can be moved back
                    if (next_ideal == bucket || next_ideal == next_bucket(bucket)) {
                        // Move it back
                        EntryAllocTraits::construct(alloc_, &entry,
                            std::move(entries_[next].key), std::move(entries_[next].value),
                            entries_[next].probe_dist - 1);
                        entry.state = State::Occupied;
                        EntryAllocTraits::destroy(alloc_, &entries_[next]);
                        entries_[next].state = State::Deleted;
                        bucket = next;
                        next = next_bucket(bucket);
                    } else {
                        break;
                    }
                }
                return true;
            }
            bucket = next_bucket(bucket);
            ++probe_dist;
            if (probe_dist >= capacity_) return false;
        }
    }

public:
    // Iterator
    class iterator {
        Entry* ptr_;
        Entry* end_;
        void advance() {
            while (ptr_ != end_ && ptr_->state != State::Occupied) ++ptr_;
        }
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::pair<const Key, Value>;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = value_type&;

        iterator() = default;
        iterator(Entry* ptr, Entry* end) : ptr_(ptr), end_(end) { advance(); }

        reference operator*() const noexcept { return *reinterpret_cast<pointer>(ptr_); }
        pointer operator->() const noexcept { return reinterpret_cast<pointer>(ptr_); }

        iterator& operator++() { ++ptr_; advance(); return *this; }
        iterator operator++(int) { auto tmp = *this; ++*this; return tmp; }

        friend bool operator==(const iterator& a, const iterator& b) noexcept { return a.ptr_ == b.ptr_; }
        friend bool operator!=(const iterator& a, const iterator& b) noexcept { return a.ptr_ != b.ptr_; }
    };

    class const_iterator {
        const Entry* ptr_;
        const Entry* end_;
        void advance() {
            while (ptr_ != end_ && ptr_->state != State::Occupied) ++ptr_;
        }
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = std::pair<const Key, Value>;
        using difference_type = std::ptrdiff_t;
        using pointer = const value_type*;
        using reference = const value_type&;

        const_iterator() = default;
        const_iterator(const Entry* ptr, const Entry* end) : ptr_(ptr), end_(end) { advance(); }
        const_iterator(const iterator& other) : ptr_(other.ptr_), end_(other.end_) {}

        reference operator*() const noexcept { return *reinterpret_cast<pointer>(ptr_); }
        pointer operator->() const noexcept { return reinterpret_cast<pointer>(ptr_); }

        const_iterator& operator++() { ++ptr_; advance(); return *this; }
        const_iterator operator++(int) { auto tmp = *this; ++*this; return tmp; }

        friend bool operator==(const const_iterator& a, const const_iterator& b) noexcept { return a.ptr_ == b.ptr_; }
        friend bool operator!=(const const_iterator& a, const const_iterator& b) noexcept { return a.ptr_ != b.ptr_; }
    };

    // Constructors
    RobinHoodHashMap() : alloc_(memory::get_tls_arena_ptr()) {
        allocate_entries(INITIAL_CAPACITY);
    }

    explicit RobinHoodHashMap(const Allocator& alloc) : alloc_(alloc) {
        allocate_entries(INITIAL_CAPACITY);
    }

    RobinHoodHashMap(size_type n, const Hash& hash = Hash{}, const KeyEqual& eq = KeyEqual{},
                     const Allocator& alloc = Allocator{memory::get_tls_arena_ptr()})
        : hasher_(hash), key_equal_(eq), alloc_(alloc) {
        allocate_entries(std::max(n * 2, INITIAL_CAPACITY));
    }

    ~RobinHoodHashMap() { destroy_entries(); }

    RobinHoodHashMap(const RobinHoodHashMap&) = delete;
    RobinHoodHashMap& operator=(const RobinHoodHashMap&) = delete;

    RobinHoodHashMap(RobinHoodHashMap&& other) noexcept
        : entries_(other.entries_), capacity_(other.capacity_), size_(other.size_),
          deleted_(other.deleted_), hasher_(std::move(other.hasher_)),
          key_equal_(std::move(other.key_equal_)), alloc_(std::move(other.alloc_)) {
        other.entries_ = nullptr;
        other.capacity_ = 0;
        other.size_ = 0;
        other.deleted_ = 0;
    }

    RobinHoodHashMap& operator=(RobinHoodHashMap&& other) noexcept {
        if (this != &other) {
            destroy_entries();
            entries_ = other.entries_;
            capacity_ = other.capacity_;
            size_ = other.size_;
            deleted_ = other.deleted_;
            hasher_ = std::move(other.hasher_);
            key_equal_ = std::move(other.key_equal_);
            alloc_ = std::move(other.alloc_);
            other.entries_ = nullptr;
            other.capacity_ = 0;
            other.size_ = 0;
            other.deleted_ = 0;
        }
        return *this;
    }

    // Iterators
    [[nodiscard]] iterator begin() noexcept { return {entries_, entries_ + capacity_}; }
    [[nodiscard]] const_iterator begin() const noexcept { return {entries_, entries_ + capacity_}; }
    [[nodiscard]] const_iterator cbegin() const noexcept { return begin(); }
    [[nodiscard]] iterator end() noexcept { return {entries_ + capacity_, entries_ + capacity_}; }
    [[nodiscard]] const_iterator end() const noexcept { return {entries_ + capacity_, entries_ + capacity_}; }
    [[nodiscard]] const_iterator cend() const noexcept { return end(); }

    // Capacity
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
    [[nodiscard]] size_type size() const noexcept { return size_; }
    [[nodiscard]] size_type max_size() const noexcept { return EntryAllocTraits::max_size(alloc_); }
    [[nodiscard]] size_type capacity() const noexcept { return capacity_; }
    [[nodiscard]] size_type bucket_count() const noexcept { return capacity_; }
    [[nodiscard]] float load_factor() const noexcept { return capacity_ ? float(size_) / capacity_ : 0.0f; }
    [[nodiscard]] float max_load_factor() const noexcept { return MAX_LOAD_FACTOR; }

    void reserve(size_type n) {
        if (n > size_ / MAX_LOAD_FACTOR) {
            rehash(std::max(n, capacity_));
        }
    }

    // Element access
    template <typename K>
    Value& operator[](K&& k) {
        auto [entry, inserted] = insert_impl(std::forward<K>(k), Value{});
        return entry->second;
    }

    template <typename K>
    Value& at(const K& k) {
        auto* entry = find_impl(k);
        if (!entry) throw std::out_of_range("RobinHoodHashMap::at");
        return entry->second;
    }

    template <typename K>
    const Value& at(const K& k) const {
        auto* entry = find_impl(k);
        if (!entry) throw std::out_of_range("RobinHoodHashMap::at");
        return entry->second;
    }

    // Lookup
    template <typename K>
    [[nodiscard]] iterator find(const K& k) noexcept {
        auto* entry = find_impl(k);
        return entry ? iterator(reinterpret_cast<Entry*>(entry), entries_ + capacity_) : end();
    }

    template <typename K>
    [[nodiscard]] const_iterator find(const K& k) const noexcept {
        auto* entry = find_impl(k);
        return entry ? const_iterator(reinterpret_cast<const Entry*>(entry), entries_ + capacity_) : end();
    }

    template <typename K>
    [[nodiscard]] bool contains(const K& k) const noexcept {
        return find_impl(k) != nullptr;
    }

    template <typename K>
    [[nodiscard]] size_type count(const K& k) const noexcept {
        return contains(k) ? 1 : 0;
    }

    template <typename K>
    [[nodiscard]] std::pair<iterator, iterator> equal_range(const K& k) noexcept {
        auto it = find(k);
        return {it, it == end() ? end() : std::next(it)};
    }

    template <typename K>
    [[nodiscard]] std::pair<const_iterator, const_iterator> equal_range(const K& k) const noexcept {
        auto it = find(k);
        return {it, it == end() ? end() : std::next(it)};
    }

    // Insert
    template <typename K, typename V>
    std::pair<iterator, bool> insert(K&& k, V&& v) {
        auto [entry, inserted] = insert_impl(std::forward<K>(k), std::forward<V>(v));
        return {iterator(entry, entries_ + capacity_), inserted};
    }

    template <typename K, typename V>
    std::pair<iterator, bool> emplace(K&& k, V&& v) {
        return insert(std::forward<K>(k), std::forward<V>(v));
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        return insert(std::forward<Args>(args)...);
    }

    // Erase
    template <typename K>
    size_type erase(const K& k) {
        return erase_impl(k) ? 1 : 0;
    }

    iterator erase(iterator pos) {
        auto next = std::next(pos);
        erase_impl(pos->first);
        return next;
    }

    iterator erase(const_iterator pos) {
        return erase(iterator(const_cast<Entry*>(pos.ptr_), pos.end_));
    }

    iterator erase(const_iterator first, const_iterator last) {
        while (first != last) {
            first = erase(first);
        }
        return iterator(const_cast<Entry*>(last.ptr_), last.end_);
    }

    void clear() noexcept {
        for (size_type i = 0; i < capacity_; ++i) {
            if (entries_[i].state == State::Occupied) {
                EntryAllocTraits::destroy(alloc_, &entries_[i]);
                entries_[i].state = State::Empty;
            }
        }
        size_ = 0;
        deleted_ = 0;
    }

    void swap(RobinHoodHashMap& other) noexcept {
        std::swap(entries_, other.entries_);
        std::swap(capacity_, other.capacity_);
        std::swap(size_, other.size_);
        std::swap(deleted_, other.deleted_);
        std::swap(hasher_, other.hasher_);
        std::swap(key_equal_, other.key_equal_);
        if constexpr (EntryAllocTraits::propagate_on_container_swap::value) {
            std::swap(alloc_, other.alloc_);
        }
    }

    [[nodiscard]] hasher hash_function() const noexcept { return hasher_; }
    [[nodiscard]] key_equal key_eq() const noexcept { return key_equal_; }
    [[nodiscard]] allocator_type get_allocator() const noexcept { return alloc_; }
};

} // namespace detail

// Public API - khash_t compatible
template <typename Key,
          typename Value,
          typename Hash = std::hash<Key>,
          typename KeyEqual = std::equal_to<>,
          typename Allocator = memory::PolymorphicAllocator<std::pair<const Key, Value>>>
class HashMap : public detail::RobinHoodHashMap<Key, Value, Hash, KeyEqual, Allocator> {
    using Base = detail::RobinHoodHashMap<Key, Value, Hash, KeyEqual, Allocator>;
public:
    using Base::Base;
    using Base::operator=;
    using typename Base::iterator;
    using typename Base::const_iterator;
    using typename Base::size_type;
    using typename Base::key_type;
    using typename Base::mapped_type;
    using typename Base::value_type;

    // khash_t compatibility
    using khiter_t = iterator;
    using khint_t = size_type;

    // khash API
    [[nodiscard]] khiter_t kh_begin() noexcept { return this->begin(); }
    [[nodiscard]] khiter_t kh_end() noexcept { return this->end(); }
    [[nodiscard]] const_iterator kh_begin() const noexcept { return this->cbegin(); }
    [[nodiscard]] const_iterator kh_end() const noexcept { return this->cend(); }

    [[nodiscard]] bool kh_exist(khiter_t it) const noexcept { return it != this->end(); }
    [[nodiscard]] Key& kh_key(khiter_t it) noexcept { return it->first; }
    [[nodiscard]] const Key& kh_key(khiter_t it) const noexcept { return it->first; }
    [[nodiscard]] Value& kh_val(khiter_t it) noexcept { return it->second; }
    [[nodiscard]] const Value& kh_val(khiter_t it) const noexcept { return it->second; }

    template <typename K>
    [[nodiscard]] khiter_t kh_get(const K& k) noexcept { return this->find(k); }

    template <typename K>
    [[nodiscard]] khint_t kh_put(const K& k, int* ret = nullptr) {
        auto [it, inserted] = this->insert(k, Value{});
        if (ret) *ret = inserted ? 0 : 1;
        return inserted ? 0 : 1; // khash returns 0 for new, 1 for existing
    }

    template <typename K>
    khint_t kh_del(const K& k) {
        return this->erase(k) ? 0 : 1; // khash returns 0 for success
    }

    void kh_clear() noexcept { this->clear(); }
    void kh_destroy() noexcept { this->clear(); }
    void kh_resize(size_type n) { this->reserve(n); }
};

// HashSet - thin wrapper
template <typename Key,
          typename Hash = std::hash<Key>,
          typename KeyEqual = std::equal_to<>,
          typename Allocator = memory::PolymorphicAllocator<Key>>
class HashSet : public detail::RobinHoodHashMap<Key, std::byte, Hash, KeyEqual, Allocator> {
    using Base = detail::RobinHoodHashMap<Key, std::byte, Hash, KeyEqual, Allocator>;
public:
    using Base::Base;
    using typename Base::iterator;
    using typename Base::const_iterator;
    using typename Base::size_type;
    using typename Base::key_type;
    using typename Base::value_type;

    template <typename K>
    std::pair<iterator, bool> insert(K&& k) {
        return Base::insert_impl(std::forward<K>(k), std::byte{});
    }

    template <typename... Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        return insert(std::forward<Args>(args)...);
    }
};

// Convenience aliases
template <typename Key, typename Value>
using PmrHashMap = HashMap<Key, Value, std::hash<Key>, std::equal_to<>,
                           memory::PolymorphicAllocator<std::pair<const Key, Value>>>;

template <typename Key>
using PmrHashSet = HashSet<Key, std::hash<Key>, std::equal_to<>,
                           memory::PolymorphicAllocator<Key>>;

// String-specific optimizations
template <typename Value>
using StringHashMap = HashMap<PmrString, Value, std::hash<PmrString>, std::equal_to<>,
                              memory::PolymorphicAllocator<std::pair<const PmrString, Value>>>;

using StringHashSet = HashSet<core::PmrString, std::hash<core::PmrString>, std::equal_to<>,
                              memory::PolymorphicAllocator<core::PmrString>>;

} // namespace bwa::core