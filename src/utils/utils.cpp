#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/index/fm_index.hpp>
#include <chrono>
#include <random>
#include <iostream>

namespace bwa::utils {

// Timer utility
class Timer {
    std::chrono::high_resolution_clock::time_point start_;
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}
    void reset() { start_ = std::chrono::high_resolution_clock::now(); }
    [[nodiscard]] double elapsed_ms() const noexcept {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }
    [[nodiscard]] double elapsed_s() const noexcept {
        return elapsed_ms() / 1000.0;
    }
};

// Random sequence generator
class RandomSequence {
    std::mt19937_64 rng_;
    std::uniform_int_distribution<int> dist_{0, 3};
public:
    explicit RandomSequence(uint64_t seed = 42) : rng_(seed) {}

    void generate(index::PackedSequence& seq, size_t len, double n_rate = 0.0) {
        seq.resize(len);
        for (size_t i = 0; i < len; ++i) {
            if (n_rate > 0 && (rng_() % 10000) < static_cast<int>(n_rate * 10000)) {
                seq.set(i, 4); // N
            } else {
                seq.set(i, dist_(rng_));
            }
        }
    }

    void generate_string(core::PmrString& s, size_t len, double n_rate = 0.0) {
        s.clear();
        s.reserve(len);
        static constexpr char BASES[] = "ACGTN";
        for (size_t i = 0; i < len; ++i) {
            if (n_rate > 0 && (rng_() % 10000) < static_cast<int>(n_rate * 10000)) {
                s.kputc('N');
            } else {
                s.kputc(BASES[dist_(rng_)]);
            }
        }
    }
};

// Memory usage reporter
struct MemoryStats {
    size_t arena_used = 0;
    size_t arena_capacity = 0;
    size_t vectors_total = 0;
    size_t strings_total = 0;
    size_t hash_maps_total = 0;
};

MemoryStats get_memory_stats() {
    MemoryStats stats;
    // Would iterate over all arenas in real implementation
    return stats;
}

// Format bytes as human-readable
core::PmrString format_bytes(size_t bytes) {
    core::PmrString s;
    const char* units[] = {"B", "KB", "MB", "GB", "TB"};
    int unit = 0;
    double val = static_cast<double>(bytes);
    while (val >= 1024 && unit < 4) {
        val /= 1024;
        ++unit;
    }
    s.ksprintf("{:.2f} {}", val, units[unit]);
    return s;
}

} // namespace bwa::utils