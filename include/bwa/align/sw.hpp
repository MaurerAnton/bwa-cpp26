#pragma once

#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/arena.hpp>
#include <span>
#include <cstdint>
#include <array>
#include <algorithm>
#include <limits>
#include <concepts>

namespace bwa::align {

// Scoring parameters
struct Scoring {
    int match = 1;       // Match score
    int mismatch = -4;   // Mismatch penalty
    int gap_open = -6;   // Gap open penalty
    int gap_ext = -1;    // Gap extension penalty
    int clip_pen = -5;   // Clipping penalty (for local alignment)
    int unpaired = -17;  // Unpaired penalty (for paired-end)

    [[nodiscard]] static Scoring bwa_mem_default() noexcept {
        return {1, -4, -6, -1, -5, -17};
    }

    [[nodiscard]] static Scoring bwa_mem_online() noexcept {
        return {1, -6, -6, -1, -5, -17};
    }
};

// Alignment result
struct Alignment {
    int32_t score = 0;          // Alignment score
    int32_t ref_begin = 0;      // Reference start (0-based)
    int32_t ref_end = 0;        // Reference end (exclusive)
    int32_t query_begin = 0;    // Query start (0-based)
    int32_t query_end = 0;      // Query end (exclusive)
    int32_t n_cigar = 0;        // Number of CIGAR operations
    std::vector<uint32_t> cigar;     // CIGAR operations (encoded as (len<<4)|op)
    int32_t n_mismatch = 0;     // Number of mismatches
    int32_t n_gap_open = 0;     // Number of gap opens
    int32_t n_gap_ext = 0;      // Number of gap extensions

    void clear() noexcept {
        score = 0;
        ref_begin = ref_end = 0;
        query_begin = query_end = 0;
        n_cigar = 0;
        cigar.clear();
        n_mismatch = n_gap_open = n_gap_ext = 0;
    }

    [[nodiscard]] bool empty() const noexcept { return n_cigar == 0; }
};

// CIGAR operations
enum class CigarOp : uint8_t {
    Match = 0,   // M (alignment match)
    Ins = 1,     // I (insertion to reference)
    Del = 2,     // D (deletion from reference)
    Skip = 3,    // N (skipped region from reference)
    SoftClip = 4, // S (soft clipping)
    HardClip = 5, // H (hard clipping)
    Pad = 6,     // P (padding)
    Equal = 7,   // = (sequence match)
    Diff = 8     // X (sequence mismatch)
};

[[nodiscard]] inline constexpr uint32_t encode_cigar(int len, CigarOp op) noexcept {
    return (static_cast<uint32_t>(len) << 4) | static_cast<uint32_t>(op);
}

[[nodiscard]] inline constexpr int cigar_len(uint32_t cigar) noexcept {
    return static_cast<int>(cigar >> 4);
}

[[nodiscard]] inline constexpr CigarOp cigar_op(uint32_t cigar) noexcept {
    return static_cast<CigarOp>(cigar & 0xF);
}

[[nodiscard]] inline char cigar_char(CigarOp op) noexcept {
    static constexpr char MAP[] = "MIDNSHP=X";
    return MAP[static_cast<uint8_t>(op)];
}

// SIMD-accelerated Smith-Waterman (global, local, extension)
// Replaces ksw.c with std::simd (C++26) or Highway fallback
namespace detail {

// 16-bit score matrix for SIMD
using score_t = int16_t;
// Placeholder for SIMD type - will use std::simd when available
using simd_score_t = std::array<int16_t, 8>;

// SIMD max for 8 lanes
inline simd_score_t simd_max_8(simd_score_t a, simd_score_t b) noexcept {
    simd_score_t r;
    for (int i = 0; i < 8; ++i) r[i] = (a[i] > b[i]) ? a[i] : b[i];
    return r;
}

// Load 8 int16 from memory
inline simd_score_t simd_load_8(const int16_t* p) noexcept {
    simd_score_t v;
    for (int i = 0; i < 8; ++i) v[i] = p[i];
    return v;
}

// Store 8 int16 to memory
inline void simd_store_8(int16_t* p, simd_score_t v) noexcept {
    for (int i = 0; i < 8; ++i) p[i] = v[i];
}

} // namespace detail

// Banded Smith-Waterman for alignment extension
// query: query sequence (packed 2-bit)
// ref: reference sequence (packed 2-bit)
// Returns optimal alignment within band
Alignment sw_extend(const Scoring& sc,
                    std::span<const uint8_t> query,
                    std::span<const uint8_t> ref,
                    int32_t query_start, int32_t ref_start,
                    int32_t band_width = 32,
                    int32_t max_score_drop = 100);

// Global alignment (Needleman-Wunsch with affine gaps)
Alignment sw_global(const Scoring& sc,
                    std::span<const uint8_t> query,
                    std::span<const uint8_t> ref,
                    int32_t band_width = 0); // 0 = full DP

// Local alignment (Smith-Waterman with affine gaps)
Alignment sw_local(const Scoring& sc,
                   std::span<const uint8_t> query,
                   std::span<const uint8_t> ref,
                   int32_t band_width = 0);

// Semi-global: query aligned globally, reference locally (for extension)
Alignment sw_semi_global(const Scoring& sc,
                         std::span<const uint8_t> query,
                         std::span<const uint8_t> ref,
                         int32_t band_width = 0);

// Batch alignment for multiple queries
void sw_batch(const Scoring& sc,
              std::span<const std::span<const uint8_t>> queries,
              std::span<const uint8_t> ref,
              std::span<Alignment> results,
              int32_t band_width = 32);

// SIMD-accelerated batch (8 queries at a time)
void sw_batch_simd(const Scoring& sc,
                   std::span<const std::span<const uint8_t>> queries,
                   std::span<const uint8_t> ref,
                   std::span<Alignment> results,
                   int32_t band_width = 32);

} // namespace bwa::align