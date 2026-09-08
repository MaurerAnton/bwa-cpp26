#pragma once

#include <bwa/core/vector.hpp>
#include <cstdint>
#include <cstddef>

// Forward declaration
namespace bwa::index {
    class PackedSequence;
}

namespace bwa::index::detail {

// SA-IS: Linear-time Suffix Array Construction
// Based on: Nong, Zhang, Chan (2009) "Two Efficient Algorithms for Linear Time Suffix Array Construction"
// Implementation reference: https://github.com/algorithm-ninja/sais

namespace sais {

// Build suffix array using SA-IS algorithm - O(n) time
core::Vector<uint32_t> build_suffix_array(const PackedSequence& seq,
                                          bwa::memory::Arena& arena) noexcept;

} // namespace sais

} // namespace bwa::index::detail