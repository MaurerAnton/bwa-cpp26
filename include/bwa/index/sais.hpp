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
// CONTRACT (must hold before re-enabling in FMIndex::build, which currently
// uses the verified brute-force builder): return all n+1 suffixes INCLUDING
// the empty suffix (value n) sorted first, with codes
// sentinel=0 < N=1 < A=2 < C=3 < G=4 < T=5.
core::Vector<uint32_t> build_suffix_array(const PackedSequence& seq,
                                          bwa::memory::Arena& arena) noexcept;

} // namespace sais

} // namespace bwa::index::detail