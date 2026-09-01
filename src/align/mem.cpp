#include <bwa/align/mem.hpp>
#include <bwa/index/fm_index.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/vector.hpp>
#include <algorithm>

namespace bwa::align {

Vector<MEM> MEMFinder::find(const std::span<const uint8_t>& query,
                             memory::Arena& arena) const {
    Vector<MEM> mems(arena);
    find(query, mems);
    return mems;
}

void MEMFinder::find(const std::span<const uint8_t>& query,
                     Vector<MEM>& mems) const {
    find_strand(query, mems, true);
}

// Additional: bidirectional MEM finding
Vector<MEM> find_mems_bidirectional(const index::FMIndex& fwd_idx,
                                     const index::FMIndex& rev_idx,
                                     const std::span<const uint8_t>& query,
                                     int min_len,
                                     int max_occ,
                                     memory::Arena& arena) {
    Vector<MEM> mems(arena);

    MEMFinder fwd_finder(fwd_idx, min_len, max_occ, 1);
    MEMFinder rev_finder(rev_idx, min_len, max_occ, 1);

    fwd_finder.find(query, mems);

    // Reverse complement query and find on reverse index
    PackedSequence rc_query;
    rc_query.resize(query.size());
    for (size_t i = 0; i < query.size(); ++i) {
        uint8_t b = query[query.size() - 1 - i];
        rc_query.set(i, b == 0 ? 3 : b == 3 ? 0 : b == 1 ? 2 : b == 2 ? 1 : 4);
    }

    Vector<MEM> rev_mems(arena);
    rev_finder.find(rc_query.words(), rev_mems);

    // Convert reverse positions back to forward coordinates
    for (MEM& mem : rev_mems) {
        mem.is_forward = false;
        // Adjust reference position for reverse strand
        // mem.ref_pos = fwd_idx.length() - (mem.ref_pos + mem.len);
        mems.push_back(mem);
    }

    return mems;
}

} // namespace bwa::align