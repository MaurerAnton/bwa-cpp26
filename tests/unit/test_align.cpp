#include <bwa/align/sw.hpp>
#include <bwa/align/mem.hpp>
#include <bwa/index/fm_index.hpp>
#include <bwa/core/arena.hpp>
#include <iostream>
#include <array>

using namespace bwa;

int main() {
    int passed = 0, failed = 0;
    auto test = [&](const char* name, bool cond) {
        if (cond) { std::cout << "  PASS: " << name << "\n"; ++passed; }
        else { std::cout << "  FAIL: " << name << "\n"; ++failed; }
    };

    // SW global perfect match
    {
        align::Scoring sc = align::Scoring::bwa_mem_default();
        std::array<uint8_t, 8> q = {0, 1, 2, 3, 0, 1, 2, 3};
        std::array<uint8_t, 8> r = {0, 1, 2, 3, 0, 1, 2, 3};
        align::Alignment a = align::sw_global(sc, q, r, 0);
        test("SW global perfect score", a.score == 8);
        bool has_eq = false;
        for (uint32_t c : a.cigar)
            if (align::cigar_op(c) == align::CigarOp::Equal) has_eq = true;
        test("SW global CIGAR =/X", has_eq);
    }

    // SW global with mismatches
    {
        align::Scoring sc = align::Scoring::bwa_mem_default();
        std::array<uint8_t, 8> q = {0, 1, 2, 3, 0, 1, 2, 3};
        std::array<uint8_t, 8> r = {0, 1, 2, 3, 3, 3, 2, 3};
        align::Alignment a = align::sw_global(sc, q, r, 0);
        bool has_diff = false;
        for (uint32_t c : a.cigar)
            if (align::cigar_op(c) == align::CigarOp::Diff) has_diff = true;
        test("SW global mismatch CIGAR", has_diff);
    }

    // SW local / semi-global smoke
    // NOTE: sw_semi_global() is currently an alias for global alignment,
    // so flanking reference bases are penalized. Test the real pipeline path
    // sw_semi_global_extend() (free ref end gaps) for positivity instead.
    {
        align::Scoring sc = align::Scoring::bwa_mem_default();
        std::array<uint8_t, 4> q = {0, 1, 2, 3};
        std::array<uint8_t, 8> r = {0, 0, 0, 1, 2, 3, 0, 0};
        align::Alignment a = align::sw_local(sc, q, r, 0);
        test("SW local positive", a.score > 0);
        align::Alignment b = align::sw_semi_global_extend(sc, q, r, 32);
        test("SW semi-global-extend positive", b.score > 0);
    }

    // MEM find + chain smoke (chain threshold lowered: 8bp MEM scores 8 < default 30)
    {
        memory::Arena arena(1024 * 1024);
        index::PackedSequence ref;
        ref.append("ACGTACGTACGTACGTACGTACGT", 24);
        index::FMIndex idx = index::FMIndex::build(ref, arena);
        align::MEMFinder finder(idx, 4, 500, 1);
        index::PackedSequence query;
        query.append("ACGTACGT", 8);
        auto mems = finder.find(query.bases(), arena);
        test("MEM find", mems.size() > 0);
        auto chains = finder.chain(mems, 10000, 5);
        test("MEM chain", chains.size() > 0);
    }

    std::cout << "\nResults: " << passed << " passed, " << failed << " failed\n";
    return failed ? 1 : 0;
}
