#include <bwa/index/fm_index.hpp>
#include <bwa/align/mem.hpp>
#include <bwa/io/seq_io.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/string.hpp>
#include <iostream>
#include <fstream>

using namespace bwa;
using namespace bwa::index;

int main() {
    int passed = 0, failed = 0;

    auto test = [&](const char* name, bool cond) {
        if (cond) { std::cout << "  PASS: " << name << "\n"; ++passed; }
        else { std::cout << "  FAIL: " << name << "\n"; ++failed; }
    };

    // Build small index from N-free sequence
    {
        memory::Arena arena(10 * 1024 * 1024);
        PackedSequence seq;
        seq.append("ACGTACGTACGTACGTACGTACGTAC", 26);

        FMIndex idx = FMIndex::build(seq, arena);
        test("FMIndex build", idx.size() == 26);
        // BWT has n+1 rows (including sentinel row)
        test("FMIndex BWT size", idx.bwt_size() == 27);

        // Test rank
        test("FMIndex rank", idx.rank(0, 10) >= 0);

        // Test C array
        test("FMIndex C array", idx.C(0) == 0);

        // Test backward search
        PackedSequence query;
        query.append("ACGT", 4);
        auto [l, r] = idx.backward_search(query);
        test("FMIndex backward search", l < r);

        // Test count
        size_t cnt = idx.count(query);
        test("FMIndex count", cnt > 0);

        // Test locate
        if (l < r) {
            auto pos = idx.locate(l);
            test("FMIndex locate", pos.has_value());
        }
    }

    // N-containing reference: Ns are separators for extension, but exact
    // ACGT matches are found wherever they occur (positions 0,4,12,16);
    // query Ns never match
    {
        memory::Arena arena(10 * 1024 * 1024);
        PackedSequence seq;
        seq.append("ACGTACGTNNNNACGTACGT", 20);

        FMIndex idx = FMIndex::build(seq, arena);
        test("FMIndex N build", idx.size() == 20);
        // Total counts include sentinel + all 5 symbols
        test("FMIndex N counts", idx.total(0) + idx.total(1) + idx.total(2) +
             idx.total(3) + idx.total(4) + idx.total(5) == 21);

        PackedSequence query;
        query.append("ACGT", 4);
        auto [l, r] = idx.backward_search(query);
        test("FMIndex N backward search", l < r);
        test("FMIndex N count", idx.count(query) == static_cast<size_t>(r - l));
        // All 4 exact occurrences (neighbors don't matter for exact search)
        test("FMIndex N count value", r - l == 4);
        if (l < r) {
            auto pos = idx.locate(l);
            test("FMIndex N locate", pos.has_value() && *pos < 20);
            // Every reported occurrence must really match
            bool all_match = true;
            for (size_t p = l; p < r; ++p) {
                auto q = idx.locate(p);
                if (!q || *q + 4 > 20) { all_match = false; break; }
                for (int k = 0; k < 4; ++k) {
                    if (seq.get(*q + k) != query.get(k)) { all_match = false; break; }
                }
                if (!all_match) break;
            }
            test("FMIndex N locate correct", all_match && (r - l) == 4);
        }

        PackedSequence nquery;
        nquery.append("NNNN", 4);
        auto [nl, nr] = idx.backward_search(nquery);
        test("FMIndex N query empty", nl >= nr);
    }

    // MultiFMIndex
    {
        memory::Arena arena(10 * 1024 * 1024);
        MultiFMIndex multi;

        PackedSequence seq1;
        seq1.append("AAAAACCCCC", 10);
        multi.add_sequence(seq1, arena);

        PackedSequence seq2;
        seq2.append("GGGGGTTTTT", 10);
        multi.add_sequence(seq2, arena);

        test("MultiFMIndex num seq", multi.num_sequences() == 2);
        test("MultiFMIndex total len", multi.total_length() == 20);

        auto [idx, pos] = multi.locate_sequence(5);
        test("MultiFMIndex locate seq 0", idx == 0 && pos == 5);
        auto [idx2, pos2] = multi.locate_sequence(15);
        test("MultiFMIndex locate seq 1", idx2 == 1 && pos2 == 5);
    }

    // SMEM Iterator
    {
        memory::Arena arena(10 * 1024 * 1024);
        PackedSequence seq;
        seq.append("ACGTACGTACGTACGTACGTACGTACGTAC", 30);
        FMIndex idx = FMIndex::build(seq, arena);

        PackedSequence query;
        query.append("ACGTACGT", 8);
        align::SMEMIterator iter(idx, query.bases(), 4, 500);

        int count = 0;
        while (iter.has_next()) {
            auto mem = iter.next();
            if (mem) ++count;
        }
        test("SMEMIterator", count > 0);
    }

    // MEMFinder chaining
    {
        memory::Arena arena(10 * 1024 * 1024);
        PackedSequence seq;
        seq.append("AAAAACCCCCGGGGGTTTTTAAAAACCCCC", 30);
        FMIndex idx = FMIndex::build(seq, arena);

        align::MEMFinder finder(idx, 4, 500, 1);
        PackedSequence query;
        query.append("AAAAACCCCC", 10);

        auto mems = finder.find(query.bases(), arena);
        test("MEMFinder find", mems.size() > 0);

        auto chains = finder.chain(mems, 10000, 5);
        test("MEMFinder chain", chains.size() > 0);
        if (!chains.empty()) {
            test("Chain score", chains[0].score > 0);
        }
    }

    std::cout << "\nResults: " << passed << " passed, " << failed << " failed\n";
    return failed ? 1 : 0;
}