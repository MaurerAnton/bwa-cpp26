#include <bwa/align/sw.hpp>
#include <bwa/align/mem.hpp>
#include <bwa/index/fm_index.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/pipeline.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <array>
#include <initializer_list>

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

    // SW local / semi-global.
    // sw_semi_global aligns the query end-to-end with free reference flanks
    // (unlike sw_global, which penalizes them).
    {
        align::Scoring sc = align::Scoring::bwa_mem_default();
        std::array<uint8_t, 4> q = {0, 1, 2, 3};
        std::array<uint8_t, 8> r = {0, 0, 0, 1, 2, 3, 0, 0};
        align::Alignment a = align::sw_local(sc, q, r, 0);
        test("SW local positive", a.score > 0);
        align::Alignment b = align::sw_semi_global_extend(sc, q, r, 32);
        test("SW semi-global-extend positive", b.score > 0);
        // Exact fit: 4 matches at ref offsets [2, 6), flanks free
        align::Alignment c = align::sw_semi_global(sc, q, r, 0);
        test("SW semi-global score", c.score == 4);
        test("SW semi-global span",
             c.ref_begin == 2 && c.ref_end == 6 && c.query_begin == 0 &&
             c.query_end == 4);
        bool all_equal = !c.cigar.empty();
        for (uint32_t ci : c.cigar) {
            if (align::cigar_op(ci) != align::CigarOp::Equal) all_equal = false;
        }
        test("SW semi-global CIGAR", all_equal);
        // Banded extension agrees with full semi-global here
        test("SW semi-global banded agrees", b.score == c.score);
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

    // DP-optimal chaining on hand-built MEMs (no FM-index involved).
    {
        memory::Arena arena(1024 * 1024);
        index::PackedSequence ref;
        ref.append("ACGTACGTACGTACGTACGTACGT", 24);
        index::FMIndex idx = index::FMIndex::build(ref, arena);
        align::MEMFinder finder(idx, 4, 500, 1);

        auto make_mem = [](int32_t qp, int32_t rp, int32_t len,
                           int32_t score, int32_t ref = 0) {
            align::MEM m;
            m.query_pos = qp;
            m.ref_pos = rp;
            m.ref_id = ref;
            m.len = len;
            m.score = score;
            m.is_forward = true;
            return m;
        };
        auto make_vec = [&](std::initializer_list<align::MEM> ms) {
            core::Vector<align::MEM> v(&arena);
            for (const auto& m : ms) v.push_back(m);
            return v;
        };

        // 1. Colinear seeds join with gap penalty max(qgap, rgap)
        {
            auto mems = make_vec({make_mem(0, 100, 10, 10), make_mem(12, 112, 10, 10)});
            auto chains = finder.chain(mems, 10000, 5);
            test("DP colinear count", chains.size() == 1);
            test("DP colinear score", chains.size() == 1 && chains[0].score == 18);
            test("DP colinear members", chains.size() == 1 && chains[0].mems.size() == 2);
            test("DP colinear bounds",
                 chains.size() == 1 && chains[0].query_begin == 0 &&
                 chains[0].query_end == 22 && chains[0].ref_begin == 100 &&
                 chains[0].ref_end == 122 && chains[0].ref_id == 0);
        }

        // 2. Greedy-failure regression: a far jump must not drag a good
        // prefix below the threshold (greedy scored -50 and dropped both).
        {
            auto mems = make_vec({make_mem(0, 0, 10, 10), make_mem(10, 10, 10, 10),
                                  make_mem(20, 100, 10, 10)});
            auto chains = finder.chain(mems, 10000, 5);
            test("DP rescue count", chains.size() == 2);
            test("DP rescue best", chains.size() == 2 && chains[0].score == 20 &&
                 chains[0].mems.size() == 2);
            test("DP rescue second", chains.size() == 2 && chains[1].score == 10);
        }

        // 3. Same coordinates on different references stay separate
        {
            auto mems = make_vec({make_mem(0, 0, 10, 10, 0), make_mem(0, 0, 10, 10, 1)});
            auto chains = finder.chain(mems, 10000, 5);
            test("DP refsplit count", chains.size() == 2);
            test("DP refsplit ids",
                 chains.size() == 2 && chains[0].ref_id != chains[1].ref_id);
        }

        // 4. Query-overlapping seeds do not link (strict colinearity)
        {
            auto mems = make_vec({make_mem(0, 0, 10, 10), make_mem(5, 5, 10, 10)});
            auto chains = finder.chain(mems, 10000, 5);
            test("DP overlap split", chains.size() == 2);
        }

        // 5. Empty input, empty output
        {
            auto mems = make_vec({});
            auto chains = finder.chain(mems, 10000, 5);
            test("DP empty", chains.size() == 0);
        }
    }

    // Multi-reference alignment: reads must map to the correct contig.
    // Regression test: the aligner previously searched only fm_index()[0],
    // making all other references invisible.
    {
        const char* fa_path = "/tmp/bwa_test_multiref.fa";
        std::string chrA_seq, chrB_seq;
        for (int i = 0; i < 32; ++i) { chrA_seq += "ACGT"; chrB_seq += "TGCA"; }
        {
            std::ofstream fa(fa_path);
            fa << ">chrA\n" << chrA_seq << "\n";
            fa << ">chrB\n" << chrB_seq << "\n";
        }
        Index idx = Index::build(fa_path);
        test("MultiRef num references", idx.num_references() == 2);
        Aligner aligner(idx);

        io::SeqRecord readA, readB;
        readA.name = "readA";
        readA.seq = "ACGTACGTACGTACGTACGT";  // chrA only
        readA.qual = "IIIIIIIIIIIIIIIIIIII";
        readB.name = "readB";
        readB.seq = "TGCATGCATGCATGCATGCA";  // chrB only
        readB.qual = "IIIIIIIIIIIIIIIIIIII";

        AlignmentResult resA = aligner.align(readA);
        test("MultiRef readA mapped", resA.mapped);
        test("MultiRef readA rname", resA.primary.rname == "chrA");
        test("MultiRef readA pos", resA.primary.pos >= 1 &&
             resA.primary.pos <= 128 - 20 + 1);
        test("MultiRef readA cigar", !resA.primary.cigar.empty());

        AlignmentResult resB = aligner.align(readB);
        test("MultiRef readB mapped", resB.mapped);
        test("MultiRef readB rname", resB.primary.rname == "chrB");
        test("MultiRef readB pos", resB.primary.pos >= 1 &&
             resB.primary.pos <= 128 - 20 + 1);

        std::remove(fa_path);
    }

    // N-containing reference: flanking matches must survive the N-run,
    // and the save/load roundtrip must preserve N-ness (v2 format).
    {
        const char* fa_path = "/tmp/bwa_test_nref.fa";
        const char* idx_prefix = "/tmp/bwa_test_nref_idx";
        {
            std::ofstream fa(fa_path);
            fa << ">chrN\nACGTACGTACGTACGTNNNNNNNNGATTACAGATTACA\n";
        }
        Index idx = Index::build(fa_path);
        idx.save(idx_prefix);
        Index loaded = Index::load(idx_prefix);
        test("NRef roundtrip refs", loaded.num_references() == 1);
        Aligner aligner(loaded);

        io::SeqRecord readL, readR, readN;
        readL.name = "readL";
        readL.seq = "GTACGTACGTAC";  // left flank, offset 2
        readL.qual = "IIIIIIIIIIII";
        readR.name = "readR";
        readR.seq = "TTACAGATTACA";  // right flank
        readR.qual = "IIIIIIIIIIII";
        readN.name = "readN";
        readN.seq = "NNNNNNNNNNNN";  // pure N: must not map
        readN.qual = "IIIIIIIIIIII";

        AlignmentResult resL = aligner.align(readL);
        test("NRef left mapped", resL.mapped);
        test("NRef left rname", resL.primary.rname == "chrN");
        // Perfect 12-mer: only a full-match placement scores 12
        test("NRef left perfect", resL.best_score == 12);
        test("NRef left pos", resL.primary.pos == 3);

        AlignmentResult resR = aligner.align(readR);
        test("NRef right mapped", resR.mapped);
        test("NRef right rname", resR.primary.rname == "chrN");
        // Unique in right flank (0-based 26) -> 1-based 27
        test("NRef right perfect", resR.best_score == 12);
        test("NRef right pos", resR.primary.pos == 27);

        AlignmentResult resN = aligner.align(readN);
        test("NRef N-read unmapped", !resN.mapped);

        std::remove(fa_path);
        std::remove("/tmp/bwa_test_nref_idx.meta");
        std::remove("/tmp/bwa_test_nref_idx.bwt");
        std::remove("/tmp/bwa_test_nref_idx.sa");
        std::remove("/tmp/bwa_test_nref_idx.occ");
        std::remove("/tmp/bwa_test_nref_idx.pac");
    }

    // MAPQ model: faithful port of BWA mem_approx_mapq_se (bwamem.c).
    // Hand-computed expectations (match=1, mismatch=4):
    //   unique 20-mer, min_seed 9: 6.02*(20-9) = 66.2 -> cap 60
    //   unique 20-mer, min_seed 19 (BWA default): 6.02*1 = 6
    //   tie 20/20: 0; sub>score: 0
    //   100bp 100-vs-90: tmp=log50/log100=0.8495, 6.02*10*0.7216=43
    //   same + 1 suboptimal hit: 43-(4.343*log2)=40
    //   100bp weak 30-vs-9, identity 0.86: 6.02*21*0.7396^2... = 50
    //   chain-suboptimal dominates sub: 6.02*(20-18) = 12
    {
        auto sig = [](int32_t score, int32_t sub, int32_t span, int32_t min_seed = 9,
                      int32_t sub_n = 0, int32_t csub = 0) {
            Aligner::MapqSignals s;
            s.score = score;
            s.sub = sub;
            s.csub = csub;
            s.sub_n = sub_n;
            s.seedcov = span;
            s.query_span = span;
            s.ref_span = span;
            s.match = 1;
            s.mismatch = 4;
            s.min_seed_len = min_seed;
            return s;
        };
        test("MAPQ unique capped", Aligner::approx_mapq_se(sig(20, 0, 20)) == 60);
        test("MAPQ BWA-default seed", Aligner::approx_mapq_se(sig(20, 0, 20, 19)) == 6);
        test("MAPQ tie", Aligner::approx_mapq_se(sig(20, 20, 20)) == 0);
        test("MAPQ sub exceeds", Aligner::approx_mapq_se(sig(10, 15, 20)) == 0);
        test("MAPQ long near-tie", Aligner::approx_mapq_se(sig(100, 90, 100)) == 43);
        test("MAPQ sub_n penalty",
             Aligner::approx_mapq_se(sig(100, 90, 100, 9, 1)) == 40);
        test("MAPQ weak identity", Aligner::approx_mapq_se(sig(30, 0, 100)) == 50);
        test("MAPQ csub dominates",
             Aligner::approx_mapq_se(sig(20, 5, 20, 9, 0, 18)) == 12);
        // Monotonicity: more suboptimal hits never raise MAPQ
        uint8_t q0 = Aligner::approx_mapq_se(sig(100, 90, 100, 9, 0));
        uint8_t q3 = Aligner::approx_mapq_se(sig(100, 90, 100, 9, 3));
        test("MAPQ monotone sub_n", q3 <= q0);
    }

    std::cout << "\nResults: " << passed << " passed, " << failed << " failed\n";
    return failed ? 1 : 0;
}
