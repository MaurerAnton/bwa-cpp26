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

    // SW extension with clipping (BWA ksw_extend contract): a query with an
    // unalignable prefix/suffix must align only its core and report the
    // aligned query span, not fail or force end-to-end alignment.
    {
        align::Scoring sc = align::Scoring::bwa_mem_default();
        // R = 21 bases; query = 19 A's + R. A's are absent from R, so the
        // only positive-scoring alignment is R itself.
        std::vector<uint8_t> R = {2, 0, 1, 2, 3, 2, 3, 0, 1, 1, 0,
                                  2, 2, 1, 0, 3, 1, 2, 2, 0, 3};
        std::vector<uint8_t> Q(19, 0);
        Q.insert(Q.end(), R.begin(), R.end());
        auto a = align::sw_semi_global_extend(
            sc, std::span<const uint8_t>(Q.data(), Q.size()),
            std::span<const uint8_t>(R.data(), R.size()), 32);
        test("SW clip score", a.score == 21);
        test("SW clip query span",
             a.query_begin == 19 && a.query_end == 40);
        test("SW clip ref span", a.ref_begin == 0 && a.ref_end == 21);
        bool has_soft = false;
        for (uint32_t ci : a.cigar) {
            if (align::cigar_op(ci) == align::CigarOp::SoftClip) has_soft = true;
        }
        test("SW clip CIGAR has S", has_soft);

        // Wide reference window (flanks on both sides) must find the same
        // core: this regressed when the banded DP was indexed wrongly.
        std::vector<uint8_t> wide(20, 1);
        wide.insert(wide.end(), R.begin(), R.end());
        wide.insert(wide.end(), 20, 3);
        auto w = align::sw_semi_global_extend(
            sc, std::span<const uint8_t>(Q.data(), Q.size()),
            std::span<const uint8_t>(wide.data(), wide.size()), 32);
        test("SW wide-window score", w.score == 21);
        test("SW wide-window span",
             w.query_begin == 19 && w.query_end == 40 && w.ref_begin == 20 &&
             w.ref_end == 41);
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

    // Paired-end: FR orientation gives PROPER_PAIR with signed TLEN and
    // '=' RNEXT; an unmapped mate yields UNMAP/MUNMAP handling.
    {
        const char* fa_path = "/tmp/bwa_test_pair.fa";
        const std::string r1seg = "AATTACATAACATACACGTC";  // pos 31 (1-based)
        const std::string r2fwd = "GCTGTGTCCACCCCATCGGA";  // pos 121 (1-based)
        // RC(r2fwd) maps reverse; fillers are homopolymers without 9-mer
        // seeds from either segment.
        {
            std::ofstream fa(fa_path);
            fa << ">chrP\n" << std::string(30, 'G') << r1seg
               << std::string(70, 'C') << r2fwd << std::string(30, 'G') << "\n";
        }
        Index idx = Index::build(fa_path);
        Aligner aligner(idx);

        io::SeqRecord pe1, pe2;
        pe1.name = "pe1/1";
        pe1.seq = r1seg;
        pe1.qual = "IIIIIIIIIIIIIIIIIIII";
        pe2.name = "pe1/2";
        // reverse complement of r2fwd
        std::string r2rc;
        for (auto it = r2fwd.rbegin(); it != r2fwd.rend(); ++it) {
            char c = *it;
            r2rc += (c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A');
        }
        pe2.seq = r2rc;
        pe2.qual = "IIIIIIIIIIIIIIIIIIII";

        auto [pres1, pres2] = aligner.align_pair(pe1, pe2);
        test("Pair both mapped", pres1.mapped && pres2.mapped);
        test("Pair read1 forward proper",
             (pres1.primary.flag &
              (AlnRecord::F_PAIRED | AlnRecord::F_PROPER_PAIR | AlnRecord::F_READ1)) ==
                 (AlnRecord::F_PAIRED | AlnRecord::F_PROPER_PAIR | AlnRecord::F_READ1));
        test("Pair read1 not reversed",
             (pres1.primary.flag & AlnRecord::F_REVERSE) == 0);
        test("Pair read1 pos", pres1.primary.pos == 31);
        test("Pair read2 reversed proper",
             (pres2.primary.flag &
              (AlnRecord::F_PAIRED | AlnRecord::F_PROPER_PAIR |
               AlnRecord::F_READ2 | AlnRecord::F_REVERSE)) ==
                 (AlnRecord::F_PAIRED | AlnRecord::F_PROPER_PAIR |
                  AlnRecord::F_READ2 | AlnRecord::F_REVERSE));
        test("Pair read2 pos", pres2.primary.pos == 121);
        test("Pair rnext", pres1.primary.rnext == "=" && pres2.primary.rnext == "=");
        test("Pair pnext",
             pres1.primary.pnext == 121 && pres2.primary.pnext == 31);
        // 5' ends: 0-based 30 fwd, 139 rev -> inclusive TLEN 110
        test("Pair tlen", pres1.primary.tlen == 110 && pres2.primary.tlen == -110);

        // Unmapped mate (pure N): stays emitted with UNMAP, mate gets MUNMAP.
        io::SeqRecord peN;
        peN.name = "peU/2";
        peN.seq = "NNNNNNNNNNNNNNNNNNNN";
        peN.qual = "IIIIIIIIIIIIIIIIIIII";
        auto [ures1, ures2] = aligner.align_pair(pe1, peN);
        test("Pair-N first mapped", ures1.mapped && !ures2.mapped);
        test("Pair-N mate unmapped flag",
             (ures1.primary.flag & AlnRecord::F_MUNMAP) != 0);
        test("Pair-N unmapped record",
             (ures2.primary.flag &
              (AlnRecord::F_PAIRED | AlnRecord::F_UNMAP | AlnRecord::F_READ2)) ==
                 (AlnRecord::F_PAIRED | AlnRecord::F_UNMAP | AlnRecord::F_READ2));
        test("Pair-N unmapped star",
             ures2.primary.rname == "*" && ures2.primary.cigar.empty());
        test("Pair-N mate coords",
             ures2.primary.rnext == "chrP" && ures2.primary.pnext == ures1.primary.pos);

        std::remove(fa_path);
    }

    // Mate rescue: a long end too diverged to seed (errors every 4th base
    // defeat exact and 1-mismatch 9-mer seeds) must still be recovered by
    // targeted extension around its mapped mate.
    {
        const char* fa_path = "/tmp/bwa_test_rescue.fa";
        // Deterministic pseudo-random reference, 10000 bp.
        std::string ref;
        ref.reserve(10000);
        uint32_t x = 987654321u;
        for (int i = 0; i < 10000; ++i) {
            x = x * 1664525u + 1013904223u;
            ref += "ACGT"[(x >> 16) & 3];
        }
        {
            std::ofstream fa(fa_path);
            fa << ">chrR\n" << ref << "\n";
        }
        Index idx = Index::build(fa_path);
        Aligner aligner(idx);

        // Mate: perfect 150-mer at offset 1000 (forward).
        io::SeqRecord mate;
        mate.name = "rescue/1";
        mate.seq = ref.substr(1000, 150);
        mate.qual = std::string(150, 'I');

        // Diverged 300-mer from offset 1500: RC (reverse strand) with every
        // 4th base substituted, so no 9-mer seed (exact or 1-mismatch) hits.
        std::string seg = ref.substr(1500, 300);
        std::string rc;
        rc.reserve(300);
        for (auto it = seg.rbegin(); it != seg.rend(); ++it) {
            char c = *it;
            rc += (c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A');
        }
        for (size_t i = 0; i < rc.size(); i += 4) {
            rc[i] = rc[i] == 'A' ? 'C' : 'A';
        }
        io::SeqRecord div;
        div.name = "rescue/2";
        div.seq = rc;
        div.qual = std::string(300, 'I');

        // Sanity: the diverged end alone must not map (else this isn't a
        // rescue test at all).
        AlignmentResult solo = aligner.align(div);
        test("Rescue solo unmapped", !solo.mapped);

        auto [rres1, rres2] = aligner.align_pair(mate, div);
        test("Rescue mate mapped", rres1.mapped);
        test("Rescue recovered", rres2.mapped);
        if (rres2.mapped) {
            test("Rescue rname", rres2.primary.rname == "chrR");
            // Reverse-strand placement at the true locus (1-based 1501).
            test("Rescue reverse",
                 (rres2.primary.flag & AlnRecord::F_REVERSE) != 0);
            test("Rescue pos",
                 rres2.primary.pos >= 1490 && rres2.primary.pos <= 1510);
            // FR orientation, correct order: proper pair. TLEN is 799, not
            // 800, because the mutated first base soft-clips (1S), shifting
            // the reverse 5' end one base left.
            test("Rescue proper",
                 (rres1.primary.flag & AlnRecord::F_PROPER_PAIR) != 0 &&
                 (rres2.primary.flag & AlnRecord::F_PROPER_PAIR) != 0);
            test("Rescue tlen",
                 rres1.primary.tlen == 799 && rres2.primary.tlen == -799);
        }

        std::remove(fa_path);
    }

    // Chimeric (split) read: two segments from different references must
    // produce a soft-clipped primary plus a hard-clipped supplementary with
    // reciprocal SA:Z tags.
    {
        const char* fa_path = "/tmp/bwa_test_chimera.fa";
        const std::string segA = "AATTACATAACATACACGTC";  // chr1 pos 31
        const std::string segB = "GCTGTGTCCACCCCATCGGA";  // chr2 pos 31
        {
            std::ofstream fa(fa_path);
            fa << ">chrX\n" << std::string(30, 'G') << segA
               << std::string(30, 'G') << "\n";
            fa << ">chrY\n" << std::string(30, 'C') << segB
               << std::string(30, 'C') << "\n";
        }
        Index idx = Index::build(fa_path);
        Aligner aligner(idx);

        io::SeqRecord chim;
        chim.name = "chim1";
        chim.seq = segA + segB;
        chim.qual = std::string(40, 'I');

        AlignmentResult res = aligner.align(chim);
        test("Chimera mapped", res.mapped);
        test("Chimera supplementary", res.supplementary.size() == 1);
        test("Chimera primary clip",
             res.primary.cigar.size() == 2 &&
             align::cigar_op(res.primary.cigar.back()) ==
                 align::CigarOp::SoftClip);
        if (res.supplementary.size() == 1) {
            const AlnRecord& supp = res.supplementary[0];
            test("Chimera supp flag",
                 (supp.flag & AlnRecord::F_SUPPLEMENTARY) != 0 &&
                 (supp.flag & AlnRecord::F_SECONDARY) == 0);
            test("Chimera supp ref", supp.rname != res.primary.rname);
            test("Chimera supp hard clip",
                 !supp.cigar.empty() &&
                 align::cigar_op(supp.cigar.front()) ==
                     align::CigarOp::HardClip);
            // Hard-clipped bases are omitted from SEQ.
            test("Chimera supp seq trimmed", supp.seq.size() == 21);
            bool p_sa = false, s_sa = false;
            for (const auto& t : res.primary.tags)
                if (t.first == "SA:Z") p_sa = true;
            for (const auto& t : supp.tags)
                if (t.first == "SA:Z") s_sa = true;
            test("Chimera SA tags", p_sa && s_sa);
        }

        std::remove(fa_path);
    }

    // -x read-type presets set seed/scoring values; unknown names fail.
    {
        Config c = Config::default_mem();
        test("Preset ont2d", c.apply_preset("ont2d") && c.min_seed_len == 14 &&
                             c.scoring.mismatch == -1 && c.scoring.clip_pen == 0);
        c = Config::default_mem();
        test("Preset pacbio", c.apply_preset("pacbio") && c.min_seed_len == 17 &&
                              c.scoring.gap_open == -1);
        c = Config::default_mem();
        test("Preset intractable", c.apply_preset("intractable") &&
                                    c.min_seed_len == 19 &&
                                    c.scoring.mismatch == -4);
        c = Config::default_mem();
        test("Preset unknown", !c.apply_preset("bogus") && c.min_seed_len == 9);
    }

    // Asymmetric chimera: a 25 bp segment on chrY (ref 1) and a 20 bp
    // segment on chrX (ref 0). Score picks chrY as primary by default;
    // -5 promotes the smaller-coordinate chrX hit to primary (soft-clipped,
    // full SEQ) and demotes chrY to hard-clipped supplementary. -Y keeps
    // soft clips (and full SEQ) on the supplementary instead.
    {
        const char* fa_path = "/tmp/bwa_test_chimera5.fa";
        const std::string segA = "AATTACATAACATACACGTCTAGCT";  // 25 bp, chrY
        const std::string segB = "GCTGTGTCCACCCCATCGGA";       // 20 bp, chrX
        {
            std::ofstream fa(fa_path);
            fa << ">chrX\n" << std::string(30, 'G') << segB
               << std::string(30, 'G') << "\n";
            fa << ">chrY\n" << std::string(30, 'C') << segA
               << std::string(30, 'C') << "\n";
        }
        Index idx = Index::build(fa_path);

        io::SeqRecord asym;
        asym.name = "asym";
        asym.seq = segA + segB;
        asym.qual = std::string(45, 'I');

        auto has_hard_clip = [](const AlnRecord& a) {
            if (a.cigar.empty()) return false;
            auto op = align::cigar_op(a.cigar.front());
            if (op == align::CigarOp::HardClip) return true;
            op = align::cigar_op(a.cigar.back());
            return op == align::CigarOp::HardClip;
        };

        // Default: chrY (higher SW score) is primary.
        {
            Aligner aligner(idx);
            AlignmentResult res = aligner.align(asym);
            test("Asym mapped", res.mapped && res.supplementary.size() == 1);
            if (res.mapped && res.supplementary.size() == 1) {
                test("Asym primary", res.primary.rname == "chrY");
                test("Asym supp", res.supplementary[0].rname == "chrX");
                test("Asym supp hard clip",
                     has_hard_clip(res.supplementary[0]));
                test("Asym supp trimmed",
                     res.supplementary[0].seq.size() == 20);
            }
        }

        // -5: chrX (smaller ref id) becomes primary with full SEQ.
        {
            Config cfg = Config::default_mem();
            cfg.smallest_coord_primary = true;
            Aligner aligner(idx, cfg);
            AlignmentResult res = aligner.align(asym);
            test("Minus5 mapped", res.mapped && res.supplementary.size() == 1);
            if (res.mapped && res.supplementary.size() == 1) {
                test("Minus5 primary", res.primary.rname == "chrX");
                test("Minus5 primary full seq", res.primary.seq.size() == 45);
                test("Minus5 supp", res.supplementary[0].rname == "chrY");
                test("Minus5 supp hard clip",
                     has_hard_clip(res.supplementary[0]));
                test("Minus5 supp trimmed",
                     res.supplementary[0].seq.size() == 25);
            }
        }

        // -Y: supplementary keeps soft clips and full SEQ.
        {
            Config cfg = Config::default_mem();
            cfg.soft_clip_supplementary = true;
            Aligner aligner(idx, cfg);
            AlignmentResult res = aligner.align(asym);
            test("SoftY mapped", res.mapped && res.supplementary.size() == 1);
            if (res.mapped && res.supplementary.size() == 1) {
                test("SoftY no hard clip",
                     !has_hard_clip(res.supplementary[0]));
                test("SoftY full seq",
                     res.supplementary[0].seq.size() == 45);
            }
        }

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

    // Long read (>224 bp, i.e. beyond PmrString's inline capacity) must keep
    // its SEQ/QUAL intact. Regression: the aligner used to reset the TLS
    // arena while the reader's heap-backed buffers were still live, so long
    // reads were overwritten by alignment scratch data.
    {
        const char* fa_path = "/tmp/bwa_test_longread.fa";
        const char* idx_prefix = "/tmp/bwa_test_longread_idx";
        // Deterministic pseudo-random reference (LCG), 2000 bp.
        std::string ref;
        ref.reserve(2000);
        uint32_t x = 12345;
        for (int i = 0; i < 2000; ++i) {
            x = x * 1664525u + 1013904223u;
            ref += "ACGT"[(x >> 16) & 3];
        }
        {
            std::ofstream fa(fa_path);
            fa << ">longchr\n" << ref << "\n";
        }
        Index idx = Index::build(fa_path);
        idx.save(idx_prefix);
        Index loaded = Index::load(idx_prefix);
        Aligner aligner(loaded);

        // 800 bp read from offset 500, with 8 substitutions plus a few
        // indels (long-read-like error profile).
        std::string read_seq = ref.substr(500, 800);
        for (int k = 0; k < 8; ++k) {
            size_t p = 40 + static_cast<size_t>(k) * 90;
            read_seq[p] = read_seq[p] == 'A' ? 'C' : 'A';
        }
        read_seq.insert(200, "TT");          // insertion
        read_seq.erase(400, 1);              // deletion
        read_seq.insert(600, "A");           // insertion
        read_seq.erase(700, 1);              // deletion
        io::SeqRecord lr;
        lr.name = "longread";
        lr.seq = read_seq;
        lr.qual = std::string(read_seq.size(), 'I');

        AlignmentResult res = aligner.align(lr);
        test("LongRead mapped", res.mapped);
        test("LongRead rname", res.primary.rname == "longchr");
        test("LongRead pos", res.primary.pos == 501);
        // SEQ/QUAL must be byte-identical to the input (corruption guard).
        test("LongRead seq intact", res.primary.seq == read_seq);
        test("LongRead qual intact", res.primary.qual == lr.qual);
        // The extension must cover almost the whole read (regression: the
        // window/diagonal and scoring bugs clipped most of it).
        int32_t aligned = 0;
        for (uint32_t c : res.primary.cigar) {
            auto op = align::cigar_op(c);
            if (op == align::CigarOp::Match || op == align::CigarOp::Equal ||
                op == align::CigarOp::Diff || op == align::CigarOp::Ins) {
                aligned += align::cigar_len(c);
            }
        }
        test("LongRead span", aligned >= static_cast<int32_t>(read_seq.size()) * 95 / 100);

        std::remove(fa_path);
        std::remove("/tmp/bwa_test_longread_idx.meta");
        std::remove("/tmp/bwa_test_longread_idx.bwt");
        std::remove("/tmp/bwa_test_longread_idx.sa");
        std::remove("/tmp/bwa_test_longread_idx.occ");
        std::remove("/tmp/bwa_test_longread_idx.pac");
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
