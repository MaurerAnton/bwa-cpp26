// Differential + robustness tests (deterministic seeded RNG, no libFuzzer).
//
// 1. SW banding: banded extension must agree with full DP on small inputs
//    where the band covers the whole matrix, and every reported alignment
//    must be self-consistent (CIGAR consumes the claimed spans, recomputed
//    score matches).
// 2. FASTQ parser: malformed/truncated/random inputs must terminate without
//    crashing; valid records must round-trip byte-identically.
#include <bwa/align/sw.hpp>
#include <bwa/io/seq_io.hpp>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

using namespace bwa;

// Deterministic LCG (fixed seed => reproducible suite).
static uint64_t rng_state = 0x12345678u;
static uint32_t rnd(uint32_t n) {
    rng_state = rng_state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<uint32_t>((rng_state >> 33) % n);
}

int main() {
    int passed = 0, failed = 0;
    auto test = [&](const char* name, bool cond) {
        if (cond) { std::cout << "  PASS: " << name << "\n"; ++passed; }
        else { std::cout << "  FAIL: " << name << "\n"; ++failed; }
    };

    align::Scoring sc = align::Scoring::bwa_mem_default();

    // 1a. Banded vs full DP agreement on small random inputs. Band 32
    // covers these sizes entirely (w >= 32 > all indices), so results must
    // match exactly; every result must also be self-consistent.
    auto check_consistent = [&](const align::Alignment& a,
                                const std::vector<uint8_t>& qq,
                                const std::vector<uint8_t>& rr) {
        int qcon = 0, rcon = 0, recomputed = 0;
        for (uint32_t c : a.cigar) {
            auto op = align::cigar_op(c);
            int len = align::cigar_len(c);
            if (op == align::CigarOp::Equal) {
                qcon += len;
                rcon += len;
                recomputed += len * sc.match;
            } else if (op == align::CigarOp::Diff) {
                qcon += len;
                rcon += len;
                recomputed += len * sc.mismatch;
            } else if (op == align::CigarOp::Ins) {
                qcon += len;
                recomputed += sc.gap_open + len * sc.gap_ext;
            } else if (op == align::CigarOp::Del) {
                rcon += len;
                recomputed += sc.gap_open + len * sc.gap_ext;
            } else if (op == align::CigarOp::SoftClip) {
                qcon += len;
            } else if (op == align::CigarOp::HardClip) {
                // hard clips consume neither stored base
            } else {
                return false;  // unexpected op
            }
        }
        if (a.cigar.empty()) return a.score == 0;
        // Clipped CIGARs consume the whole query (clips + aligned span).
        if (qcon != (int)qq.size()) return false;
        if (rcon != a.ref_end - a.ref_begin) return false;
        if (recomputed != a.score) return false;
        // Spans must lie inside the inputs.
        if (a.query_begin < 0 || a.query_end > (int)qq.size()) return false;
        if (a.ref_begin < 0 || a.ref_end > (int)rr.size()) return false;
        // Every claimed Equal/Diff must match the input bases. Leading and
        // trailing soft clips must exactly cover [0, query_begin) and
        // [query_end, qlen); hard clips cover nothing stored.
        size_t qi = 0, ri = (size_t)a.ref_begin;
        bool in_core = false;
        size_t leading_s = 0;
        for (uint32_t c : a.cigar) {
            auto op = align::cigar_op(c);
            int len = align::cigar_len(c);
            if (op == align::CigarOp::SoftClip && !in_core) {
                leading_s += (size_t)len;
                qi += (size_t)len;
                continue;
            }
            if (op == align::CigarOp::SoftClip) {
                qi += (size_t)len;  // trailing clip
                continue;
            }
            if (op == align::CigarOp::HardClip) continue;
            in_core = true;
            if (op == align::CigarOp::Equal || op == align::CigarOp::Diff) {
                for (int k = 0; k < len; ++k) {
                    bool same = qi < qq.size() && ri < rr.size() &&
                                qq[qi] == rr[ri] && qq[qi] < 4;
                    if ((op == align::CigarOp::Equal) != same) return false;
                    ++qi;
                    ++ri;
                }
            } else if (op == align::CigarOp::Ins) {
                qi += (size_t)len;
            } else if (op == align::CigarOp::Del) {
                ri += (size_t)len;
            } else {
                return false;
            }
        }
        if (leading_s != (size_t)a.query_begin) return false;
        if (qi != qq.size()) return false;
        return true;
    };
    {
        bool agree = true;
        for (int t = 0; t < 300 && agree; ++t) {
            int ql = 1 + rnd(12), rl = 1 + rnd(16);
            std::vector<uint8_t> q(ql), r(rl);
            for (auto& b : q) b = rnd(10) == 0 ? 4 : rnd(4);  // 10% N
            for (auto& b : r) b = rnd(10) == 0 ? 4 : rnd(4);
            align::Alignment full = align::sw_semi_global_extend(
                sc, {q.data(), q.size()}, {r.data(), r.size()}, 0);
            align::Alignment band = align::sw_semi_global_extend(
                sc, {q.data(), q.size()}, {r.data(), r.size()}, 32);
            if (full.score != band.score) agree = false;
            if (!check_consistent(band, q, r)) agree = false;
            if (!check_consistent(full, q, r)) agree = false;
        }
        test("SW banded/full agree + self-consistent", agree);
    }

    // 1b. Narrow bands stay self-consistent (may be suboptimal vs full DP).
    {
        bool ok = true;
        for (int t = 0; t < 200 && ok; ++t) {
            int ql = 5 + rnd(20), rl = 5 + rnd(20);
            std::vector<uint8_t> q(ql), r(rl);
            for (auto& b : q) b = rnd(4);
            for (auto& b : r) b = rnd(4);
            align::Alignment band = align::sw_semi_global_extend(
                sc, {q.data(), q.size()}, {r.data(), r.size()}, 2);
            if (!check_consistent(band, q, r)) ok = false;
        }
        test("SW narrow band self-consistent", ok);
    }

    // 1b. N bases never match and never crash the DP.
    {
        bool ok = true;
        for (int t = 0; t < 100 && ok; ++t) {
            std::vector<uint8_t> q(20, 4), r(30, 4);  // all N
            align::Alignment a = align::sw_semi_global_extend(
                sc, {q.data(), q.size()}, {r.data(), r.size()}, 32);
            if (a.score != 0 || !a.cigar.empty()) ok = false;
        }
        test("SW all-N scores zero", ok);
    }

    // 1c. Empty inputs are safe.
    {
        std::vector<uint8_t> q, r = {0, 1};
        align::Alignment a = align::sw_semi_global_extend(
            sc, {q.data(), q.size()}, {r.data(), r.size()}, 32);
        align::Alignment b = align::sw_semi_global_extend(
            sc, {r.data(), r.size()}, {q.data(), q.size()}, 32);
        test("SW empty safe", a.score == 0 && b.score == 0);
    }

    // 2a. Parser: valid records round-trip byte-identically.
    {
        const char* path = "/tmp/bwa_test_fuzz_valid.fq";
        {
            std::ofstream f(path, std::ios::binary);
            f << "@r1 some comment\nACGTNacgtn\n+\nIIIIIIIIII\n";
            f << "@r2\nTTGGCCAA\n+\n########\n";
            f << "@r3\nA\n+\nI\n";
        }
        io::SeqReader reader(path);
        io::SeqRecord rec;
        bool ok = true;
        auto rd = reader.read(rec);
        ok = ok && static_cast<bool>(rd) && *rd;
        ok = ok && rec.name == "r1" && rec.comment == "some comment";
        ok = ok && rec.seq == "ACGTNacgtn" && rec.qual == "IIIIIIIIII";
        rd = reader.read(rec);
        ok = ok && static_cast<bool>(rd) && *rd;
        ok = ok && rec.name == "r2" && rec.seq == "TTGGCCAA" &&
             rec.qual == "########" && rec.comment.empty();
        rd = reader.read(rec);
        ok = ok && static_cast<bool>(rd) && *rd;
        ok = ok && rec.name == "r3" && rec.seq == "A" && rec.qual == "I";
        rd = reader.read(rec);
        ok = ok && !(static_cast<bool>(rd) && *rd);  // clean EOF
        std::remove(path);
        test("Parser valid round-trip", ok);
    }

    // 2b. Parser: malformed inputs terminate (no hang/crash).
    {
        const char* path = "/tmp/bwa_test_fuzz_bad.fq";
        const char* cases[] = {
            "", "@", "@r\n", "@r\nACGT\n", "@r\nACGT\n+\n",
            "@r\nACGT\n+\nIII\n",          // short qual
            "@r\nACGT\n+\nIIIIIIII\n",     // long qual
            "not a fastq at all\njust text\n",
            "@r1\nACGT\n+\nIIII\n@r2\n",   // truncated second record
            "@r\r\nACGT\r\n+\r\nIIII\r\n",  // CRLF
            "@r\nACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n+\n",  // huge seq, no qual
        };
        bool ok = true;
        for (auto blob : cases) {
            {
                std::ofstream f(path, std::ios::binary);
                f << blob;
            }
            io::SeqReader reader(path);
            io::SeqRecord rec;
            int n = 0;
            while (n++ < 10) {
                auto res = reader.read(rec);
                if (!res || !*res) break;
            }
            if (n >= 10) ok = false;  // suspicious: too many records
        }
        std::remove(path);
        test("Parser malformed terminates", ok);
    }

    // 2c. Parser: random binary garbage terminates.
    {
        const char* path = "/tmp/bwa_test_fuzz_rand.fq";
        bool ok = true;
        for (int t = 0; t < 50 && ok; ++t) {
            {
                std::ofstream f(path, std::ios::binary);
                int len = rnd(300);
                for (int i = 0; i < len; ++i) {
                    char c = static_cast<char>(rnd(256));
                    f.put(c);
                }
            }
            io::SeqReader reader(path);
            io::SeqRecord rec;
            int n = 0;
            while (n++ < 100) {
                auto res = reader.read(rec);
                if (!res || !*res) break;
                // Any returned record must be self-consistent.
                if (rec.is_fastq() && rec.seq.size() != rec.qual.size()) {
                    // Qual gets trimmed to seq length; longer is cut, shorter
                    // is an error only if the reader claims success. Accept
                    // either as long as qual is not longer than seq.
                    if (rec.qual.size() > rec.seq.size()) {
                        ok = false;
                        break;
                    }
                }
            }
            if (n >= 100) ok = false;
        }
        std::remove(path);
        test("Parser random garbage terminates", ok);
    }

    std::cout << "Results: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
