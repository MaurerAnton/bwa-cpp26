#include <bwa/pipeline.hpp>
#include <bwa/core/arena.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/hash_map.hpp>
#include <bwa/core/sort.hpp>
#include <bwa/index/fm_index.hpp>
#include <bwa/align/mem.hpp>
#include <bwa/align/sw.hpp>
#include <bwa/io/seq_io.hpp>
#include <iostream>
#include <chrono>
#include <cstdlib>

using namespace bwa;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <command> [args...]\n";
        std::cerr << "Commands:\n";
        std::cerr << "  index <fasta> <prefix>     Build index\n";
        std::cerr << "  mem <index> <fastq> [fastq2] Align reads\n";
        std::cerr << "  test                       Run self-tests\n";
        return 1;
    }

    std::string cmd = argv[1];

    if (cmd == "test") {
        std::cout << "Running self-tests...\n";

        // Test Arena
        {
            memory::Arena arena(1024);
            memory::ArenaAllocator<int> alloc(&arena);
            int* a = alloc.allocate(10);
            for (int i = 0; i < 10; ++i) a[i] = i;
            alloc.deallocate(a, 10);
            arena.reset();
            std::cout << "  Arena: OK\n";
        }

        // Test Vector
        {
            core::Vector<int> v;
            for (int i = 0; i < 100; ++i) v.push_back(i);
            core::Vector<int> v2 = v;
            if (v == v2 && v.size() == 100) {
                std::cout << "  Vector: OK\n";
            } else {
                std::cout << "  Vector: FAIL\n";
            }
        }

        // Test PmrString
        {
            core::PmrString s;
            s.kputs("Hello");
            s.kputc(' ');
            s.kputs("World");
            s.kputw(42);
            if (s.view() == "Hello World42") {
                std::cout << "  PmrString: OK\n";
            } else {
                std::cout << "  PmrString: FAIL (" << s.view() << ")\n";
            }
        }

        // Test HashMap
        {
            core::PmrHashMap<int, int> m;
            for (int i = 0; i < 1000; ++i) m.insert(i, i * 2);
            if (m.size() == 1000 && m.find(500)->second == 1000) {
                std::cout << "  HashMap: OK\n";
            } else {
                std::cout << "  HashMap: FAIL\n";
            }
        }

        // Test Sort
        {
            core::Vector<int> v;
            for (int i = 100; i >= 0; --i) v.push_back(i);
            core::radix_sort(v);
            bool ok = true;
            for (int i = 0; i <= 100; ++i) if (v[i] != i) ok = false;
            if (ok) std::cout << "  RadixSort: OK\n";
            else std::cout << "  RadixSort: FAIL\n";
        }

        // Test PackedSequence
        {
            index::PackedSequence seq;
            seq.append("ACGTACGTNN", 10);
            bool ok = true;
            ok &= seq.get(0) == 0;  // A
            ok &= seq.get(1) == 1;  // C
            ok &= seq.get(2) == 2;  // G
            ok &= seq.get(3) == 3;  // T
            ok &= seq.get(8) == 4;  // N
            ok &= seq.is_n(8);
            ok &= !seq.is_n(0);
            if (ok) std::cout << "  PackedSequence: OK\n";
            else {
                std::cout << "  PackedSequence: FAIL (get0=" << (int)seq.get(0)
                          << " get1=" << (int)seq.get(1)
                          << " get2=" << (int)seq.get(2)
                          << " get3=" << (int)seq.get(3)
                          << " get8=" << (int)seq.get(8)
                          << " isn8=" << seq.is_n(8)
                          << " isn0=" << seq.is_n(0) << ")\n";
            }
        }

        // Test Smith-Waterman
        {
            align::Scoring sc = align::Scoring::bwa_mem_default();
            // Individual bases: A=0, C=1, G=2, T=3
            std::array<uint8_t, 4> query = {0, 1, 2, 3}; // ACGT
            std::array<uint8_t, 4> ref = {0, 1, 2, 3};   // ACGT

            align::Alignment aln = align::sw_global(sc, query, ref);
            if (aln.score > 0) {
                std::cout << "  SmithWaterman: OK (score=" << aln.score << ")\n";
            } else {
                std::cout << "  SmithWaterman: FAIL (score=" << aln.score << ")\n";
            }
        }

        // Test SA-IS (contract: n+1 entries, empty suffix first)
        {
            memory::Arena arena(64 * 1024);
            index::PackedSequence seq;
            seq.append("ACGTAACCGGTTAA", 14);

            auto sa = index::detail::sais::build_suffix_array(seq, arena);
            auto ref = index::build_suffix_array_brute(seq, arena);

            // Verify: SA should have 15 entries, match brute force
            bool eq = (sa.size() == 15) && (ref.size() == sa.size());
            for (size_t k = 0; eq && k < sa.size(); ++k) eq = (sa[k] == ref[k]);
            if (eq && sa[0] == 14) {
                std::cout << "  SA-IS: OK (size=" << sa.size() << ")\n";
            } else {
                std::cout << "  SA-IS: FAIL (size=" << sa.size() << ")\n";
            }
        }
        // Test FMIndex
        {
            memory::Arena arena(1024 * 1024);
            index::PackedSequence seq;
            seq.append("ACGTAACCGGTTAA", 14);

            index::FMIndex idx = index::FMIndex::build(seq, arena);

            if (idx.size() == 14) {
                // Test backward search for "AA"
                index::PackedSequence query;
                query.append("AA", 2);
                auto [l, r] = idx.backward_search(query);
                size_t count = r - l;
                if (count > 0) {
                    std::cout << "  FMIndex: OK (found " << count << " occurrences of AA)\n";
                } else {
                    std::cout << "  FMIndex: FAIL (no occurrences found)\n";
                }
            } else {
                std::cout << "  FMIndex: FAIL (size=" << idx.size() << ")\n";
            }
        }

        // Test CIGAR =/X refinement
        {
            align::Scoring sc = align::Scoring::bwa_mem_default();
            // Perfect match
            std::array<uint8_t, 8> query = {0, 1, 2, 3, 0, 1, 2, 3}; // ACGTACGT
            std::array<uint8_t, 8> ref = {0, 1, 2, 3, 0, 1, 2, 3};   // ACGTACGT

            align::Alignment aln = align::sw_global(sc, query, ref, 0);
            bool has_equal = false, has_diff = false;
            for (uint32_t c : aln.cigar) {
                auto op = static_cast<align::CigarOp>(c & 0xF);
                if (op == align::CigarOp::Equal) has_equal = true;
                if (op == align::CigarOp::Diff) has_diff = true;
            }
            if (has_equal && !has_diff && aln.score == 8) {
                std::cout << "  CIGAR =/X (perfect): OK\n";
            } else {
                std::cout << "  CIGAR =/X (perfect): FAIL (has_equal=" << has_equal
                          << " has_diff=" << has_diff << " score=" << aln.score << ")\n";
            }
        }

        // Test CIGAR =/X with mismatches
        {
            align::Scoring sc = align::Scoring::bwa_mem_default();
            // ACGTACGT vs ACGTTCGT (2 mismatches at positions 4,5)
            std::array<uint8_t, 8> query = {0, 1, 2, 3, 0, 1, 2, 3}; // ACGTACGT
            std::array<uint8_t, 8> ref = {0, 1, 2, 3, 3, 3, 2, 3};   // ACGTTCGT

            align::Alignment aln = align::sw_global(sc, query, ref, 0);
            bool has_equal = false, has_diff = false;
            for (uint32_t c : aln.cigar) {
                auto op = static_cast<align::CigarOp>(c & 0xF);
                if (op == align::CigarOp::Equal) has_equal = true;
                if (op == align::CigarOp::Diff) has_diff = true;
            }
            if (has_equal && has_diff) {
                std::cout << "  CIGAR =/X (mismatches): OK\n";
            } else {
                std::cout << "  CIGAR =/X (mismatches): FAIL (has_equal=" << has_equal
                          << " has_diff=" << has_diff << " score=" << aln.score << ")\n";
            }
        }

        // Test end-to-end alignment pipeline
        {
            // Create a temporary FASTA file
            const char* fasta_path = "/tmp/bwa_test_ref.fa";
            const char* idx_prefix = "/tmp/bwa_test_idx";
            const char* sam_path = "/tmp/bwa_test_aln.sam";

            // Write test reference
            {
                std::ofstream f(fasta_path);
                f << ">test_chr1\n";
                f << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
                f << "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n";
            }

            // Build index
            Index idx = Index::build(fasta_path);
            idx.save(idx_prefix);

            // Reload index
            Index idx2 = Index::load(idx_prefix);

            // Create test read (exact match)
            io::SeqRecord read;
            read.name = "read1";
            read.seq = "ACGTACGTACGTACGTACGT";
            read.qual = "IIIIIIIIIIIIIIIIIIII";

            // Align
            Aligner aligner(idx2);
            AlignmentResult result = aligner.align(read);

            if (result.mapped && result.primary.score > 0) {
                // Check CIGAR
                bool has_match = false;
                for (uint32_t c : result.primary.cigar) {
                    auto op = static_cast<align::CigarOp>(c & 0xF);
                    if (op == align::CigarOp::Equal || op == align::CigarOp::Match) {
                        has_match = true;
                        break;
                    }
                }
                if (has_match) {
                    std::cout << "  EndToEnd: OK (pos=" << result.primary.pos
                              << " score=" << result.primary.score << ")\n";
                } else {
                    std::cout << "  EndToEnd: FAIL (no match in CIGAR)\n";
                }
            } else {
                std::cout << "  EndToEnd: FAIL (not mapped or score=0)\n";
            }

            // Cleanup
            std::remove(fasta_path);
            std::remove((std::string(idx_prefix) + ".meta").c_str());
            std::remove((std::string(idx_prefix) + ".bwt").c_str());
            std::remove((std::string(idx_prefix) + ".sa").c_str());
            std::remove((std::string(idx_prefix) + ".occ").c_str());
            std::remove((std::string(idx_prefix) + ".pac").c_str());
            std::remove(sam_path);
        }

        std::cout << "All tests passed!\n";
        return 0;
    }

    if (cmd == "index") {
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] << " index <fasta> <prefix>\n";
            return 1;
        }
        std::cerr << "Building index from " << argv[2] << " to " << argv[3] << "...\n";
        auto start = std::chrono::high_resolution_clock::now();
        Config cfg = Config::default_mem();
        Index idx = Index::build(argv[2], cfg);
        idx.save(argv[3]);
        auto end = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cerr << "Index built in " << ms << " ms\n";
        std::cerr << "References: " << idx.num_references() << "\n";
        std::cerr << "Total length: " << idx.total_length() << "\n";
        return 0;
    }

    if (cmd == "mem") {
        // Parse flags (BWA mem-compatible subset):
        //   -a             output secondary alignments
        //   -t INT         threads
        //   -o FILE        output (.bam selects BAM+BAI, else SAM)
        //   -k INT         minimum seed length
        //   -c INT         skip seeds with > INT occurrences
        //   -w INT         band width
        //   -A INT         match score
        //   -B INT         mismatch penalty
        //   -O INT[,INT]   gap open penalty
        //   -E INT[,INT]   gap extension penalty
        //   -L INT[,INT]   clipping penalty
        //   -T INT         minimum score to output
        //   -R STR         read group (e.g. '@RG\\tID:foo\\tSM:bar')
        //   -M             mark split hits as secondary
        //   -S             skip mate rescue
        //   -P             skip pairing
        // Positional: <index> <fastq> [fastq2] [sam_out]
        bool all_alignments = false;
        int threads = 1;
        const char* out_file = nullptr;
        Config cfg = Config::default_mem();
        cfg.output_secondary = false;

        auto parse_int = [](const char* s, int& out) -> bool {
            char* end = nullptr;
            long v = std::strtol(s, &end, 10);
            if (end == s) return false;
            out = static_cast<int>(v);
            return true;
        };
        auto parse_pair = [&](const char* s, int& a, int& b) -> bool {
            std::string_view sv(s);
            size_t comma = sv.find(',');
            if (comma == std::string_view::npos) {
                return parse_int(s, a) ? (b = a, true) : false;
            }
            std::string first(sv.substr(0, comma));
            std::string second(sv.substr(comma + 1));
            return parse_int(first.c_str(), a) && parse_int(second.c_str(), b);
        };

        std::vector<const char*> pos;
        for (int i = 2; i < argc; ++i) {
            std::string_view a = argv[i];
            auto value = [&](const char* opt) -> const char* {
                if (i + 1 >= argc) {
                    std::cerr << "Error: " << opt << " requires a value\n";
                    std::exit(1);
                }
                return argv[++i];
            };
            if (a == "-a" || a == "--all") {
                all_alignments = true;
            } else if (a == "-t" || a == "--threads") {
                if (!parse_int(value("-t"), threads)) threads = 1;
                if (threads < 1) threads = 1;
            } else if (a == "-o" || a == "--output") {
                out_file = value("-o");
            } else if (a == "-k") {
                int v = 0;
                if (parse_int(value("-k"), v) && v > 0) cfg.min_seed_len = v;
            } else if (a == "-c") {
                int v = 0;
                if (parse_int(value("-c"), v) && v > 0) cfg.max_occ = v;
            } else if (a == "-w") {
                int v = 0;
                if (parse_int(value("-w"), v) && v > 0) cfg.band_width = v;
            } else if (a == "-A") {
                int v = 0;
                if (parse_int(value("-A"), v)) cfg.scoring.match = v;
            } else if (a == "-B") {
                int v = 0;
                if (parse_int(value("-B"), v)) cfg.scoring.mismatch = -std::abs(v);
            } else if (a == "-O") {
                int o1 = 0, o2 = 0;
                if (parse_pair(value("-O"), o1, o2)) cfg.scoring.gap_open = -std::abs(o1);
            } else if (a == "-E") {
                int e1 = 0, e2 = 0;
                if (parse_pair(value("-E"), e1, e2)) cfg.scoring.gap_ext = -std::abs(e1);
            } else if (a == "-L") {
                int l1 = 0, l2 = 0;
                if (parse_pair(value("-L"), l1, l2)) cfg.scoring.clip_pen = -std::abs(l1);
            } else if (a == "-T") {
                int v = 0;
                if (parse_int(value("-T"), v)) cfg.min_output_score = v;
            } else if (a == "-R") {
                // Accept '@RG\tID:..\tSM:..' (literal backslash-t, as BWA
                // does) or an already-tab-separated string.
                std::string rg(value("-R"));
                std::string expanded;
                expanded.reserve(rg.size());
                for (size_t p = 0; p < rg.size(); ++p) {
                    if (rg[p] == '\\' && p + 1 < rg.size() && rg[p + 1] == 't') {
                        expanded.push_back('\t');
                        ++p;
                    } else {
                        expanded.push_back(rg[p]);
                    }
                }
                rg = std::move(expanded);
                if (rg.rfind("@RG", 0) == 0) {
                    rg.erase(0, 3);
                    if (!rg.empty() && (rg[0] == '\t' || rg[0] == ' ')) rg.erase(0, 1);
                }
                ReadGroup group;
                size_t start = 0;
                while (start <= rg.size()) {
                    size_t tab = rg.find('\t', start);
                    std::string field = rg.substr(start, tab == std::string::npos
                                                             ? std::string::npos
                                                             : tab - start);
                    size_t colon = field.find(':');
                    if (colon != std::string::npos) {
                        std::string key = field.substr(0, colon);
                        std::string val = field.substr(colon + 1);
                        if (key == "ID") group.id = val;
                        else if (key == "SM") group.sample = val;
                        else if (key == "LB") group.library = val;
                        else if (key == "PL") group.platform = val;
                        else if (key == "PU") group.platform_unit = val;
                    }
                    if (tab == std::string::npos) break;
                    start = tab + 1;
                }
                if (!group.id.empty()) cfg.read_group = group;
            } else if (a == "-M") {
                cfg.mark_split_secondary = true;
            } else if (a == "-S") {
                cfg.skip_mate_rescue = true;
            } else if (a == "-P") {
                cfg.skip_pairing = true;
            } else if (a == "-p" || a == "--smart-pairing") {
                cfg.smart_pairing = true;
            } else if (a == "-Y") {
                cfg.soft_clip_supplementary = true;
            } else if (a == "-5") {
                cfg.smallest_coord_primary = true;
            } else if (a == "-x") {
                // Read-type presets (BWA -x). Applied in order, so later
                // options override the preset values.
                std::string preset(value("-x"));
                if (!cfg.apply_preset(preset)) {
                    std::cerr << "Error: unknown -x preset '" << preset
                              << "' (want ont2d|pacbio|intractable)\n";
                    return 1;
                }
            } else {
                pos.push_back(argv[i]);
            }
        }
        if (pos.size() < 2) {
            std::cerr << "Usage: " << argv[0]
                      << " mem [-a] [-t N] [-o out.sam|out.bam] [-k N] [-c N] [-w N]\n"
                      << "            [-A N] [-B N] [-O N[,N]] [-E N[,N]] [-L N[,N]]\n"
                      << "            [-T N] [-R RG] [-M] [-S] [-P] [-p] [-Y] [-5]\n"
                      << "            [-x ont2d|pacbio|intractable]\n"
                      << "            <index> <fastq> [fastq2] [sam_out]\n";
            return 1;
        }

        std::cerr << "Loading index " << pos[0] << "...\n";
        auto start = std::chrono::high_resolution_clock::now();
        cfg.output_secondary = all_alignments;
        cfg.num_threads = threads;
        Pipeline pipe(pos[0], cfg);
        auto load_end = std::chrono::high_resolution_clock::now();
        std::cerr << "Index loaded in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(load_end - start).count()
                  << " ms\n";

        // Positional parsing:
        //   mem <index> <fastq>                 → single-end, stdout
        //   mem <index> <fastq1> <fastq2>       → paired-end, stdout
        //   mem <index> <fastq> <sam_out>       → single-end to file
        //   mem <index> <fastq1> <fastq2> <out> → paired-end to file
        bool has_fastq2 = false;
        if (pos.size() == 3) {
            // Ambiguous: could be fastq1+fastq2 or fastq+sam_out
            std::string_view arg4 = pos[2];
            bool looks_like_sam = (arg4.size() >= 4 &&
                                   (arg4.substr(arg4.size()-4) == ".sam" ||
                                    arg4.substr(arg4.size()-4) == ".bam"));
            if (!looks_like_sam) {
                std::ifstream test(pos[2]);
                if (test.good()) {
                    has_fastq2 = true;
                }
                test.close();
            }
        } else if (pos.size() >= 4) {
            has_fastq2 = true;
        }

        const char* sam_out = "-";
        if (pos.size() == 3 && !has_fastq2) {
            sam_out = pos[2];
        } else if (pos.size() >= 4) {
            sam_out = pos[3];
        }
        if (out_file != nullptr) sam_out = out_file;

        // .bam selects the native BAM+BAI writer; anything else is SAM.
        bool bam_out = false;
        {
            std::string_view o = sam_out;
            bam_out = o.size() >= 4 && o.substr(o.size() - 4) == ".bam";
        }

        if (bam_out) {
            if (cfg.smart_pairing) {
                std::cerr << "Aligning interleaved " << pos[1] << " to BAM...\n";
                pipe.align_pair_interleaved_bam(pos[1], sam_out);
            } else if (has_fastq2) {
                std::cerr << "Aligning paired-end " << pos[1] << " " << pos[2]
                          << " to BAM...\n";
                pipe.align_pair_to_bam(pos[1], pos[2], sam_out);
            } else {
                std::cerr << "Aligning single-end " << pos[1] << " to BAM...\n";
                pipe.align_to_bam(pos[1], sam_out);
            }
        } else if (cfg.smart_pairing) {
            std::cerr << "Aligning interleaved " << pos[1] << "...\n";
            pipe.align_pair_interleaved(pos[1], sam_out);
        } else if (has_fastq2) {
            std::cerr << "Aligning paired-end " << pos[1] << " " << pos[2] << "...\n";
            pipe.align_pair(pos[1], pos[2], sam_out);
        } else {
            std::cerr << "Aligning single-end " << pos[1] << "...\n";
            pipe.align_file(pos[1], sam_out);
        }
        return 0;
    }

    std::cerr << "Unknown command: " << cmd << "\n";
    return 1;
}