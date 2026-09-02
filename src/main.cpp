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

        // Test SA-IS
        {
            memory::Arena arena(64 * 1024);
            index::PackedSequence seq;
            seq.append("ACGTAACCGGTTAA", 14);

            auto sa = index::build_suffix_array_sais(seq, arena);

            // Verify: SA should have 14 entries
            if (sa.size() == 14) {
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

        // Test FASTA reading
        {
            io::SeqReader reader("/tmp/test_debug.fq");
            if (!reader.is_open()) {
                std::cout << "  FASTA read: FAIL (can't open file)\n";
            } else {
                io::SeqRecord rec;
                int count = 0;
                while (true) {
                    auto result = reader.read(rec);
                    if (!result) break;
                    if (!*result) break;
                    count++;
                    std::cout << "  FASTA read: record " << count
                              << " name='" << rec.name.view() << "'"
                              << " seq_len=" << rec.seq.size()
                              << " seq='" << rec.seq.view() << "'"
                              << " is_fasta=" << rec.is_fasta() << "\n";
                }
                std::cout << "  FASTA read: " << count << " records\n";
            }
        }

        std::cout << "All tests passed!\n";
        return 0;
    }

    if (cmd == "index") {
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] << " index <fasta> <prefix>\n";
            return 1;
        }
        std::cout << "Building index from " << argv[2] << " to " << argv[3] << "...\n";
        auto start = std::chrono::high_resolution_clock::now();
        Config cfg = Config::default_mem();
        Index idx = Index::build(argv[2], cfg);
        idx.save(argv[3]);
        auto end = std::chrono::high_resolution_clock::now();
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Index built in " << ms << " ms\n";
        std::cout << "References: " << idx.num_references() << "\n";
        std::cout << "Total length: " << idx.total_length() << "\n";
        return 0;
    }

    if (cmd == "mem") {
        if (argc < 4) {
            std::cerr << "Usage: " << argv[0] << " mem <index> <fastq> [fastq2]\n";
            return 1;
        }
        std::cout << "Loading index " << argv[2] << "...\n";
        auto start = std::chrono::high_resolution_clock::now();
        Pipeline pipe(argv[2]);
        auto load_end = std::chrono::high_resolution_clock::now();
        std::cout << "Index loaded in "
                  << std::chrono::duration_cast<std::chrono::milliseconds>(load_end - start).count()
                  << " ms\n";

        if (argc == 4) {
            std::cout << "Aligning single-end " << argv[3] << "...\n";
            pipe.align_file(argv[3]);
        } else if (argc == 5) {
            std::cout << "Aligning paired-end " << argv[3] << " " << argv[4] << "...\n";
            pipe.align_pair(argv[3], argv[4]);
        }
        return 0;
    }

    std::cerr << "Unknown command: " << cmd << "\n";
    return 1;
}