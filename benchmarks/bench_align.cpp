#include <bwa/align/sw.hpp>
#include <bwa/index/fm_index.hpp>
#include <bwa/core/arena.hpp>
#include <chrono>
#include <iostream>
#include <array>
#include <random>

using namespace bwa;

int main() {
    // SW benchmark: 1000 x 100bp global alignments
    {
        align::Scoring sc = align::Scoring::bwa_mem_default();
        std::array<uint8_t, 100> q{}, r{};
        std::mt19937 rng(42);
        for (auto& b : q) b = rng() % 4;
        for (auto& b : r) b = rng() % 4;
        auto t0 = std::chrono::high_resolution_clock::now();
        int total = 0;
        for (int i = 0; i < 1000; ++i)
            total += align::sw_global(sc, q, r, 0).score;
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "sw_global 1000x100bp: " << ms << " ms (checksum=" << total << ")\n";
    }
    // FM-index build benchmark
    {
        memory::Arena arena(64 * 1024 * 1024);
        index::PackedSequence seq;
        std::string s(10000, 'A');
        for (size_t i = 0; i < s.size(); ++i) s[i] = "ACGT"[i % 4];
        seq.append(s.c_str(), s.size());
        auto t0 = std::chrono::high_resolution_clock::now();
        index::FMIndex idx = index::FMIndex::build(seq, arena);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "fm_build 10kb: " << ms << " ms size=" << idx.size() << "\n";
    }
    return 0;
}
