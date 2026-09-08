#include <bwa/core/vector.hpp>
#include <bwa/core/sort.hpp>
#include <bwa/core/arena.hpp>
#include <chrono>
#include <iostream>
#include <random>

using namespace bwa;

int main() {
    // radix_sort benchmark
    {
        core::Vector<uint32_t> v;
        v.reserve(1'000'000);
        std::mt19937 rng(42);
        for (int i = 0; i < 1'000'000; ++i) v.push_back(rng());
        auto t0 = std::chrono::high_resolution_clock::now();
        core::radix_sort(v);
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "radix_sort 1M u32: " << ms << " ms sorted="
                  << core::is_sorted(v.begin(), v.end()) << "\n";
    }
    // arena benchmark
    {
        memory::Arena arena(64 * 1024 * 1024);
        auto t0 = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100000; ++i) {
            volatile void* p = arena.allocate(64, 8);
            (void)p;
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "arena 100k x 64B: " << ms << " ms\n";
    }
    return 0;
}
