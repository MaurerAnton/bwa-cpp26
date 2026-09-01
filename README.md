# BWA-CPP26: Modern C++26 Rewrite of BWA

A memory-safe, high-performance rewrite of BWA (Burrows-Wheeler Aligner) in C++26.

## Features

- **C++26 Standard**: Uses `std::simd`, `std::execution::par`, `std::expected`, `std::mdspan`, `std::format`
- **Memory Safety**: PMR allocators, arena allocation, bounds-checked containers, no raw pointers in hot paths
- **Sanitizer-Friendly**: Built with ASan/UBSan/MSan support
- **Performance**: SIMD-accelerated kernels, cache-friendly data structures, parallel algorithms
- **Modular Design**: Clean separation of index, alignment, I/O, and core utilities

## Architecture

```
bwa-cpp26/
├── include/bwa/
│   ├── core/           # Memory, containers, strings, sorting, hashing
│   ├── index/          # FM-index, BWT, suffix arrays
│   ├── align/          # MEM finding, Smith-Waterman, chaining
│   ├── io/             # FASTQ/FASTA/BAM parsing, streaming
│   └── pipeline.hpp    # High-level API
├── src/
│   ├── core/           # Core implementations
│   ├── index/          # Index building/loading
│   ├── align/          # Alignment kernels
│   ├── io/             # I/O implementations
│   ├── utils/          # Utilities
│   └── main.cpp        # CLI entry point
├── tests/
│   ├── unit/           # Unit tests
│   ├── integration/    # Integration tests
│   └── fuzz/           # Fuzz targets
└── benchmarks/         # Performance benchmarks
```

## Core Replacements

| C Component | C++26 Replacement |
|-------------|-------------------|
| `kvec.h` | `core::Vector<T>` (PMR, SSO, bounds-checked) |
| `kstring.h` | `core::PmrString` (SSO, `to_chars`, `format_to`) |
| `khash.h` | `core::HashMap<K,V>` (Robin Hood, heterogeneous lookup) |
| `kseq.h` | `io::SeqReader` (streaming, gzip, mmap) |
| `ksort.h` | `core::radix_sort`, `core::parallel_sort` |
| `kthread.h` | `std::jthread`, `std::execution::par` |
| `bwt.c` | `index::FMIndex` (SA-IS, sampled SA, occ table) |
| `ksw.c` | `align::sw_*` (SIMD banded DP) |

## Building

### Requirements
- C++26 compiler (GCC 13+, Clang 16+, MSVC 19.40+)
- CMake 3.28+
- zlib (for gzip support)
- Optional: [Highway](https://github.com/google/highway) for portable SIMD
- Optional: [Abseil](https://abseil.io) for `flat_hash_map`

### Build Commands

```bash
# Configure
cmake -B build -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF

# Build
cmake --build build -j$(nproc)

# Run tests
cd build && ctest --output-on-failure

# Run benchmarks
./build/bwa-bench

# Run self-test
./build/bwa test
```

### Sanitizer Build (for development)
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DSANITIZE=ON
cmake --build build -j$(nproc)
```

## Usage

```bash
# Build index
./bwa index reference.fasta reference.idx

# Single-end alignment
./bwa mem reference.idx reads.fq > aln.sam

# Paired-end alignment
./bwa mem reference.idx reads_1.fq reads_2.fq > aln.sam
```

## Memory Model

- **Thread-local arenas**: Each thread gets a 64KB arena for allocations
- **PMR polymorphic allocators**: All containers use `std::pmr::polymorphic_allocator`
- **No global state**: Fully thread-safe, no static initialization order issues
- **Deterministic cleanup**: RAII everywhere, no memory leaks possible

## Testing

```bash
# Unit tests
./build/bwa-tests

# Fuzzing (with libFuzzer)
clang++ -fsanitize=fuzzer -std=c++26 -Iinclude tests/fuzz/fuzz_parser.cpp -o fuzz_parser
./fuzz_parser corpus/
```

## Performance Targets

- **Index building**: < 2x original BWA (SA-IS vs prefix-doubling)
- **Alignment throughput**: > 1.5x original BWA-MEM (SIMD + better memory layout)
- **Memory usage**: < 1.2x original (packed structures, no fragmentation)
- **Latency**: P99 < 10ms per read (vs ~50ms in original)

## License

AGPL-3.0-only

## Contributing

1. All code must pass sanitizers (ASan/UBSan/MSan)
2. All new code must have unit tests
3. Performance regressions > 5% require justification
4. Follow C++ Core Guidelines and C++26 best practices