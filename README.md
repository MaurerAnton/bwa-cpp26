# BWA-CPP26

A modern C++26 implementation of the Burrows-Wheeler Aligner (BWA-MEM).

## Features

- **Modern C++26**: Uses C++26 features for safety and performance
- **Memory Safe**: Custom arena allocators, no manual memory management
- **AGPL-3.0-only Licensed**: Same license as original BWA
- **BWA-MEM Compatible**: Implements the BWA-MEM algorithm for short and long reads
- **Multi-threading**: Parallel alignment with configurable thread count
- **SAM/BAM Output**: Standard alignment output formats with BGZF compression

## Algorithm

- **MEM Finding**: Uses FM-index with backward search (SMEM algorithm)
- **Chaining**: DP-based optimal chaining with gap penalties
- **Alignment**: Banded Smith-Waterman with affine gap costs
- **Suffix Array**: SA-IS linear-time construction (O(n))

## Building

```bash
# Requires: GCC 16.2.1+, zlib
g++ -std=c++26 -Iinclude -O2 -o bwa-cpp26 \
    src/main.cpp src/align/sw.cpp src/index/fm_index.cpp \
    src/pipeline.cpp src/io/seq_io.cpp src/utils/utils.cpp \
    src/core/arena.cpp -lz
```

Or using CMake:
```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

## Usage

### Build Index
```bash
bwa-cpp26 index reference.fa index_prefix
```

### Align Reads (Single-end)
```bash
bwa-cpp26 mem index_prefix reads.fq > alignments.sam
```

### Align Reads (Paired-end)
```bash
bwa-cpp26 mem index_prefix reads_1.fq reads_2.fq > alignments.sam
```

### Output all alignments (secondary hits)
```bash
bwa-cpp26 mem -a index_prefix reads.fq > alignments.sam
```

Flags: `-a` outputs secondary alignments (suppressed by default, like
BWA); `-t <n>` sets the thread count.

### Output to BAM
```bash
bwa-cpp26 mem index_prefix reads.fq | samtools view -b - > alignments.bam
```

Or directly (native BGZF BAM + coordinate-sorted BAI):
```bash
bwa-cpp26 mem index_prefix reads.fq -o alignments.bam
```
The `.bam` extension selects the native writer; any other name is SAM.
BAM output is coordinate-sorted and gets a `<file>.bam.bai` index.

## Configuration

```cpp
bwa::Config cfg = bwa::Config::default_mem();
// Or use presets:
bwa::Config cfg = bwa::Config::fast();
bwa::Config cfg = bwa::Config::accurate();
bwa::Config cfg = bwa::Config::long_reads();  // ONT/PacBio

bwa::Pipeline pipeline("index_prefix", cfg);
pipeline.align_file("reads.fq", "output.sam");
```

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_seed_len` | 19 | Minimum seed length |
| `max_occ` | 500 | Max seed occurrences |
| `band_width` | 32 | DP band width |
| `min_chain_score` | 30 | Minimum chain score |
| `max_gap` | 10000 | Max gap in chaining |
| `num_threads` | 1 | Thread count |

## API Example

```cpp
#include <bwa/pipeline.hpp>
#include <bwa/io/seq_io.hpp>

int main() {
    // Build index
    bwa::Index idx = bwa::Index::build("reference.fa", bwa::Config::default_mem());
    idx.save("index_prefix");
    
    // Load and align
    bwa::Pipeline pipeline("index_prefix");
    
    // Single read
    bwa::io::SeqRecord read;
    read.name = "read1";
    read.seq = "ACGTACGTACGT";
    read.qual = "IIIIIIIIIIII";
    
    auto result = pipeline.aligner().align(read);
    if (result.mapped) {
        std::cout << "Aligned to " << result.primary.rname 
                  << " at " << result.primary.pos << "\n";
    }
    
    // Batch alignment
    bwa::io::BatchSeqReader reader("reads.fq", 1000);
    std::vector<bwa::AlignmentResult> results(1000);
    size_t n = reader.fill();
    pipeline.aligner().align_batch({reader.buffer().data(), n}, results);
}
```

## File Formats

### Index Files
- `prefix.meta` - Metadata (reference names, lengths, MD5)
- `prefix.bwt` - Burrows-Wheeler Transform
- `prefix.sa` - Sampled suffix array
- `prefix.occ` - Occurrence table
- `prefix.pac` - Packed reference sequences

### Output Formats
- **SAM**: Text format with @HD, @SQ, @PG headers
- **BAM**: Binary format with BGZF compression

## Testing

```bash
# Run self-tests
./bwa-cpp26 test

# Run integration test
bash tests/integration/test_pipeline.sh
```

## License

AGPL-3.0-only - Same as original BWA.

## Differences from Original BWA

| Feature | BWA | BWA-CPP26 |
|---------|-----|-----------|
| Language | C | C++26 |
| Memory | Manual | Arena allocators |
| SA Construction | SA-IS (C) | SA-IS (C++26) |
| Threading | pthreads | std::thread |
| BAM | HTStlib | Native BGZF |

## Roadmap

- [x] Memory-mapped index files for large genomes
- [ ] BAI index generation for BAM
- [x] Real data testing (E. coli 4.6 MB: 2000/2000 simulated 150-mers mapped, ~97% MAPQ60)
- [ ] Performance benchmarking
- [x] Supplementary alignment (SA tag)
- [x] Proper paired-end (FR proper-pair flags, TLEN, mate rescue)
- [ ] Base quality recalibration

## Performance notes

- Always benchmark with a Release build (`cmake -B build-release
  -DCMAKE_BUILD_TYPE=Release`): the default unoptimized build is ~3x
  slower at indexing and ~10x slower at alignment.
- Reference points (single thread): E. coli K-12 (4.6 MB) index ~43 s
  Release; 2000 simulated 150-mers with 1% error align in ~100 s
  (~50 ms/read) with 0 unmapped.

## Contributing

See CONTRIBUTING.md for guidelines.
