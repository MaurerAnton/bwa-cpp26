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
- [x] BAI index generation for BAM
- [x] Real data testing (E. coli 4.6 MB: 2000/2000 simulated 150-mers mapped, ~97% MAPQ60)
- [x] Performance benchmarking (see below)
- [x] Supplementary alignment (SA tag)
- [x] Proper paired-end (FR proper-pair flags, TLEN, mate rescue)
- [x] Insert-size distribution estimation (sampled, trimmed moments)
- [ ] Base quality recalibration

## Performance notes

- Always benchmark with a non-sanitized Release build:
  `cmake -B build-release -DCMAKE_BUILD_TYPE=Release -DSANITIZE=OFF`.
  Sanitizers are enabled by default for unoptimized builds and cost
  ~10-15x at runtime.
- Reference points (E. coli K-12, 4.6 MB, 4-core machine): index build
  ~12 s; 2000 simulated 150-mers with 1% error align in ~5.6 s
  (~2.8 ms/read) single-threaded; 0 unmapped, 1978/2000 within 5 bp of
  truth, 1983 MAPQ60.
- 1000 simulated 350 bp-insert pairs: 98.9% flagged proper pair with
  TLEN = +/-350 (insert size is estimated from the first 512 pairs).
- Long reads: a 5 kb read with 5% substitutions + 0.5% indels aligns
  end-to-end (4979/5003 bases, MAPQ 60) in <1 s single-threaded.

## Command-line options

```
bwa mem [-a] [-t N] [-o out.sam|out.bam] [-k N] [-c N] [-w N]
        [-A N] [-B N] [-O N[,N]] [-E N[,N]] [-L N[,N]]
        [-T N] [-R RG] [-M] [-S] [-P] [-p] [-Y] [-5]
        [-x ont2d|pacbio|intractable]
        <index> <fastq> [fastq2] [sam_out]
```

| Flag | Meaning | Default |
|------|---------|---------|
| `-a` | output secondary alignments | off |
| `-t N` | threads | 1 |
| `-o F` | output file (`.bam` → BAM+BAI) | stdout SAM |
| `-k N` | minimum seed length | 9 |
| `-c N` | skip seeds with > N occurrences | 500 |
| `-w N` | band width | 32 |
| `-A/-B` | match/mismatch score | 1/-4 |
| `-O/-E` | gap open/extension penalty | -6/-1 |
| `-L` | clipping penalty | -5 |
| `-T N` | minimum output score | 0 (off) |
| `-R STR` | read group (`@RG\tID:..\tSM:..`) | none |
| `-M` | mark split hits as secondary | off |
| `-S` | skip mate rescue | off |
| `-P` | skip pairing | off |
| `-p` | interleaved pairs in one FASTQ | off |
| `-Y` | soft-clip (not hard-clip) split hits | off |
| `-5` | smallest-coordinate split hit is primary | off |
| `-x PRESET` | read-type preset (`ont2d`, `pacbio`, `intractable`) | short-read defaults |

## Contributing

See CONTRIBUTING.md for guidelines.
