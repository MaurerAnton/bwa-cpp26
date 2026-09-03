#!/bin/bash
# Test with real sequencing data
# Downloads a small dataset and tests alignment

set -e

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
WORK_DIR="$TEST_DIR/work"
mkdir -p "$WORK_DIR"

echo "=== Real Data Test ==="
echo "Working directory: $WORK_DIR"

# Download a small E. coli reference and reads
cd "$WORK_DIR"

# E. coli K-12 reference (small, ~4.6Mbp)
if [ ! -f ecoli_ref.fa ]; then
    echo "Downloading E. coli reference..."
    curl -sL -o ecoli_ref.fa.gz "https://ftp.ensemblgenomes.org/pub/bacteria/release-57/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz" 2>/dev/null || \
    curl -sL -o ecoli_ref.fa.gz "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_000913.3&rettype=fasta&retmode=text" 2>/dev/null || \
    wget -q -O ecoli_ref.fa.gz "https://ftp.ensemblgenomes.org/pub/bacteria/release-57/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.toplevel.fa.gz" || \
    echo "Warning: Could not download E. coli reference"
fi

if [ -f ecoli_ref.fa.gz ]; then
    gunzip -kf ecoli_ref.fa.gz 2>/dev/null || true
fi

if [ ! -f ecoli_ref.fa ]; then
    # Create a synthetic test as fallback
    echo "Creating synthetic test data..."
    cat > ecoli_ref.fa << 'EOF'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
fi

# Create synthetic reads with known positions
cat > ecoli_reads.fq << 'EOF'
@read1
ACGTACGTAC
+
IIIIIIIIII
@read2
CGTACGTACG
+
IIIIIIIIII
@read3
GTACGTACGT
+
IIIIIIIIII
@read4
ACGTACGTACGTACGTACG
+
IIIIIIIIIIIIIIIIII
@read5
TACGTACGTACGTACGTAC
+
IIIIIIIIIIIIIIIIII
EOF

echo ""
echo "=== Building Index ==="
"../../build_test" index ecoli_ref.fa ecoli_idx 2>&1

echo ""
echo "=== Aligning Reads ==="
"../../build_test" mem ecoli_idx ecoli_reads.fq 2>&1 | head -20

echo ""
echo "=== Test Complete ==="
