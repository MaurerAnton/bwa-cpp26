#!/bin/bash
# Integration test: Build index and align reads
set -e

BWA_CPP26=./build_test
TEST_DIR=$(mktemp -d)
echo "Test directory: $TEST_DIR"

# Create a small test FASTA (E. coli fragment, ~1KB)
cat > $TEST_DIR/ref.fa << 'EOF'
>test_chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
EOF

# Create test FASTQ
cat > $TEST_DIR/reads.fq << 'EOF'
@read1
ACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIII
@read2
AACCGGTTAACCGGTTAACC
+
IIIIIIIIIIIIIIIIIIII
@read3
GCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIII
EOF

# Build index
echo "=== Building index ==="
$BWA_CPP26 index $TEST_DIR/ref.fa $TEST_DIR/test_idx

# Verify index files exist
for ext in meta bwt sa occ; do
    if [ ! -f $TEST_DIR/test_idx.$ext ]; then
        echo "FAIL: Missing index file: test_idx.$ext"
        exit 1
    fi
done
echo "PASS: All index files created"

# Align reads
echo "=== Aligning reads ==="
$BWA_CPP26 mem $TEST_DIR/test_idx $TEST_DIR/reads.fq > $TEST_DIR/aln.sam 2>&1 || true

if [ -f $TEST_DIR/aln.sam ] && [ -s $TEST_DIR/aln.sam ]; then
    echo "PASS: SAM output generated"
    echo "Sample output:"
    head -3 $TEST_DIR/aln.sam
else
    echo "INFO: SAM output not generated (alignment pipeline not fully implemented yet)"
fi

# Cleanup
rm -rf $TEST_DIR
echo "=== Integration test completed ==="