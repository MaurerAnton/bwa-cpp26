#!/bin/bash
# Integration test: build index and align reads end to end.
# Usage: [BWA_CPP26=path/to/bwa] bash tests/integration/test_pipeline.sh
set -eu

BWA_CPP26="${BWA_CPP26:-./build/src/bwa}"
if [ ! -x "$BWA_CPP26" ]; then
    echo "FAIL: bwa binary not found at $BWA_CPP26 (set BWA_CPP26=...)"
    exit 1
fi

TEST_DIR=$(mktemp -d)
trap 'rm -rf "$TEST_DIR"' EXIT
echo "Test directory: $TEST_DIR"

# Small two-region reference: ACGT repeat (1-128) + AACCGGTT repeat (129-256)
cat > "$TEST_DIR/ref.fa" << 'EOF'
>test_chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
AACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTTAACCGGTT
EOF

# Two mappable reads (one per region) + one absent read (no SAM record expected)
cat > "$TEST_DIR/reads.fq" << 'EOF'
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

fail() { echo "FAIL: $1"; exit 1; }

echo "=== Building index ==="
"$BWA_CPP26" index "$TEST_DIR/ref.fa" "$TEST_DIR/test_idx" > /dev/null

for ext in meta bwt sa occ pac; do
    [ -f "$TEST_DIR/test_idx.$ext" ] || fail "missing index file test_idx.$ext"
done
echo "PASS: all index files created"

echo "=== Aligning reads ==="
"$BWA_CPP26" mem "$TEST_DIR/test_idx" "$TEST_DIR/reads.fq" "$TEST_DIR/aln.sam" \
    > /dev/null 2> "$TEST_DIR/stderr.txt"
[ -s "$TEST_DIR/aln.sam" ] || { cat "$TEST_DIR/stderr.txt"; fail "empty SAM output"; }

grep -q '^@HD' "$TEST_DIR/aln.sam" || fail "missing @HD header"
grep -q '^@SQ.*SN:test_chr1.*LN:256' "$TEST_DIR/aln.sam" || fail "missing @SQ header"

# read1/read2: mapped to test_chr1, MAPQ 60, perfect 20= CIGAR
check_mapped() {
    local line
    line=$(awk -v q="$1" '$1 == q' "$TEST_DIR/aln.sam")
    [ -n "$line" ] || fail "$1 has no SAM record"
    echo "$line" | awk -v q="$1" -v pos="$2" \
        '$3 == "test_chr1" && $4 == pos && $5 == 60 && $6 == "20=" {
             exit 0
         }
         { print "unexpected record for " q ": " $0 > "/dev/stderr"; exit 1 }' \
        || fail "$1 mapping incorrect"
}
check_mapped read1 45
echo "PASS: read1 mapped correctly"
check_mapped read2 169
echo "PASS: read2 mapped correctly"

# read3 has no match in the reference: must not appear as mapped
if awk '$1 == "read3" && $3 != "*" && and($2, 4) == 0' "$TEST_DIR/aln.sam" | grep -q .; then
    fail "read3 unexpectedly mapped"
fi
echo "PASS: absent read not mapped"

echo "=== Integration test completed ==="
