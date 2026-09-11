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

# read1/read2 come from perfect tandem repeats, so every placement is
# equally valid and MAPQ must be low (ambiguous). Assert a valid placement
# inside the repeat block with an intact 20= CIGAR/SEQ, not an exact locus.
check_repeat_mapped() {
    local line
    line=$(awk -v q="$1" '$1 == q && $2 != 256 && $2 != 2048' "$TEST_DIR/aln.sam" | head -n 1)
    [ -n "$line" ] || fail "$1 has no primary SAM record"
    echo "$line" | awk -v q="$1" -v lo="$2" -v hi="$3" -v seq="$4" \
        '$3 == "test_chr1" && $4 >= lo && $4 <= hi && $5 <= 10 && $6 == "20=" \
         && $10 == seq {
             exit 0
          }
          { print "unexpected record for " q ": " $0 > "/dev/stderr"; exit 1 }' \
        || fail "$1 mapping incorrect"
}
check_repeat_mapped read1 1 128 ACGTACGTACGTACGTACGT
echo "PASS: read1 mapped correctly (repeat, low MAPQ)"
check_repeat_mapped read2 129 256 AACCGGTTAACCGGTTAACC
echo "PASS: read2 mapped correctly (repeat, low MAPQ)"

# read3 has no match in the reference: must not appear as mapped
if awk '$1 == "read3" && $3 != "*" && and($2, 4) == 0' "$TEST_DIR/aln.sam" | grep -q .; then
    fail "read3 unexpectedly mapped"
fi
echo "PASS: absent read not mapped"

echo "=== BAM/BAI output ==="
"$BWA_CPP26" mem "$TEST_DIR/test_idx" "$TEST_DIR/reads.fq" \
    -o "$TEST_DIR/aln.bam" > /dev/null 2> "$TEST_DIR/bam_stderr.txt"
[ -s "$TEST_DIR/aln.bam" ] || { cat "$TEST_DIR/bam_stderr.txt"; fail "empty BAM output"; }
[ -s "$TEST_DIR/aln.bam.bai" ] || fail "missing BAI output"
# The BAM must contain the same mapped records as the SAM.
python3 "$(dirname "$0")/validate_bam.py" "$TEST_DIR/aln.bam" "$TEST_DIR/aln.sam" \
    || fail "BAM/BAI structural validation failed"
echo "PASS: BAM/BAI validated"

echo "=== Paired-end (insert-size estimation) ==="
python3 - "$TEST_DIR" << 'PY'
import random, sys
d = sys.argv[1]
random.seed(7)
ref = ''.join(random.choice('ACGT') for _ in range(2000))
open(d + '/pair_ref.fa', 'w').write('>pair_chr\n' + ref + '\n')
L, ins = 50, 300
f1 = open(d + '/p1.fq', 'w')
f2 = open(d + '/p2.fq', 'w')
rc = str.maketrans('ACGT', 'TGCA')
for i in range(30):
    p = 100 + i * 50
    a = ref[p:p + L]
    b = ref[p + ins - L:p + ins]
    r2 = b.translate(rc)[::-1]
    f1.write(f'@pair{i}/1\n{a}\n+\n' + 'I' * L + '\n')
    f2.write(f'@pair{i}/2\n{r2}\n+\n' + 'I' * L + '\n')
f1.close()
f2.close()
PY
"$BWA_CPP26" index "$TEST_DIR/pair_ref.fa" "$TEST_DIR/pair_idx" > /dev/null
"$BWA_CPP26" mem "$TEST_DIR/pair_idx" "$TEST_DIR/p1.fq" "$TEST_DIR/p2.fq" \
    "$TEST_DIR/pair.sam" > /dev/null 2> "$TEST_DIR/pair_stderr.txt"
[ -s "$TEST_DIR/pair.sam" ] || { cat "$TEST_DIR/pair_stderr.txt"; fail "empty paired SAM"; }
python3 - "$TEST_DIR/pair.sam" << 'PY'
import sys
proper = 0
total = 0
tlen_ok = 0
for line in open(sys.argv[1]):
    if line.startswith('@'):
        continue
    f = line.rstrip('\n').split('\t')
    total += 1
    if int(f[1]) & 0x2:
        proper += 1
        if abs(int(f[8])) == 300:
            tlen_ok += 1
if total != 60:
    sys.exit(f"expected 60 paired records, got {total}")
if proper != 60:
    sys.exit(f"expected 60 proper pairs, got {proper}")
if tlen_ok != 60:
    sys.exit(f"expected TLEN +/-300 for all, got {tlen_ok}")
print("PASS: 30/30 pairs proper with TLEN 300")
PY
# Paired BAM must validate too.
"$BWA_CPP26" mem "$TEST_DIR/pair_idx" "$TEST_DIR/p1.fq" "$TEST_DIR/p2.fq" \
    -o "$TEST_DIR/pair.bam" > /dev/null 2>&1
python3 "$(dirname "$0")/validate_bam.py" "$TEST_DIR/pair.bam" "$TEST_DIR/pair.sam" \
    || fail "paired BAM/BAI validation failed"
echo "PASS: paired BAM/BAI validated"

echo "=== Interleaved pairs (-p) and split flags (-5/-Y) ==="
# Interleave the same pairs into one file; -p must reproduce pair.sam.
python3 - "$TEST_DIR" << 'PY'
import sys
d = sys.argv[1]
r1 = open(d + '/p1.fq').read().split('@')[1:]
r2 = open(d + '/p2.fq').read().split('@')[1:]
with open(d + '/inter.fq', 'w') as f:
    for a, b in zip(r1, r2):
        f.write('@' + a + '@' + b)
PY
"$BWA_CPP26" mem -p "$TEST_DIR/pair_idx" "$TEST_DIR/inter.fq" \
    "$TEST_DIR/inter.sam" > /dev/null 2>&1
diff <(grep -v '^@' "$TEST_DIR/pair.sam") <(grep -v '^@' "$TEST_DIR/inter.sam") \
    || fail "-p output differs from two-file pairing"
echo "PASS: -p interleaved matches paired output"
# -5 and -Y must change the chimera output without breaking validation.
"$BWA_CPP26" mem -5 -Y "$TEST_DIR/pair_idx" "$TEST_DIR/p1.fq" \
    "$TEST_DIR/inter5.sam" > /dev/null 2>&1
[ -s "$TEST_DIR/inter5.sam" ] || fail "-5/-Y run produced nothing"

echo "=== Integration test completed ==="
