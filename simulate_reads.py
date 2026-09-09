#!/usr/bin/env python3
"""
Read simulator for BWA validation.
Generates synthetic reads from a reference with known ground truth.
"""

import sys
import random
import gzip
import os
import subprocess
import hashlib
import math
import argparse
from collections import defaultdict

# IUPAC codes for N
NUCLEOTIDES = ['A', 'C', 'G', 'T']
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

def reverse_complement(seq):
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))

def mutate_read(seq, error_rate=0.01, indel_rate=0.001):
    """Introduce sequencing errors."""
    result = []
    i = 0
    while i < len(seq):
        r = random.random()
        if r < error_rate:
            # Substitution
            result.append(random.choice([b for b in NUCLEOTIDES if b != seq[i]]))
            i += 1
        elif r < error_rate + indel_rate:
            # Indel
            if random.random() < 0.5:
                # Deletion
                i += 1
            else:
                # Insertion
                result.append(random.choice(NUCLEOTIDES))
        else:
            result.append(seq[i])
            i += 1
    return ''.join(result)

def simulate_reads(reference, num_reads, read_len, error_rate=0.01, seed=42):
    """Generate simulated reads with ground truth."""
    random.seed(seed)
    ref_len = len(reference)
    
    # Precompute unique 20-mers to find mappable regions
    kmer_size = 20
    kmer_counts = {}
    for i in range(len(reference) - kmer_size + 1):
        kmer = reference[i:i+kmer_size]
        kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    
    # Find uniquely mappable positions (all 20-mers in the read are unique)
    mappable = [False] * ref_len
    for i in range(ref_len - read_len + 1):
        unique = True
        for k in range(read_len - kmer_size + 1):
            kmer = reference[i+k:i+k+kmer_size]
            if kmer_counts.get(kmer, 0) > 1:
                unique = False
                break
        if unique:
            for j in range(read_len):
                mappable[i+j] = True
    
    # Find all valid start positions
    valid_starts = [i for i in range(ref_len - read_len + 1) if mappable[i]]
    
    if not valid_starts:
        # Fallback: use all positions
        valid_starts = list(range(ref_len - read_len + 1))
        print(f"Warning: No uniquely mappable regions found, using all positions")
    
    reads = []
    ground_truth = []
    
    for _ in range(num_reads):
        pos = random.choice(valid_starts)
        forward = random.choice([True, False])
        
        if forward:
            true_seq = reference[pos:pos+read_len]
            true_pos = pos
        else:
            true_seq = reverse_complement(reference[pos:pos+read_len])
            true_pos = pos
        
        read_seq = mutate_read(true_seq, error_rate)
        
        # Quality string
        qual = ''.join(chr(min(40, max(2, int(-10 * math.log10(error_rate) + random.gauss(0, 5))))) for _ in read_seq)
        
        reads.append({
            'name': f'read_{len(reads)}',
            'seq': read_seq,
            'qual': qual,
            'forward': forward
        })
        ground_truth.append({
            'pos': true_pos,
            'forward': forward,
            'true_seq': true_seq
        })
    
    return reads, ground_truth

def write_fastq(reads, path):
    with open(path, 'w') as f:
        for r in reads:
            f.write(f"@{r['name']}\n{r['seq']}\n+\n{r['qual']}\n")

def write_fasta(ref, path, name='ref'):
    with open(path, 'w') as f:
        f.write(f'>{name}\n')
        for i in range(0, len(ref), 80):
            f.write(ref[i:i+80] + '\n')

def run_aligner(binary, index_prefix, fastq, output_sam, threads=1):
    """Run aligner and return SAM output."""
    cmd = [binary, 'mem', '-t', str(threads), index_prefix, fastq]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    with open(output_sam, 'w') as f:
        f.write(result.stdout)
    return result.returncode, result.stderr

def parse_sam(sam_path):
    """Parse SAM file into aligned records."""
    records = []
    with open(sam_path) as f:
        for line in f:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            records.append({
                'qname': fields[0],
                'flag': int(fields[1]),
                'rname': fields[2],
                'pos': int(fields[3]),
                'mapq': int(fields[4]),
                'cigar': fields[5],
                'seq': fields[9]
            })
    return records

def evaluate_alignment(ground_truth, sam_records, ref_name):
    """Compare alignments against ground truth."""
    results = {
        'total': len(ground_truth),
        'mapped': 0,
        'correct_pos': 0,
        'correct_strand': 0,
        'correct_cigar': 0,
        'mapq_sum': 0,
        'unmapped': 0,
        'wrong_pos': 0,
        'wrong_strand': 0
    }
    
    # Build lookup by read name
    sam_by_name = {r['qname']: r for r in sam_records}
    
    for i, gt in enumerate(ground_truth):
        read_name = f'read_{i}'
        sam = sam_by_name.get(read_name)
        
        if not sam or sam['flag'] & 4:  # unmapped
            results['unmapped'] += 1
            continue
        
        results['mapped'] += 1
        results['mapq_sum'] += sam['mapq']
        
        # Check position (1-based in SAM)
        if sam['pos'] == gt['pos'] + 1:
            results['correct_pos'] += 1
        else:
            results['wrong_pos'] += 1
        
        # Check strand (flag 16 = reverse)
        sam_forward = not (sam['flag'] & 16)
        if sam_forward == gt['forward']:
            results['correct_strand'] += 1
        else:
            results['wrong_strand'] += 1
        
        # CIGAR check (perfect match = all '=')
        if sam['cigar'] == f'{len(ground_truth[0]["true_seq"])}=':
            results['correct_cigar'] += 1
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Validate BWA implementation against ground truth')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--num-reads', type=int, default=10000)
    parser.add_argument('--read-len', type=int, default=150)
    parser.add_argument('--error-rate', type=float, default=0.01)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--bwa-cpp', default='./build/src/bwa')
    parser.add_argument('--bwa-sys', default='bwa')
    parser.add_argument('--out-dir', default='/tmp/bwa_validation')
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Load reference
    ref_seq = ''
    with open(args.ref) as f:
        for line in f:
            if not line.startswith('>'):
                ref_seq += line.strip()
    
    print(f"Reference length: {len(ref_seq)}")
    print(f"Generating {args.num_reads} reads of length {args.read_len}...")
    
    # Simulate reads
    reads, ground_truth = simulate_reads(ref_seq, args.num_reads, args.read_len, args.error_rate, args.seed)
    print(f"Generated {len(reads)} reads")
    
    # Write FASTA
    fastq_path = os.path.join(args.out_dir, 'reads.fq')
    write_fastq(reads, fastq_path)
    
    # Write reference FASTA (use the same reference file for indexing)
    ref_fa = args.ref
    import shutil
    shutil.copy2(args.ref, os.path.join(args.out_dir, 'ref.fa'))
    ref_fa = os.path.join(args.out_dir, 'ref.fa')
    
    # Build indexes
    idx_prefix = os.path.join(args.out_dir, 'idx')
    print("Building index with bwa-cpp26...")
    subprocess.run([args.bwa_cpp, 'index', ref_fa, idx_prefix], check=True)
    
    # Also build with system bwa if available
    has_sys_bwa = False
    try:
        subprocess.run([args.bwa_sys, 'index', ref_fa, idx_prefix + '_sys'], check=True, capture_output=True)
        has_sys_bwa = True
    except:
        print("System bwa not available, skipping comparison")
    
    # Align with bwa-cpp26
    print("Aligning with bwa-cpp26...")
    sam_cpp = os.path.join(args.out_dir, 'aln_cpp.sam')
    run_aligner(args.bwa_cpp, idx_prefix, fastq_path, sam_cpp)
    
    sam_records_cpp = parse_sam(sam_cpp)
    results_cpp = evaluate_alignment(ground_truth, sam_records_cpp, 'cpp')
    
    # Align with system bwa if available
    results_sys = None
    if has_sys_bwa:
        print("Aligning with system bwa...")
        sam_sys = os.path.join(args.out_dir, 'aln_sys.sam')
        run_aligner(args.bwa_sys, idx_prefix + '_sys', fastq_path, sam_sys)
        sam_records_sys = parse_sam(sam_sys)
        results_sys = evaluate_alignment(ground_truth, sam_records_sys, 'sys')
    
    # Print results
    print("\n=== RESULTS ===")
    for name, res in [('bwa-cpp26', results_cpp), ('system bwa', results_sys)]:
        if not res:
            print(f"{name}: N/A")
            continue
        print(f"\n{name}:")
        print(f"  Total reads:      {res['total']}")
        print(f"  Mapped:           {res['mapped']} ({100*res['mapped']/res['total']:.1f}%)")
        if res['mapped'] > 0:
            print(f"  Correct position: {res['correct_pos']} ({100*res['correct_pos']/res['mapped']:.1f}%)")
            print(f"  Correct strand:   {res['correct_strand']} ({100*res['correct_strand']/res['mapped']:.1f}%)")
            print(f"  Perfect CIGAR:    {res['correct_cigar']} ({100*res['correct_cigar']/res['mapped']:.1f}%)")
            print(f"  Avg MAPQ:         {res['mapq_sum']/max(1,res['mapped']):.1f}")
            print(f"  Unmapped:         {res['unmapped']}")
            print(f"  Wrong position:   {res['wrong_pos']}")
            print(f"  Wrong strand:     {res['wrong_strand']}")
        else:
            print(f"  Unmapped:         {res['unmapped']}")
            print(f"  Wrong position:   {res['wrong_pos']}")
            print(f"  Wrong strand:     {res['wrong_strand']}")

if __name__ == '__main__':
    main()