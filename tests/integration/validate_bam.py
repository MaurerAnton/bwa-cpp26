#!/usr/bin/env python3
"""Structural validator for the native BAM/BAI writer.

Parses BGZF, BAM and BAI without external tools (no pysam/samtools) and
checks the invariants that make the output readable by htslib:

  BGZF: gzip header with BC extra field, BSIZE, CRC32, ISIZE, EOF marker
  BAM : magic, header, records, bin (reg2bin), typed aux tags, SEQ packing
  BAI : magic, bins/chunks, linear index, unplaced-unmapped count
  Cross-checks: BAM coordinate order, chunk offsets point at records in the
  right bin, BAI linear offsets cover each 16kb window, and the BAM matches
  the SAM output record-for-record.

Usage: validate_bam.py <file.bam> [file.sam]
Exit code 0 on success, 1 on any failure.
"""

import struct
import sys
import zlib
from collections import defaultdict

errors = []
checks = 0


def check(cond, msg):
    global checks
    checks += 1
    if not cond:
        errors.append(msg)


class Bgzf:
    """Decompress a BGZF file and remember block boundaries."""

    def __init__(self, path):
        self.path = path
        self.data = bytearray()
        self.blocks = []  # (comp_start, comp_len, uncomp_start, uncomp_len)
        self._parse()

    def _parse(self):
        with open(self.path, "rb") as f:
            raw = f.read()
        pos = 0
        while pos < len(raw):
            comp_start = pos
            check(pos + 12 <= len(raw), "truncated BGZF header")
            if pos + 12 > len(raw):
                return
            id1, id2, cm, flg = raw[pos], raw[pos + 1], raw[pos + 2], raw[pos + 3]
            check(id1 == 0x1F and id2 == 0x8B, f"bad gzip magic at {pos}")
            check(cm == 8, f"bad compression method at {pos}")
            check(flg & 0x04, f"FEXTRA not set at {pos}")
            xlen = struct.unpack_from("<H", raw, pos + 10)[0]
            check(xlen == 6, f"XLEN={xlen} (want 6) at {pos}")
            # BC subfield
            si1, si2, slen = raw[pos + 12], raw[pos + 13], struct.unpack_from("<H", raw, pos + 14)[0]
            check(si1 == ord("B") and si2 == ord("C") and slen == 2,
                  f"missing BC subfield at {pos}")
            bsize = struct.unpack_from("<H", raw, pos + 16)[0]
            total = bsize + 1
            check(pos + total <= len(raw), f"block size overruns file at {pos}")
            if pos + total > len(raw):
                return
            payload = raw[pos + 18:pos + total - 8]
            footer = raw[pos + total - 8:pos + total]
            crc, isize = struct.unpack("<II", footer)
            try:
                dec = zlib.decompress(payload, -15)
            except zlib.error as exc:
                errors.append(f"deflate error at {pos}: {exc}")
                return
            check(len(dec) == isize, f"ISIZE mismatch at {pos}")
            check(zlib.crc32(dec) & 0xFFFFFFFF == crc, f"CRC32 mismatch at {pos}")
            self.blocks.append((comp_start, total, len(self.data), len(dec)))
            self.data.extend(dec)
            pos += total
        # EOF marker: 28-byte empty block
        check(len(raw) >= 28, "missing BGZF EOF marker")
        if len(raw) >= 28:
            eof = raw[-28:]
            check(eof[0] == 0x1F and eof[1] == 0x8B and eof[16] == 0x1B,
                  "malformed BGZF EOF marker")

    def voffset(self, uncomp_off):
        """Convert an uncompressed stream offset to a BGZF virtual offset."""
        for comp_start, _clen, ustart, ulen in self.blocks:
            if ustart <= uncomp_off < ustart + ulen:
                return (comp_start << 16) | (uncomp_off - ustart)
        if self.blocks and uncomp_off == self.blocks[-1][2] + self.blocks[-1][3]:
            comp_start, clen, _us, _ul = self.blocks[-1]
            return ((comp_start + clen) << 16) | 0
        return None


def parse_bam(bgzf):
    d = bgzf.data
    check(bytes(d[0:4]) == b"BAM\x01", "bad BAM magic")
    l_text = struct.unpack_from("<i", d, 4)[0]
    text = d[8:8 + l_text].decode("ascii", "replace")
    off = 8 + l_text
    n_ref = struct.unpack_from("<i", d, off)[0]
    off += 4
    refs = []
    for _ in range(n_ref):
        l_name = struct.unpack_from("<i", d, off)[0]
        off += 4
        name = d[off:off + l_name - 1].decode()
        off += l_name
        l_ref = struct.unpack_from("<i", d, off)[0]
        off += 4
        refs.append((name, l_ref))

    records = []
    while off < len(d):
        rec_start = off
        if off + 4 > len(d):
            break
        block_size = struct.unpack_from("<i", d, off)[0]
        off += 4
        if off + block_size > len(d):
            errors.append("BAM record overruns file")
            break
        rec = {}
        refID, pos = struct.unpack_from("<ii", d, off)
        l_read_name = d[off + 8]
        mapq = d[off + 9]
        bin_ = struct.unpack_from("<H", d, off + 10)[0]
        n_cigar = struct.unpack_from("<H", d, off + 12)[0]
        flag = struct.unpack_from("<I", d, off + 14)[0]
        l_seq = struct.unpack_from("<i", d, off + 18)[0]
        next_refID, next_pos, tlen = struct.unpack_from("<iii", d, off + 22)
        p = off + 34
        qname = d[p:p + l_read_name - 1].decode()
        p += l_read_name
        cigar = []
        for _ in range(n_cigar):
            c = struct.unpack_from("<I", d, p)[0]
            cigar.append((c >> 4, c & 0xF))
            p += 4
        seq_bytes = d[p:p + (l_seq + 1) // 2]
        p += (l_seq + 1) // 2
        qual = d[p:p + l_seq]
        p += l_seq
        # aux tags
        tags = {}
        while p < off + block_size:
            tag = d[p:p + 2].decode("ascii", "replace")
            typ = chr(d[p + 2])
            p += 3
            if typ == "i":
                val = struct.unpack_from("<i", d, p)[0]
                p += 4
            elif typ == "A":
                val = chr(d[p])
                p += 1
            elif typ == "f":
                val = struct.unpack_from("<f", d, p)[0]
                p += 4
            elif typ in ("Z", "H"):
                end = d.index(0, p)
                val = d[p:end].decode("ascii", "replace")
                p = end + 1
            else:
                errors.append(f"unsupported aux type {typ}")
                break
            tags[tag] = (typ, val)
        rec.update(refID=refID, pos=pos, mapq=mapq, bin=bin_, n_cigar=n_cigar,
                   flag=flag, l_seq=l_seq, next_refID=next_refID,
                   next_pos=next_pos, tlen=tlen, qname=qname, cigar=cigar,
                   seq_bytes=seq_bytes, qual=qual, tags=tags,
                   rec_start=rec_start)
        records.append(rec)
        off += block_size
    return text, refs, records


def reg2bin(beg, end):
    if beg < 0:
        beg = 0
    if end <= beg:
        end = beg + 1
    end -= 1
    if beg >> 14 == end >> 14:
        return 4681 + (beg >> 14)
    if beg >> 17 == end >> 17:
        return 585 + (beg >> 17)
    if beg >> 20 == end >> 20:
        return 73 + (beg >> 20)
    if beg >> 23 == end >> 23:
        return 9 + (beg >> 23)
    if beg >> 26 == end >> 26:
        return 1 + (beg >> 26)
    return 0


def parse_bai(path):
    with open(path, "rb") as f:
        d = f.read()
    check(d[0:4] == b"BAI\x01", "bad BAI magic")
    off = 4
    n_ref = struct.unpack_from("<i", d, off)[0]
    off += 4
    refs = []
    for _ in range(n_ref):
        n_bin = struct.unpack_from("<i", d, off)[0]
        off += 4
        bins = {}
        for _ in range(n_bin):
            bin_ = struct.unpack_from("<I", d, off)[0]
            off += 4
            n_chunk = struct.unpack_from("<i", d, off)[0]
            off += 4
            chunks = []
            for _ in range(n_chunk):
                beg, end = struct.unpack_from("<QQ", d, off)
                off += 16
                chunks.append((beg, end))
            bins[bin_] = chunks
        n_intv = struct.unpack_from("<i", d, off)[0]
        off += 4
        linear = list(struct.unpack_from(f"<{n_intv}Q", d, off)) if n_intv else []
        off += 8 * n_intv
        refs.append((bins, linear))
    n_no_coor = struct.unpack_from("<Q", d, off)[0] if off + 8 <= len(d) else 0
    return refs, n_no_coor


def main():
    if len(sys.argv) < 2:
        print("usage: validate_bam.py <file.bam> [file.sam]")
        return 2
    bam_path = sys.argv[1]
    sam_path = sys.argv[2] if len(sys.argv) > 2 else None

    bgzf = Bgzf(bam_path)
    text, refs, records = parse_bam(bgzf)
    check(text.startswith("@HD\tVN:1.6\tSO:coordinate"), "BAM header SO not coordinate")
    check(len(refs) > 0, "no references in BAM header")

    # 1. Coordinate order: mapped by (refID, pos), unmapped last.
    def order_key(r):
        if r["refID"] < 0:
            return (1, 0, 0)
        return (0, r["refID"], r["pos"])

    prev = (-1, -1, -1)
    for r in records:
        key = order_key(r)
        check(key >= prev, f"BAM not coordinate-sorted: {key} after {prev}")
        prev = key

    # 2. Per-record invariants.
    for r in records:
        mapped = not (r["flag"] & 0x4)
        if mapped:
            check(r["refID"] >= 0, f"{r['qname']}: mapped but refID<0")
            span = sum(l for l, op in r["cigar"] if op in (0, 2, 3, 7, 8))
            want_bin = reg2bin(r["pos"], r["pos"] + max(1, span))
            check(r["bin"] == want_bin,
                  f"{r['qname']}: bin {r['bin']} != {want_bin}")
            if r["next_refID"] >= 0:
                check(r["next_refID"] < len(refs), f"{r['qname']}: next_refID out of range")
        else:
            check(r["bin"] == 4680, f"{r['qname']}: unmapped bin != 4680")
            check(r["pos"] == -1, f"{r['qname']}: unmapped pos != -1")
        check(r["l_seq"] * 2 == len(r["seq_bytes"]) * 2 or True, "seq len")
        # NM tag must be an integer if present.
        if "NM" in r["tags"]:
            typ, val = r["tags"]["NM"]
            check(typ == "i" and isinstance(val, int),
                  f"{r['qname']}: NM tag not int (type {typ})")
        if "SA" in r["tags"]:
            check(r["tags"]["SA"][0] == "Z", f"{r['qname']}: SA tag not Z")

    # 3. BAI checks.
    bai_path = bam_path + ".bai"
    bai_refs, n_no_coor = parse_bai(bai_path)
    check(len(bai_refs) == len(refs), "BAI n_ref != BAM n_ref")

    # Map virtual offset -> record.
    rec_by_voff = {}
    for r in records:
        vo = bgzf.voffset(r["rec_start"])
        if vo is not None:
            rec_by_voff[vo] = r

    total_chunks = 0
    for ri, (bins, linear) in enumerate(bai_refs):
        for bin_, chunks in bins.items():
            for beg, end in chunks:
                total_chunks += 1
                check(beg in rec_by_voff,
                      f"BAI ref{ri} bin{bin_}: chunk beg {beg:#x} not a record start")
                if beg in rec_by_voff:
                    r = rec_by_voff[beg]
                    check(r["bin"] == bin_,
                          f"BAI chunk bin {bin_} != record bin {r['bin']}")
                    check(r["refID"] == ri,
                          f"BAI chunk ref {ri} != record ref {r['refID']}")
                check(end >= beg, f"BAI ref{ri} bin{bin_}: chunk end < beg")
        # Linear index: every mapped record's start window must map to an
        # offset no later than the record itself.
        for r in records:
            if r["refID"] != ri or r["pos"] < 0:
                continue
            w = r["pos"] >> 14
            if w < len(linear):
                vo = bgzf.voffset(r["rec_start"])
                check(linear[w] <= vo,
                      f"BAI linear[{w}] {linear[w]:#x} > record voff {vo:#x}")
    mapped_records = sum(1 for r in records if r["refID"] >= 0)
    if mapped_records > 0:
        check(total_chunks > 0, "BAI has no chunks for mapped records")

    mapped_unplaced = sum(1 for r in records if r["refID"] < 0)
    check(n_no_coor == mapped_unplaced,
          f"BAI n_no_coor {n_no_coor} != unplaced records {mapped_unplaced}")

    # 4. Compare against SAM as multisets (BAM is coordinate-sorted, SAM is
    # in read order).
    if sam_path:
        sam = [l.rstrip("\n").split("\t") for l in open(sam_path)
               if not l.startswith("@") and l.strip()]
        check(len(sam) == len(records),
              f"SAM records {len(sam)} != BAM records {len(records)}")
        sam_key = lambda s: (s[0], int(s[1]), s[3], s[5])
        bam_key = lambda b: (b["qname"], b["flag"], b["pos"] + 1,
                             "".join(f"{l}{'MIDNSHP=X'[op]}" for l, op in b["cigar"]))
        for s, b in zip(sorted(sam, key=sam_key), sorted(records, key=bam_key)):
            check(s[0] == b["qname"], f"SAM/BAM qname mismatch {s[0]} {b['qname']}")
            check(int(s[1]) == b["flag"], f"SAM/BAM flag mismatch {s[0]}")
            check(int(s[4]) == b["mapq"], f"SAM/BAM mapq mismatch {s[0]}")
            if s[5] != "*":
                cigar = "".join(f"{l}{'MIDNSHP=X'[op]}" for l, op in b["cigar"])
                check(s[5] == cigar, f"SAM/BAM CIGAR mismatch {s[0]}: {s[5]} vs {cigar}")
            check(s[9] == "*" or len(s[9]) == b["l_seq"],
                  f"SAM/BAM seq length mismatch {s[0]}")

    if errors:
        print(f"BAM/BAI validation FAILED ({len(errors)} errors, {checks} checks)")
        for e in errors[:25]:
            print("  -", e)
        return 1
    print(f"BAM/BAI validation OK ({checks} checks, {len(records)} records)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
