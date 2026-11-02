#!/usr/bin/env python3
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
# Per-base share of a treated sample's rf-count mutations, skipping sites its untreated control already mutates.
# Modified bases (tRNA m1G37, m2,2G26) and SNPs mutate in both samples and would otherwise swamp the DMS A/C signal.
import argparse
import os
import struct

EOF_MARKER = b"[eofrc]"
NIBBLE = "ACGTN"
BYTE_TO_BASES = ["".join(NIBBLE[min(n, 4)] for n in (b >> 4, b & 0xF)) for b in range(256)]


def read_rc(path):
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        fh.seek(max(size - len(EOF_MARKER), 0))
        # Newer RC files end with uint64 mapped reads, a uint16 version and the marker; older ones with the uint64 only.
        end = size - (8 + 2 + len(EOF_MARKER) if fh.read() == EOF_MARKER else 8)
        fh.seek(0)
        while fh.tell() < end:
            (id_len,) = struct.unpack("<I", fh.read(4))
            tid = fh.read(id_len).rstrip(b"\0").decode()
            (seq_len,) = struct.unpack("<I", fh.read(4))
            seq = "".join(BYTE_TO_BASES[b] for b in fh.read((seq_len + 1) // 2))[:seq_len]
            counts = struct.unpack(f"<{seq_len}I", fh.read(4 * seq_len))
            coverage = struct.unpack(f"<{seq_len}I", fh.read(4 * seq_len))
            fh.read(4)  # per-transcript read count
            yield tid, seq, counts, coverage


def untreated_sites(path, max_rate, min_coverage):
    sites = {}
    for tid, _seq, counts, coverage in read_rc(path):
        mutated = {i for i, (m, c) in enumerate(zip(counts, coverage)) if c >= min_coverage and m > max_rate * c}
        if mutated:
            sites[tid] = mutated
    return sites


def base_mutations(path, skip):
    totals = dict.fromkeys("ACGT", 0)
    for tid, seq, counts, _coverage in read_rc(path):
        masked = skip.get(tid, ())
        for i, (base, m) in enumerate(zip(seq, counts)):
            if m and base in totals and i not in masked:
                totals[base] += m
    return totals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--treated", required=True)
    parser.add_argument("--untreated", required=True)
    parser.add_argument("--max-untreated-rate", type=float, default=0.1)
    parser.add_argument("--min-untreated-coverage", type=int, default=50)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    skip = untreated_sites(args.untreated, args.max_untreated_rate, args.min_untreated_coverage)
    totals = base_mutations(args.treated, skip)
    total = sum(totals.values())
    pct = {b: (100 * totals[b] / total if total else 0.0) for b in "ACGT"}
    with open(args.output, "w") as out:
        out.write("sample\tskipped_sites\tpct_a_muts\tpct_c_muts\tpct_g_muts\tpct_u_muts\n")
        out.write(f"{args.sample}\t{sum(map(len, skip.values()))}\t"
                  + "\t".join(f"{pct[b]:.2f}" for b in "ACGT") + "\n")


if __name__ == "__main__":
    main()
