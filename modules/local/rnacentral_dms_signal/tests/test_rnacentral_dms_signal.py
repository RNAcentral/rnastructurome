import struct
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[4] / "bin" / "rnacentral_dms_signal.py"


def write_rc(path, records, eof_marker=True):
    code = {"A": 0, "C": 1, "G": 2, "T": 3}
    out = bytearray()
    for tid, seq, counts, coverage in records:
        name = tid.encode() + b"\0"
        out += struct.pack("<I", len(name)) + name + struct.pack("<I", len(seq))
        nibbles = [code[b] for b in seq] + [0] * (len(seq) % 2)
        out += bytes((nibbles[i] << 4) | nibbles[i + 1] for i in range(0, len(nibbles), 2))
        out += struct.pack(f"<{len(seq)}I", *counts) + struct.pack(f"<{len(seq)}I", *coverage)
        out += struct.pack("<I", max(coverage))
    out += struct.pack("<Q", 1000)
    if eof_marker:
        out += struct.pack("<H", 1) + b"[eofrc]"
    path.write_bytes(bytes(out))


def run(tmp_path, treated, untreated, *extra):
    output = tmp_path / "out.tsv"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--sample", "s1", "--treated", str(treated),
         "--untreated", str(untreated), "--output", str(output), *extra],
        check=True,
    )
    header, row = output.read_text().splitlines()
    return dict(zip(header.split("\t"), row.split("\t")))


def test_skips_sites_the_untreated_control_mutates(tmp_path):
    # Position 3 is a modified G, mutated in both samples; it outweighs the A/C signal until skipped.
    treated, untreated = tmp_path / "t.rc", tmp_path / "u.rc"
    write_rc(treated, [("tx1", "ACGGT", [30, 20, 0, 90, 0], [100] * 5)])
    write_rc(untreated, [("tx1", "ACGGT", [1, 1, 0, 40, 1], [100] * 5)])
    row = run(tmp_path, treated, untreated)
    assert row["skipped_sites"] == "1"
    shares = [row[f"pct_{b}_muts"] for b in "acgu"]
    assert shares == ["60.00", "40.00", "0.00", "0.00"]


def test_without_skipping_matches_plain_base_shares(tmp_path):
    treated, untreated = tmp_path / "t.rc", tmp_path / "u.rc"
    write_rc(treated, [("tx1", "ACGGT", [30, 20, 0, 90, 0], [100] * 5)], eof_marker=False)
    write_rc(untreated, [("tx1", "ACGGT", [1, 1, 0, 40, 1], [100] * 5)], eof_marker=False)
    row = run(tmp_path, treated, untreated, "--max-untreated-rate", "1")
    assert row["skipped_sites"] == "0"
    assert row["pct_g_muts"] == "64.29"


def test_thin_untreated_coverage_does_not_skip(tmp_path):
    treated, untreated = tmp_path / "t.rc", tmp_path / "u.rc"
    write_rc(treated, [("tx1", "ACG", [30, 20, 90], [100] * 3), ("tx2", "AC", [5, 5], [50, 50])])
    write_rc(untreated, [("tx1", "ACG", [0, 0, 20], [30, 30, 30])])
    row = run(tmp_path, treated, untreated)
    assert row["skipped_sites"] == "0"
    assert row["pct_g_muts"] == "60.00"
