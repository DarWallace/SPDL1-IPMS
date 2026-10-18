#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
read_raw_headers.py - header reader for Thermo .raw files

Extracts, without Thermo libraries, the acquisition information a Data
Descriptor must report in Methods: instrument model, serial number, firmware
version, creation date, instrument-method name and the method text (LC + MS),
from which the gradient, DDA Top-N, resolutions, AGC targets, injection times
and isolation window are read.

Works on:
  1) loose .raw files on disk
  2) .raw files inside a .zip (including Deflate64, produced by Windows
     Explorer for large folders; requires `pip install zipfile-deflate64`)

Only the header is read (by default the first 48 MB of each .raw), where
Thermo stores RawFileHeader + SequencerRow + AutoSamplerInfo + RawFileInfo
+ InstrumentMethod. Spectra are not read.

Usage:
    python read_raw_headers.py <file.zip | file.raw | folder> [--json out.json]
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime, timedelta

HEADER_BYTES = 48 * 1024 * 1024  # 48 MB: covers header + instrument method

# ---------------------------------------------------------------- helpers

_UTF16_RE = re.compile(rb"(?:[\x20-\x7e]\x00){4,}")
_ASCII_RE = re.compile(rb"[\x20-\x7e]{6,}")


def utf16_strings(buf):
    """Printable UTF-16LE strings of 4 or more characters."""
    return [m.group().decode("utf-16-le", "ignore") for m in _UTF16_RE.finditer(buf)]


def ascii_strings(buf):
    return [m.group().decode("ascii", "ignore") for m in _ASCII_RE.finditer(buf)]


def dedup(seq):
    seen, out = set(), []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def parse_raw_file_header(buf):
    """Thermo header: magic 0x01A1 + 'Finnigan' in UTF-16 + version + date.

    The exact layout changes between file versions; only stable fields are
    read and anything that cannot be confirmed is flagged.
    """
    info = {}
    if len(buf) < 1400:
        info["signature_ok"] = False
        return info
    info["signature_ok"] = buf[0:2] == b"\x01\xa1" and buf[2:18].decode(
        "utf-16-le", "ignore"
    ).startswith("Finnigan")
    try:
        info["file_format_version"] = int.from_bytes(buf[36:40], "little")
    except Exception:
        pass
    # Windows SYSTEMTIME timestamp (16 bytes) at offset 40 in most
    # file versions: year, month, dayofweek, day, hour, minute, second, ms
    try:
        vals = [int.from_bytes(buf[40 + 2 * i: 42 + 2 * i], "little") for i in range(8)]
        year, month, _dow, day, hour, minute, second, msec = vals
        if 1990 <= year <= 2100 and 1 <= month <= 12 and 1 <= day <= 31:
            info["creation_date"] = datetime(
                year, month, day, hour % 24, minute % 60, second % 60,
                (msec % 1000) * 1000
            ).isoformat()
    except Exception:
        pass
    # Fixed-length UTF-16 tags after the header (user, computer...)
    tags = []
    for off in (52, 564, 1076):
        try:
            t = buf[off: off + 512].decode("utf-16-le", "ignore").split("\x00")[0].strip()
            if t:
                tags.append(t)
        except Exception:
            pass
    if tags:
        info["header_tags"] = dedup(tags)
    return info


# Patterns of interest in the instrument method / header
PATTERNS = {
    "instrument_model": [
        r"\bQ\s*Exactive(?:\s+(?:Plus|HF-?X?|HF))?\b",
        r"\bOrbitrap\s+\w+(?:\s+\w+)?\b",
        r"\bLTQ[\w\- ]*\b",
        r"\bExploris\s*\d+\b",
        r"\bFusion\s+\w+\b",
        r"\bTSQ\s+\w+\b",
    ],
    "lc_system": [
        r"\bEASY-?nLC\s*\d*\b",
        r"\bUltiMate\s*\d+\b",
        r"\bnanoLC[\w\- ]*\b",
        r"\bVanquish[\w\- ]*\b",
        r"\bAccela\b",
    ],
    "serial_number": [
        r"(?:Serial\s*(?:No|Number|#)\s*[:=]?\s*)([A-Za-z0-9\-]{3,})",
        r"\bExactive Series slot #\d+\b",
    ],
    "firmware": [
        r"(?:Instrument (?:software|firmware) version|Firmware)\s*[:=]?\s*([\w\.\-]+)",
        r"\bTune\s+(?:version|Application)\s*[:=]?\s*([\w\.\- ]+)",
        r"\bXcalibur\s*[\w\.\- ]*",
    ],
    "method_file": [
        r"[A-Za-z]:\\[^\r\n\x00]{3,120}\.meth\b",
        r"\b[\w\-\. ]{3,80}\.meth\b",
    ],
    "dda_topn": [
        r"\bTop\s?N?\s*=?\s*(\d{1,2})\b",
        r"\bLoop\s*[Cc]ount\s*[:=]?\s*(\d{1,2})\b",
        r"\bNumber of (?:dependent|data dependent) scans?\s*[:=]?\s*(\d{1,2})\b",
    ],
    "resolution": [r"\bResolution\s*[:=]?\s*([\d,\.]+)\b"],
    "agc_target": [r"\bAGC\s*(?:target)?\s*[:=]?\s*([\d\.eE\+]+)\b"],
    "max_it": [r"\bMaximum (?:IT|injection time)\s*[:=]?\s*([\d\.]+)\s*(?:ms)?\b"],
    "isolation_window": [r"\bIsolation window\s*[:=]?\s*([\d\.]+)\s*(?:m/z|Th)?\b"],
    "nce": [r"\b(?:NCE|Normalized [Cc]ollision [Ee]nergy)\s*[:=]?\s*([\d\.,\s]+)\b"],
    "scan_range": [r"\bScan range\s*[:=]?\s*([\d\.]+)\s*(?:to|-|–)\s*([\d\.]+)"],
    "run_time_min": [
        r"\b(?:MS )?Run [Tt]ime\s*\(?min\)?\s*[:=]?\s*([\d\.]+)",
        r"\bTotal (?:run )?time\s*[:=]?\s*([\d\.]+)\s*min",
        r"\bGradient (?:length|time)\s*[:=]?\s*([\d\.]+)\s*min",
    ],
    "column": [r"\b(?:Column|Analytical column)\s*[:=]?\s*([^\r\n\x00]{4,80})"],
    "sample_name": [
        r"(?:Sample\s*(?:Name|ID)\s*[:=]\s*)([^\r\n\x00]{1,60})",
    ],
}


def scan_patterns(text):
    found = {}
    for key, pats in PATTERNS.items():
        hits = []
        for p in pats:
            for m in re.finditer(p, text, flags=re.IGNORECASE):
                hits.append(m.group(0).strip() if m.lastindex is None else m.group(0).strip())
        if hits:
            found[key] = dedup(hits)[:8]
    return found


def extract_method_text(strings):
    """Block of lines that appear to belong to the instrument method."""
    kw = (
        "gradient", "flow", "%b", "% b", "solvent", "loop", "nce", "agc",
        "resolution", "isolation", "microscan", "polarity", "scan range",
        "injection", "column", "min", "topn", "top n", "dda", "full ms",
        "ddms2", "dd-ms2", "charge", "dynamic exclusion",
    )
    keep = []
    for s in strings:
        low = s.lower()
        if any(k in low for k in kw) and len(s) < 200:
            keep.append(s.strip())
    return dedup(keep)


def analyse_buffer(name, buf):
    rec = {"file": name, "bytes_read": len(buf)}
    rec.update(parse_raw_file_header(buf))
    u16 = utf16_strings(buf)
    asc = ascii_strings(buf)
    all_text = "\n".join(u16 + asc)
    rec["matches"] = scan_patterns(all_text)
    rec["method_lines"] = extract_method_text(u16 + asc)[:400]
    # Most informative short strings, for manual inspection
    rec["notable_strings"] = dedup(
        [s.strip() for s in u16 if 4 <= len(s) <= 90 and re.search(r"[A-Za-z]{3}", s)]
    )[:400]
    return rec


# ---------------------------------------------------------------- inputs

def iter_zip(path, limit=HEADER_BYTES):
    try:
        import zipfile_deflate64 as zf
    except ImportError:
        import zipfile as zf
        print(
            "[warning] zipfile_deflate64 not installed; Deflate64 zips will fail. "
            "Install with: pip install zipfile-deflate64",
            file=sys.stderr,
        )
    with zf.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.lower().endswith(".raw")]
        print(f"[info] {len(names)} raw files in {os.path.basename(path)}", file=sys.stderr)
        for n in names:
            size = z.getinfo(n).file_size
            print(f"[info] reading header of {n} ({size/1e9:.2f} GB)", file=sys.stderr)
            with z.open(n) as fh:
                buf = fh.read(min(limit, size))
            yield n, buf, size


def iter_disk(path, limit=HEADER_BYTES):
    if os.path.isfile(path):
        files = [path]
    else:
        files = []
        for root, _d, fs in os.walk(path):
            for f in fs:
                if f.lower().endswith(".raw"):
                    files.append(os.path.join(root, f))
    for f in sorted(files):
        size = os.path.getsize(f)
        print(f"[info] reading header of {f} ({size/1e9:.2f} GB)", file=sys.stderr)
        with open(f, "rb") as fh:
            buf = fh.read(min(limit, size))
        yield f, buf, size


def main():
    ap = argparse.ArgumentParser(description="Thermo .raw header reader")
    ap.add_argument("path", help="zip, folder or .raw file")
    ap.add_argument("--json", help="write the full report to JSON")
    ap.add_argument("--bytes", type=int, default=HEADER_BYTES,
                    help="header bytes to read per file")
    args = ap.parse_args()

    src = iter_zip if args.path.lower().endswith(".zip") else iter_disk
    report = []
    for name, buf, size in src(args.path, args.bytes):
        rec = analyse_buffer(name, buf)
        rec["file_size_bytes"] = size
        report.append(rec)
        print("=" * 78)
        print(name, f"({size/1e9:.2f} GB)")
        print("  Thermo signature:", rec.get("signature_ok"))
        print("  creation date:", rec.get("creation_date", "unreadable"))
        for k in ("instrument_model", "lc_system", "serial_number", "firmware",
                  "method_file", "dda_topn", "run_time_min", "resolution",
                  "agc_target", "max_it", "isolation_window", "nce",
                  "scan_range", "column"):
            v = rec["matches"].get(k)
            if v:
                print(f"  {k}: {v}")
        if rec.get("header_tags"):
            print("  header tags:", rec["header_tags"])

    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(report, fh, ensure_ascii=False, indent=2)
        print(f"[info] JSON report written to {args.json}", file=sys.stderr)


if __name__ == "__main__":
    main()
