#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_hspindly_ptms.py
========================
Reproducible analysis pipeline for the immunoprecipitation / LC-MS/MS
characterisation of post-translational modifications (PTMs) of human
Spindly (SPDL1, UniProt Q96EA4) in prometaphase-arrested cells.

The script parses the grouped protein/peptide exports produced by
Thermo Proteome Discoverer (Sequest HT search engine), extracts every
phosphorylation (Phospho) and ubiquitination (di-glycine, "GG") event
mapped onto Spindly, converts peptide-local modification coordinates to
absolute residue numbers of the canonical Q96EA4 sequence, compares the
experimental conditions (nocodazole [NZ], nocodazole + okadaic acid
[NZ+OK], and the lambda-phosphatase control [NZ+PPase]) and, in parallel,
builds the co-immunoprecipitated interactome table for each run.

Design goals
------------
* Zero third-party services at run time: the Q96EA4 sequence is embedded
  (MD5-checked) so the pipeline is fully reproducible offline.
* Deterministic: same inputs -> identical outputs.
* Transparent provenance: every reported site keeps the peptide, the
  Sequest HT confidence, the XCorr and the source file.

Inputs   : the three experiment folders shipped with the study.
Outputs  : CSV tables + summary matrix + publication figures (PNG/PDF).

Usage
-----
    python analyze_hspindly_ptms.py --data-dir DATA --out-dir results

Dependencies: openpyxl, matplotlib  (both pip-installable).

License: MIT
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
from dataclasses import dataclass
from typing import Optional

import openpyxl

# ---------------------------------------------------------------------------
# 0. Reference sequence (UniProt Q96EA4, SPDLY_HUMAN, isoform 1, 605 aa, SV=2)
# ---------------------------------------------------------------------------
SPINDLY_ACC = "Q96EA4"
SPINDLY_SEQ = (
    "MEADIITNLRCRLKEAEEERLKAAQYGLQLVESQNELQNQLDKCRNEMMTMTESYEQEKY"
    "TLQREVELKSRMLESLSCECEAIKQQQKMHLEKLEEQLSRSHGQEVNELKTKIEKLKVEL"
    "DEARLSEKQLKHQVDHQKELLSCKSEELRVMSERVQESMSSEMLALQIELTEMESMKTTL"
    "KEEVNELQYRQEQLELLITNLMRQVDRLKEEKEEREKEAVSYYNALEKARVANQDLQVQL"
    "DQALQQALDPNSKGNSLFAEVEDRRAAMERQLISMKVKYQSLKKQNVFNREQMQRMKLQI"
    "ATLLQMKGSQTEFEQQERLLAMLEQKNGEIKHLLGEIRNLEKFKNLYDSMESKPSVDSGT"
    "LEDNTYYTDLLQMKLDNLNKEIESTKGELSIQRMKALFESQRALDIERKLFANERCLQLS"
    "ESENMKLRAKLDELKLKYEPEETVEVPVLKKRREVLPVDITTAKDACVNNSALGGEVYRL"
    "PPQKEETQSCPNSLEDNNLQLEKSVSIYTPVVSLSPHKNLPVDMQLKKEKKCVKLIGVPA"
    "DAEALSERSGNTPNSPRLAAESKLQTEVKEGKETSSKLEKETCKKLHPILYVSSKSTPET"
    "QCPQQ"
)
assert len(SPINDLY_SEQ) == 605, "Unexpected Spindly sequence length"
SEQ_MD5 = hashlib.md5(SPINDLY_SEQ.encode()).hexdigest()

# ---------------------------------------------------------------------------
# 1. Experiment manifest
#    Maps every source file to (cell line, treatment, PTM search, batch).
#    Treatments: NZ = nocodazole; NZ_OK = nocodazole + okadaic acid;
#                NZ_PPase = nocodazole, immunoprecipitate treated with lambda phosphatase (control).
# ---------------------------------------------------------------------------
MANIFEST = [
    # HCT116  batch 2022-03-24
    ("HCT116 NZ y NZ+OK 22.03.22/NZ_240322_vs_UniprotHuman-PHOSPHO.xlsx",     "HCT116", "NZ",       "phospho", "2022-03-24"),
    ("HCT116 NZ y NZ+OK 22.03.22/NZ-OK_240322_vs_UniprotHuman-PHOSPHO.xlsx",  "HCT116", "NZ_OK",    "phospho", "2022-03-24"),
    ("HCT116 NZ y NZ+OK 22.03.22/NZ_240322_vs_UniprotHuman-UBI.xlsx",         "HCT116", "NZ",       "ubi",     "2022-03-24"),
    ("HCT116 NZ y NZ+OK 22.03.22/NZ-OK_240322_vs_UniprotHuman-UBI.xlsx",      "HCT116", "NZ_OK",    "ubi",     "2022-03-24"),
    # RPE-1   batch 2022-02-18
    ("RPE-1 NZ y NZ+OK 22.02.22/NZ_vs_UniprotHuman-phospho.xlsx",             "RPE-1",  "NZ",       "phospho", "2022-02-18"),
    ("RPE-1 NZ y NZ+OK 22.02.22/NZ-OK_vs_UniprotHuman-phospho.xlsx",          "RPE-1",  "NZ_OK",    "phospho", "2022-02-18"),
    ("RPE-1 NZ y NZ+OK 22.02.22/NZ_vs_uniprotHuman-Ubi.xlsx",                 "RPE-1",  "NZ",       "ubi",     "2022-02-18"),
    ("RPE-1 NZ y NZ+OK 22.02.22/NZ-OK_vs_UniprotHuman-Ubi.xlsx",              "RPE-1",  "NZ_OK",    "ubi",     "2022-02-18"),
    # HCT116  phosphatase-control batch (acquired 2022-01-27/29; folder named after export date)
    ("HCT116 NZ+PPasa y NZ+OK 08.02.22/NZPPasa-WashingBeadSol_vs_UniprotHuman_PTM-P.xlsx",   "HCT116", "NZ_PPase", "phospho", "2022-01-27"),
    ("HCT116 NZ+PPasa y NZ+OK 08.02.22/NZ-OK-WashingBeadSol_vs_UniprotHuman_PTM-P.xlsx",     "HCT116", "NZ_OK",    "phospho", "2022-01-27"),
    ("HCT116 NZ+PPasa y NZ+OK 08.02.22/NZPPasa-WashingBeadSol_vs_UniprotHuman_PTM-UBI.xlsx", "HCT116", "NZ_PPase", "ubi",     "2022-01-27"),
    ("HCT116 NZ+PPasa y NZ+OK 08.02.22/NZ-OK-WashingBeadSol_vs_UniprotHuman_PTM-UBI.xlsx",   "HCT116", "NZ_OK",    "ubi",     "2022-01-27"),
]

# Column layout of the peptide-level block (Proteome Discoverer, Sequest HT)
PEP_COLS = {
    "sequence": 1,
    "modifications": 2,
    "n_protein_groups": 3,
    "master_accessions": 4,
    "missed_cleavages": 5,
    "theo_mhplus": 6,
    "confidence": 7,
    "xcorr": 8,
}
# Column layout of the protein-level block
PROT_COLS = {
    "accession": 0,
    "description": 1,
    "coverage": 2,
    "n_peptides": 3,
    "n_aas": 4,
    "mw_kda": 5,
    "score": 6,
}

MOD_RE = re.compile(r"(\d+)x([A-Za-z0-9]+)\s*\[([^\]]+)\]")
POS_RE = re.compile(r"^([A-Za-z])(\d+)$")


# ---------------------------------------------------------------------------
# 2. Data model
# ---------------------------------------------------------------------------
@dataclass
class Peptide:
    sequence: str
    modifications: str
    master_accessions: str
    confidence: Optional[str]
    xcorr: Optional[float]
    protein_acc: str
    protein_desc: str


@dataclass
class Site:
    ptm: str            # "phospho" or "ubi"
    residue: str        # amino-acid letter as reported
    position: int       # absolute position in Q96EA4
    localised: bool     # True if the search engine placed the mod on a single residue
    peptide: str
    confidence: Optional[str]
    xcorr: Optional[float]
    cell_line: str
    treatment: str
    batch: str
    source_file: str


@dataclass
class ProteinHit:
    accession: str
    description: str
    coverage: Optional[float]
    n_peptides: Optional[int]
    score: Optional[float]
    cell_line: str
    treatment: str
    batch: str
    source_file: str


# ---------------------------------------------------------------------------
# 3. Parsing
# ---------------------------------------------------------------------------
def parse_pd_export(path: str):
    """Parse a grouped Proteome Discoverer 'Proteins' export.

    Layout of the export (single 'Proteins' sheet):
      * one column header row with 'Accession' in column A;
      * protein-level rows: column A holds the accession, remaining columns
        hold Description / Coverage / #Peptides / #AAs / MW / Score;
      * a peptide-level header row: column A empty, column B == 'Sequence';
      * peptide-level rows: column A empty, column B holds the peptide
        sequence, remaining columns hold Modifications / #Protein Groups /
        Master Protein Accessions / #Missed Cleavages / Theo.MH+ /
        Confidence / XCorr.
    Peptides are attributed to the most recent protein-level row.

    Returns (proteins, peptides): the list of protein-level rows and a flat
    list of Peptide objects, each tagged with its parent protein.
    """
    wb = openpyxl.load_workbook(path, data_only=True, read_only=False)
    ws = wb["Proteins"]
    rows = list(ws.iter_rows(values_only=True))

    proteins, peptides = [], []
    cur_acc = cur_desc = None

    for r in rows:
        a = str(r[0]).strip() if (r and r[0] is not None) else ""
        b = r[1] if (r and len(r) > 1) else None

        if a == "Accession":
            continue                       # protein-block column header
        if a:                              # protein-level row (accession set)
            cur_acc = a
            cur_desc = str(r[PROT_COLS["description"]] or "").strip()
            proteins.append(r)
            continue
        if b == "Sequence":
            continue                       # peptide-level column header
        if b not in (None, "") and cur_acc:
            seq = str(b).strip()
            if not re.fullmatch(r"[A-Za-z]+", seq):
                continue
            xc = r[PEP_COLS["xcorr"]] if len(r) > PEP_COLS["xcorr"] else None
            try:
                xc = float(xc) if xc not in (None, "") else None
            except (TypeError, ValueError):
                xc = None
            peptides.append(Peptide(
                sequence=seq,
                modifications=str(r[PEP_COLS["modifications"]] or ""),
                master_accessions=str(r[PEP_COLS["master_accessions"]] or ""),
                confidence=(str(r[PEP_COLS["confidence"]]).strip()
                            if r[PEP_COLS["confidence"]] not in (None, "") else None),
                xcorr=xc,
                protein_acc=cur_acc,
                protein_desc=cur_desc,
            ))

    return proteins, peptides


# ---------------------------------------------------------------------------
# 4. Modification -> residue mapping
# ---------------------------------------------------------------------------
def peptide_offsets(peptide_seq: str, protein_seq: str):
    """Return every 0-based start offset of a peptide within the protein."""
    offsets, start = [], protein_seq.find(peptide_seq)
    while start != -1:
        offsets.append(start)
        start = protein_seq.find(peptide_seq, start + 1)
    return offsets


def extract_sites(pep: Peptide, ptm_kind: str, protein_seq: str, meta: dict):
    """Convert the Modifications string of a Spindly peptide into absolute Sites.

    ptm_kind : "phospho" -> look for 'Phospho'; "ubi" -> look for 'GG'.
    Ambiguous placements (e.g. 'Phospho [S/T]', 'GG [K/S]') are recorded as
    non-localised sites at the peptide N-terminal anchor and flagged.
    """
    target = "Phospho" if ptm_kind == "phospho" else "GG"
    sites = []
    offsets = peptide_offsets(pep.sequence, protein_seq)
    if not offsets:
        return sites          # peptide not on the canonical isoform (skip)
    off = offsets[0]          # unique in practice for Spindly tryptic peptides

    for count, mod_name, inside in MOD_RE.findall(pep.modifications):
        if mod_name != target:
            continue
        for token in [t.strip() for t in inside.split(";") if t.strip()]:
            m = POS_RE.match(token)
            if m:                                   # localised, e.g. S7
                local = int(m.group(2))
                pos = off + local                   # 1-based absolute position
                seq_aa = protein_seq[pos - 1]
                sites.append(Site(ptm_kind, seq_aa, pos, True,
                                  pep.sequence, pep.confidence, pep.xcorr,
                                  meta["cell_line"], meta["treatment"],
                                  meta["batch"], meta["source_file"]))
            else:                                   # ambiguous, e.g. S/T or K/S
                sites.append(Site(ptm_kind, token, off + 1, False,
                                  pep.sequence, pep.confidence, pep.xcorr,
                                  meta["cell_line"], meta["treatment"],
                                  meta["batch"], meta["source_file"]))
    return sites


# ---------------------------------------------------------------------------
# 5. Main driver
# ---------------------------------------------------------------------------
def run(data_dir: str, out_dir: str, min_xcorr: float = 0.0,
        high_conf_only: bool = False):
    os.makedirs(out_dir, exist_ok=True)
    all_sites = []
    interactome = []
    run_summary = []

    for rel, cell, treat, ptm, batch in MANIFEST:
        path = os.path.join(data_dir, rel)
        if not os.path.exists(path):
            print(f"[warn] missing: {rel}", file=sys.stderr)
            continue
        proteins, peptides = parse_pd_export(path)
        meta = {"cell_line": cell, "treatment": treat, "batch": batch,
                "source_file": os.path.basename(rel)}

        # ---- Spindly PTM sites ----
        spindly_peps = [p for p in peptides
                        if SPINDLY_ACC in p.master_accessions
                        or SPINDLY_ACC in p.protein_acc]
        file_sites = []
        for p in spindly_peps:
            if high_conf_only and (p.confidence or "").lower() != "high":
                continue
            if p.xcorr is not None and p.xcorr < min_xcorr:
                continue
            file_sites.extend(extract_sites(p, ptm, SPINDLY_SEQ, meta))
        all_sites.extend(file_sites)

        # ---- interactome (protein-level hits) ----
        n_prot = 0
        for pr in proteins:
            acc = (pr[PROT_COLS["accession"]] or "").strip()
            if not acc:
                continue
            n_prot += 1
            interactome.append(ProteinHit(
                accession=acc,
                description=(pr[PROT_COLS["description"]] or "").strip(),
                coverage=_num(pr[PROT_COLS["coverage"]]),
                n_peptides=_int(pr[PROT_COLS["n_peptides"]]),
                score=_num(pr[PROT_COLS["score"]]),
                cell_line=cell, treatment=treat, batch=batch,
                source_file=os.path.basename(rel)))

        run_summary.append(dict(
            file=os.path.basename(rel), cell_line=cell, treatment=treat,
            ptm=ptm, batch=batch, n_proteins=n_prot,
            n_spindly_peptides=len(spindly_peps),
            n_localised_sites=sum(1 for s in file_sites if s.localised),
            n_ambiguous_sites=sum(1 for s in file_sites if not s.localised)))
        print(f"[ok] {rel:70s} proteins={n_prot:5d} "
              f"spindly_pep={len(spindly_peps):3d} sites={len(file_sites):3d}")

    _write_outputs(out_dir, all_sites, interactome, run_summary)
    return all_sites, interactome, run_summary


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _int(v):
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# 6. Output tables & summary matrix
# ---------------------------------------------------------------------------
def _write_outputs(out_dir, sites, interactome, run_summary):
    import csv

    # 6.1 raw site table
    with open(os.path.join(out_dir, "spindly_sites_all.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ptm", "residue", "position", "site_id", "localised",
                    "peptide", "confidence", "xcorr", "cell_line",
                    "treatment", "batch", "source_file"])
        for s in sorted(sites, key=lambda x: (x.ptm, x.position)):
            site_id = f"{s.residue}{s.position}" if s.localised else f"{s.residue}~{s.position}"
            w.writerow([s.ptm, s.residue, s.position, site_id, s.localised,
                        s.peptide, s.confidence, s.xcorr, s.cell_line,
                        s.treatment, s.batch, s.source_file])

    # 6.2 localised occupancy matrices  (residue x cell_line/treatment)
    phospho = [s for s in sites if s.ptm == "phospho" and s.localised]
    _occupancy_matrix(out_dir, phospho, "phospho_site_matrix.csv")
    ubi = [s for s in sites if s.ptm == "ubi" and s.localised]
    _occupancy_matrix(out_dir, ubi, "ubi_site_matrix.csv")

    # 6.3 interactome table
    with open(os.path.join(out_dir, "interactome_all_runs.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "description", "coverage_pct", "n_peptides",
                    "score", "cell_line", "treatment", "batch", "source_file"])
        for h in interactome:
            w.writerow([h.accession, h.description, h.coverage, h.n_peptides,
                        h.score, h.cell_line, h.treatment, h.batch, h.source_file])

    # 6.4 run summary
    with open(os.path.join(out_dir, "run_summary.csv"), "w", newline="") as fh:
        cols = ["file", "cell_line", "treatment", "ptm", "batch",
                "n_proteins", "n_spindly_peptides",
                "n_localised_sites", "n_ambiguous_sites"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for row in run_summary:
            w.writerow(row)

    print(f"\n[write] tables -> {out_dir}")


def _occupancy_matrix(out_dir, sites, fname):
    import csv
    conditions = sorted({(s.cell_line, s.treatment) for s in sites})
    residues = sorted({(s.position, s.residue) for s in sites})
    idx = {(pos, aa): {} for (pos, aa) in residues}
    for s in sites:
        idx[(s.position, s.residue)].setdefault((s.cell_line, s.treatment), 0)
        idx[(s.position, s.residue)][(s.cell_line, s.treatment)] += 1
    with open(os.path.join(out_dir, fname), "w", newline="") as fh:
        w = csv.writer(fh)
        header = ["site"] + [f"{c}|{t}" for (c, t) in conditions]
        w.writerow(header)
        for (pos, aa) in residues:
            row = [f"{aa}{pos}"]
            for cond in conditions:
                row.append(idx[(pos, aa)].get(cond, 0))
            w.writerow(row)


# ---------------------------------------------------------------------------
# 7. Figures
# ---------------------------------------------------------------------------
def make_figures(sites, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    phospho = [s for s in sites if s.ptm == "phospho" and s.localised]
    if not phospho:
        print("[fig] no localised phospho sites; skipping figure")
        return

    conditions = sorted({(s.cell_line, s.treatment) for s in phospho})
    residues = sorted({(s.position, s.residue) for s in phospho})
    labels = [f"{aa}{pos}" for (pos, aa) in residues]

    mat = np.zeros((len(residues), len(conditions)), dtype=int)
    for s in phospho:
        i = residues.index((s.position, s.residue))
        j = conditions.index((s.cell_line, s.treatment))
        mat[i, j] += 1

    # 7.1 Detection heatmap (PSM counts)
    fig, ax = plt.subplots(figsize=(1.4 * len(conditions) + 3, 0.5 * len(residues) + 2))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels([f"{c}\n{t.replace('_',' + ')}" for (c, t) in conditions], fontsize=8)
    ax.set_yticks(range(len(residues)))
    ax.set_yticklabels(labels, fontsize=9)
    for i in range(len(residues)):
        for j in range(len(conditions)):
            v = mat[i, j]
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=8,
                        color="black" if v < mat.max() * 0.6 else "white")
    ax.set_title("Spindly (Q96EA4) localised phosphosites\n(PSM count per condition)", fontsize=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("PSMs", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig_phospho_heatmap.png"), dpi=300)
    fig.savefig(os.path.join(out_dir, "fig_phospho_heatmap.pdf"))
    plt.close(fig)

    # 7.2 Linear map of sites along the protein
    fig, ax = plt.subplots(figsize=(10, 2.4))
    ax.plot([1, 605], [0, 0], color="0.6", lw=6, solid_capstyle="round")
    for (a, b, col, name) in [(60, 260, "#9ecae1", "CC1"),
                              (290, 440, "#9ecae1", "CC2"),
                              (440, 520, "#fdae6b", "Spindly-box region")]:
        ax.add_patch(plt.Rectangle((a, -0.12), b - a, 0.24, color=col, alpha=0.7, zorder=1))
        ax.text((a + b) / 2, -0.32, name, ha="center", va="top", fontsize=7)
    seen = {}
    for (pos, aa) in residues:
        seen[pos] = aa
    for pos, aa in sorted(seen.items()):
        ax.plot([pos, pos], [0.06, 0.55], color="#c0392b", lw=1.2, zorder=3)
        ax.plot(pos, 0.6, "o", color="#c0392b", ms=5, zorder=3)
        ax.text(pos, 0.7, f"{aa}{pos}", ha="center", va="bottom", fontsize=8, rotation=90)
    ax.set_xlim(-10, 615)
    ax.set_ylim(-0.6, 1.4)
    ax.set_yticks([])
    ax.set_xlabel("Residue number (Q96EA4, 605 aa)")
    ax.set_title("Phosphosites mapped on human Spindly", fontsize=10)
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig_phospho_linear_map.png"), dpi=300)
    fig.savefig(os.path.join(out_dir, "fig_phospho_linear_map.pdf"))
    plt.close(fig)
    print(f"[fig] figures -> {out_dir}")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", required=True, help="folder with the three experiment sub-folders")
    ap.add_argument("--out-dir", default="results", help="output folder")
    ap.add_argument("--min-xcorr", type=float, default=0.0, help="minimum Sequest HT XCorr for a PSM")
    ap.add_argument("--high-conf-only", action="store_true", help="keep only 'High' confidence PSMs")
    args = ap.parse_args()

    print(f"# Spindly sequence Q96EA4  len=605  md5={SEQ_MD5}")
    sites, inter, summ = run(args.data_dir, args.out_dir,
                             min_xcorr=args.min_xcorr,
                             high_conf_only=args.high_conf_only)
    make_figures(sites, args.out_dir)

    loc = sorted({(s.position, s.residue) for s in sites
                  if s.ptm == "phospho" and s.localised})
    print("\n# Localised Spindly phosphosites (union across all runs):")
    print("  " + ", ".join(f"{aa}{pos}" for pos, aa in loc))
    ubi = sorted({(s.position, s.residue) for s in sites
                  if s.ptm == "ubi" and s.localised})
    print("# Localised Spindly ubiquitination (GG) sites (union):")
    print("  " + ", ".join(f"{aa}{pos}" for pos, aa in ubi))


if __name__ == "__main__":
    main()
