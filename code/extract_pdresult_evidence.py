#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
extract_pdresult_evidence.py
============================
Regenerates the PSM-level evidence reported in the Spindly (SPDL1) Data
Descriptor directly from the Proteome Discoverer result files (.pdResult)
deposited in PRIDE (PXD085140).

A .pdResult file is an SQLite database. This script reads, for every
targeted-search result file:
  * WorkflowInputFiles  -> raw file name, acquisition date, instrument
  * Workflows           -> search parameters (database, tolerances, modifications)
  * TargetPsms          -> every peptide-spectrum match with its confidence
                           class (MatchConfidence) and localisation flag
                           (PsmAmbiguity)

Outputs (in --out-dir):
  spindly_phosphosite_evidence_pdResult.csv  Medium/High phospho-PSMs on Spindly
  table3_high_unambiguous_counts.csv         High + unambiguous PSMs per site and run
  y508h_variant_counts.csv                   High PSMs covering residue 508 (Y vs H)
  run_metadata.json                          acquisition metadata and search parameters

Usage:
  python extract_pdresult_evidence.py --pride-dir PATH/TO/PXD085140 --out-dir results

The PRIDE folder may be flat or organised in sub-folders; files are found by
name pattern (SPDL1_R<n>_<date>_targeted_*.pdResult).

Dependencies: Python >= 3.9 (standard library only).
License: MIT
"""
import argparse, collections, csv, glob, json, os, re, sqlite3
import xml.etree.ElementTree as ET

# Canonical UniProt Q96EA4 sequence (605 aa); the H508 form differs only at position 508.
Q96EA4 = ("MEADIITNLRCRLKEAEEERLKAAQYGLQLVESQNELQNQLDKCRNEMMTMTESYEQEKYTLQREVELKSRMLESLSCECEAIKQQQKMHLEKLEEQLSRSHGQEVNELKTKIEKLKVELDEARLSEKQLKHQVDHQKELLSCKSEELRVMSERVQESMSSEMLALQIELTEMESMKTTLKEEVNELQYRQEQLELLITNLMRQVDRLKEEKEEREKEAVSYYNALEKARVANQDLQVQLDQALQQALDPNSKGNSLFAEVEDRRAAMERQLISMKVKYQSLKKQNVFNREQMQRMKLQIATLLQMKGSQTEFEQQERLLAMLEQKNGEIKHLLGEIRNLEKFKNLYDSMESKPSVDSGTLEDNTYYTDLLQMKLDNLNKEIESTKGELSIQRMKALFESQRALDIERKLFANERCLQLSESENMKLRAKLDELKLKYEPEETVEVPVLKKRREVLPVDITTAKDACVNNSALGGEVYRLPPQKEETQSCPNSLEDNNLQLEKSVSIYTPVVSLSPHKNLPVDMQLKKEKKCVKLIGVPADAEALSERSGNTPNSPRLAAESKLQTEVKEGKETSSKLEKETCKKLHPILYVSSKSTPETQCPQQ")
assert len(Q96EA4) == 605
Q_H508 = Q96EA4[:507] + "H" + Q96EA4[508:]
CONF = {1: "Low", 2: "Medium", 3: "High"}
UNAMBIGUOUS = 5          # PsmAmbiguity value used by Proteome Discoverer for unambiguous PSMs


def search_parameters(con):
    xml = con.execute("select WorkflowXML from Workflows where WorkflowType='Processing'").fetchone()[0]
    p = {"modifications": []}
    for node in ET.fromstring(xml).iter():
        if not node.tag.endswith("ProcessingNodeParameter"):
            continue
        name, val = node.get("FriendlyName") or "", (node.text or "").strip()
        if name in ("Protein Database", "Precursor Mass Tolerance", "Fragment Mass Tolerance",
                    "Max. Missed Cleavage Sites", "Search Against Decoy Database"):
            p[name] = val
        if "Modification" in name and 'Name="' in val:
            mod = re.search(r'Name="([^"]+)"', val).group(1)
            aa = re.search(r'AminoAcids="([^"]*)"', val).group(1)
            p["modifications"].append(f"{mod} ({aa})")
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pride-dir", required=True)
    ap.add_argument("--out-dir", default="results")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(a.pride_dir, "**", "SPDL1_R*_targeted_*.pdResult"), recursive=True))
    if not files:
        raise SystemExit("No targeted .pdResult files found under --pride-dir")

    rows, meta = [], {}
    high_unamb = collections.defaultdict(collections.Counter)
    variant = collections.defaultdict(collections.Counter)
    for f in files:
        run = re.search(r"SPDL1_(R\d+)_", os.path.basename(f)).group(1)
        con = sqlite3.connect(f)
        raw = con.execute("select FileName, CreationDate, InstrumentDescription from WorkflowInputFiles "
                          "where FileType='XcaliburRawfile'").fetchone()
        meta[run] = {"pdResult": os.path.basename(f), "raw_file_at_acquisition": raw[0].split("\\")[-1],
                     "acquired": raw[1], "instrument": raw[2], "search": search_parameters(con)}
        q = ("select Sequence, Modifications, PsmAmbiguity, MatchConfidence, XCorr, Charge, FirstScan, "
             "RetentionTime, DeltaMassInPPM from TargetPsms")
        for seq, mods, amb, mc, xc, z, scan, rt, ppm in con.execute(q):
            off, var = Q96EA4.find(seq), "Y508 (canonical)"
            if off < 0:
                off, var = Q_H508.find(seq), "H508 (rs3797713)"
                if off < 0:
                    continue
            if mc == 3 and off < 508 <= off + len(seq):
                variant[run][var] += 1
            if mc < 2 or "Phospho" not in mods:
                continue
            for aa, p in re.findall(r"([A-Z])(\d+)\(Phospho\)", mods):
                pos = off + int(p)
                site = f"{Q96EA4[pos - 1]}{pos}"
                loc = "unambiguous" if amb == UNAMBIGUOUS else "ambiguous"
                if mc == 3 and loc == "unambiguous":
                    high_unamb[site][run] += 1
                rows.append(dict(run_id=run, site=site, residue=Q96EA4[pos - 1], position=pos,
                                 confidence=CONF[mc], localisation=loc, peptide=seq, peptide_variant=var,
                                 modifications=mods, xcorr=round(xc, 3), charge=z, first_scan=scan,
                                 rt_min=round(rt, 2), delta_ppm=round(ppm, 2), source_file=os.path.basename(f)))

    runs = sorted(meta, key=lambda r: int(r[1:]))
    with open(os.path.join(a.out_dir, "spindly_phosphosite_evidence_pdResult.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
        w.writerows(sorted(rows, key=lambda d: (d["position"], d["run_id"])))
    with open(os.path.join(a.out_dir, "table3_high_unambiguous_counts.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["site"] + runs)
        for s in sorted(high_unamb, key=lambda s: int(s[1:])):
            w.writerow([s] + [high_unamb[s].get(r, 0) for r in runs])
    with open(os.path.join(a.out_dir, "y508h_variant_counts.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["run_id", "high_psms_H508", "high_psms_Y508"])
        for r in runs:
            w.writerow([r, variant[r].get("H508 (rs3797713)", 0), variant[r].get("Y508 (canonical)", 0)])
    with open(os.path.join(a.out_dir, "run_metadata.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    print(f"{len(files)} result files, {len(rows)} Medium/High phospho-PSM rows written to {a.out_dir}")


if __name__ == "__main__":
    main()
