#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
read_pdresult.py - extract acquisition and search metadata from a Thermo
Proteome Discoverer result file (.pdResult).

A .pdResult file is an SQLite database. The script reads:
  - WorkflowInputFiles : original .raw name and size, acquisition date,
    retention-time range (= method length), instrument model and firmware.
  - Workflows.WorkflowXML : all processing-node parameters (Sequest HT:
    FASTA, enzyme, tolerances, modifications) and consensus-node parameters
    (validation, confidence thresholds, filters).

Usage:
    python read_pdresult.py file.pdResult [...]  [--json out.json]

The file is opened read-only and is never modified.
"""

import argparse
import json
import os
import re
import sqlite3
import sys
import xml.etree.ElementTree as ET

# Workflow parameters reported in Methods
KEEP = [
    "SpectrumFileNames", "MassAnalyzerFilter", "MSOrderFilter",
    "ActivationTypeFilter", "MinimumPrecursorMass", "MaximumPrecursorMass",
    "FastaDatabase", "Enzyme", "MaxMissedCleavages", "MinimumPeptideLength",
    "MaximumPeptideLength", "PeptideTolerance", "FragmentTolerance",
    "MaxEqualModificationsPerPeptide", "MaxDifferentialModificationsPerPeptide",
    "DynamicModification1", "DynamicModification2", "DynamicModification3",
    "DynamicModification4", "DynamicModificationNProtein",
    "StaticModification1", "StaticModification2",
    "IonSerieA", "IonSerieB", "IonSerieC", "IonSerieX", "IonSerieY", "IonSerieZ",
    "UseFlankingIons", "MaxPeptidesReported",
    "CutOffFTCharge1High", "CutOffFTCharge2High", "CutOffFTCharge3High",
    "CutOffFTCharge4High", "CutOffFTCharge1Middle", "CutOffFTCharge2Middle",
    "CutOffFTCharge3Middle", "CutOffFTCharge4Middle",
    "DecoySearch", "ConcatenatedFdrCalculation",
    "MinModificationSiteProbability", "ChosenValidationMode",
    "TargetFPRHigh", "TargetFPRMiddle", "TargetPeptideFPRHigh",
    "TargetPeptideFPRMiddle", "ValidationBasedOn",
    "PeptideConfidenceFilter", "PeptideLengthFilter", "MaximumDeltaCn",
    "ApplyStrictParsimony", "HighConfidenceThreshold",
    "MediumConfidenceThreshold", "UseDecoyDatabase",
]

MOD_RE = re.compile(r'AminoAcids="([^"]*)"\s+Name="([^"]*)"')


def short_mod(xml_text):
    m = MOD_RE.search(xml_text or "")
    return f"{m.group(2)} ({m.group(1)})" if m else xml_text


def read_workflow_params(xml_text):
    out = []
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return out
    for node in root.iter():
        if not node.tag.endswith("WorkflowNode"):
            continue
        friendly = node.get("FriendlyName") or node.get("ProcessingNodeName")
        params = {}
        for p in node.iter():
            if p.tag.endswith("ProcessingNodeParameter"):
                name, val = p.get("Name"), (p.text or "").strip()
                if name in KEEP and val:
                    params[name] = short_mod(val) if "Modification" in name or name == "Enzyme" else val
        if params:
            out.append((friendly, params))
    return out


def read_pdresult(path):
    con = sqlite3.connect("file:" + path + "?mode=ro", uri=True)
    rec = {"pdresult": os.path.basename(path), "input_files": [], "workflows": []}

    for row in con.execute(
        "select FileName, FileSize, FileType, CreationDate, RetentionTimeRange, "
        "InstrumentDescription, SoftwareRevision from WorkflowInputFiles"
    ):
        rec["input_files"].append(dict(zip(
            ["file_name", "file_size", "file_type", "creation_date",
             "rt_range_min", "instrument", "firmware"], row)))

    for wid, wtype, sw, start, xml_text in con.execute(
        "select WorkflowID, WorkflowType, SoftwareVersion, WorkflowStartDate, WorkflowXML from Workflows"
    ):
        rec["workflows"].append({
            "workflow_id": wid, "type": wtype, "pd_version": sw,
            "start": start, "nodes": read_workflow_params(xml_text),
        })

    # Spectrum and protein counts, and engine version, from the processing log
    log = {}
    for (msg,) in con.execute("select Message from WorkflowMessages where Message is not null"):
        m = re.search(r"\(([\d,]+) spectra total\)", msg)
        if m:
            log["spectra_total"] = int(m.group(1).replace(",", ""))
        m = re.search(r"^Sent ([\d,]+) spectra from file", msg)
        if m:
            log["ms2_searched"] = int(m.group(1).replace(",", ""))
        m = re.search(r"^Filter (\d+)/(\d+) Proteins", msg)
        if m:
            log["proteins_retained"], log["proteins_total"] = int(m.group(1)), int(m.group(2))
        m = re.search(r"^Found (\d+) protein groups", msg)
        if m:
            log["protein_groups"] = int(m.group(1))
        m = re.search(r"^ISE \(([\d\.]+)", msg)
        if m:
            log["search_engine_build"] = "ISE " + m.group(1)
        if "no decoy search was performed" in msg.lower():
            log["decoy_note"] = " ".join(msg.split())
    rec["log"] = log

    # Parameters observed in the spectra submitted to the search
    obs = {}
    try:
        q = ("select min(RetentionTime), max(RetentionTime), min(Charge), max(Charge), "
             "round(min(MassOverCharge),1), round(max(MassOverCharge),1), "
             "max(IonInjectTime), count(*) from MSnSpectrumInfo")
        rt0, rt1, c0, c1, mz0, mz1, it, n = con.execute(q).fetchone()
        obs = {"rt_min_observed": rt0, "rt_max_observed": rt1,
               "charge_min_observed": c0, "charge_max_observed": c1,
               "mz_min_observed": mz0, "mz_max_observed": mz1,
               "max_ion_inject_time_ms_observed": it, "n_ms2": n}
        # MS2 per survey scan: lower bound of the DDA Top-N
        counts = {}
        for (master,) in con.execute("select MasterScanNumbers from MSnSpectrumInfo"):
            counts[master] = counts.get(master, 0) + 1
        if counts:
            obs["max_ms2_per_survey_observed"] = max(counts.values())
            obs["n_survey_scans_with_ms2"] = len(counts)
    except sqlite3.Error:
        pass
    rec["observed"] = obs
    con.close()
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+")
    ap.add_argument("--json")
    args = ap.parse_args()

    report = []
    for p in args.paths:
        rec = read_pdresult(p)
        report.append(rec)
        print("=" * 78)
        print(rec["pdresult"])
        for f in rec["input_files"]:
            print(f"  [{f['file_type']}] {f['file_name']}")
            print(f"      size     : {f['file_size']:,} B" if f["file_size"] else "")
            print(f"      acquired : {f['creation_date']}")
            print(f"      RT range : {f['rt_range_min']} min")
            print(f"      instrument: {f['instrument']}  fw {f['firmware']}")
        if rec.get("log"):
            print("  -- processing log")
            for k, v in rec["log"].items():
                print(f"     {k} = {v}")
        if rec.get("observed"):
            print("  -- observed in searched spectra")
            for k, v in rec["observed"].items():
                print(f"     {k} = {v}")
        for wf in rec["workflows"]:
            print(f"  -- workflow {wf['type']} (PD {wf['pd_version']})")
            for friendly, params in wf["nodes"]:
                print(f"     [{friendly}]")
                for k, v in params.items():
                    print(f"        {k} = {v}")
    if args.json:
        json.dump(report, open(args.json, "w", encoding="utf-8"),
                  ensure_ascii=False, indent=2)
        print(f"[info] JSON written to {args.json}", file=sys.stderr)


if __name__ == "__main__":
    main()
