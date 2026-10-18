#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_manuscript_figures.py
==========================
Builds the display items for the Spindly PTM Data Descriptor:

  Figure 1  Curated phosphosite occurrence matrix (residue x condition).
  Figure 2  Linear domain architecture of Spindly with mapped phosphosites.
  Figure 3  AlphaFold model of Spindly (UniProt Q96EA4, model v6):
              (A) C-alpha backbone coloured by per-residue pLDDT with the
                  seven phosphosites shown as spheres;
              (B) per-residue pLDDT profile with the phosphosites marked.

The curated phosphosite calls come from manual spectral validation (the
authoritative result of the study). The AlphaFold C-alpha coordinates and
pLDDT values are embedded (extracted from AF-Q96EA4-F1-model_v6) so the
figure is reproducible offline; the companion PyMOL script
`spindly_alphafold.pml` regenerates a publication-grade 3D render with
the full atomic model in PyMOL.

License: MIT
"""
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection

# ---------------------------------------------------------------------------
# Curated, manually validated phosphosites (authoritative study result)
#   from the spectral inspection of the LC-MS/MS runs.
# ---------------------------------------------------------------------------
CONDITIONS = ["RPE-1 NZ", "RPE-1 NZ+OK", "HCT116 NZ", "HCT116 NZ+OK"]
PHOSPHOSITES = ["S471", "T509", "S513", "S515", "H517", "K518", "S555"]
# 1 = phosphosite detected & localised in that condition
CURATED = {
    "S471": [0, 1, 0, 1],
    "T509": [1, 0, 0, 1],
    "S513": [1, 1, 1, 1],
    "S515": [1, 1, 1, 1],
    "H517": [0, 1, 0, 0],
    "K518": [0, 0, 1, 0],
    "S555": [0, 0, 0, 1],
}

# ---------------------------------------------------------------------------
# AlphaFold C-alpha trace (AF-Q96EA4-F1-model_v6): [resi, x, y, z, pLDDT]
# Sampled every 4th residue plus the seven phosphosites (integer-rounded).
# ---------------------------------------------------------------------------
CA = [
[1,79,21,-131,57],[5,75,19,-125,68],[9,71,19,-120,80],[13,68,20,-115,85],[17,64,21,-111,90],
[21,60,23,-107,93],[25,55,23,-104,94],[29,50,23,-100,95],[33,46,22,-96,95],[37,42,20,-91,97],
[41,40,19,-85,97],[45,38,19,-79,98],[49,37,19,-73,98],[53,35,20,-68,98],[57,32,22,-62,98],
[61,29,23,-57,98],[65,25,23,-52,98],[69,21,22,-48,98],[73,17,21,-43,98],[77,15,19,-38,98],
[81,12,18,-32,97],[85,11,17,-26,97],[89,9,18,-21,95],[93,6,19,-15,96],[97,2,20,-10,95],
[101,-2,20,-6,92],[105,-7,20,-2,94],[109,-12,19,2,94],[113,-16,17,6,94],[117,-20,14,10,94],
[121,-23,13,15,94],[125,-27,13,20,95],[129,-30,13,25,96],[133,-34,14,30,95],[137,-38,16,34,93],
[141,-43,17,38,93],[145,-48,18,41,92],[149,-52,17,45,90],[153,-57,16,49,81],[157,-58,12,57,69],
[161,-50,6,63,85],[165,-46,5,59,90],[169,-42,4,54,92],[173,-40,3,48,93],[177,-37,4,43,93],
[181,-35,5,37,93],[185,-31,7,32,94],[189,-27,8,28,96],[193,-22,8,24,96],[197,-17,7,20,95],
[201,-13,6,16,94],[205,-9,4,12,94],
[209,-6,2,7,94],[213,-2,2,2,93],[217,1,3,-3,92],[221,6,4,-7,89],[225,11,5,-10,94],
[229,17,5,-12,94],[233,23,5,-14,94],[237,28,4,-17,95],[241,33,2,-20,90],[245,38,2,-24,82],
[249,42,1,-28,58],[253,38,-4,-34,53],[257,39,-16,-32,59],[261,34,-18,-27,81],[265,30,-18,-22,88],
[269,26,-17,-17,91],[273,22,-16,-14,92],[277,17,-15,-10,94],[281,12,-15,-7,89],[285,7,-15,-4,92],
[289,2,-17,0,89],[293,-2,-18,4,87],[297,-6,-19,9,80],[301,-10,-20,14,75],[305,-13,-19,18,75],
[309,-19,-18,23,64],[313,-25,-15,25,74],[317,-31,-15,26,84],[321,-37,-15,27,87],[325,-43,-16,30,91],
[329,-48,-16,34,92],[333,-52,-15,39,93],[337,-56,-13,43,88],[341,-60,-11,47,85],[345,-64,-9,51,79],
[349,-69,-8,54,64],[353,-74,-10,62,41],[357,-82,-19,68,35],[361,-89,-30,68,42],[365,-86,-35,61,70],
[369,-81,-36,57,86],[373,-77,-37,53,92],[377,-72,-36,49,94],[381,-67,-34,46,95],[385,-62,-33,43,94],
[389,-56,-31,42,94],[393,-50,-30,41,95],[397,-44,-30,39,92],[401,-39,-30,37,92],[405,-33,-31,33,91],
[409,-28,-32,29,92],[413,-24,-31,25,89],
[417,-20,-30,21,87],[421,-15,-29,17,86],[425,-10,-27,14,82],[429,-4,-26,12,80],[433,1,-26,10,80],
[437,7,-28,8,72],[441,9,-30,-2,49],[445,11,-36,-11,41],[449,18,-35,-20,47],[453,30,-34,-22,51],
[457,41,-37,-27,45],[461,50,-44,-23,40],[465,50,-52,-12,40],[469,41,-57,-3,37],[471,34,-57,-1,39],
[473,28,-55,-2,42],[477,16,-52,-8,36],[481,4,-48,-14,38],[485,-10,-44,-15,38],[489,-21,-41,-7,30],
[493,-26,-46,4,29],[497,-23,-54,14,40],[501,-24,-52,27,33],[505,-25,-43,36,29],[509,-19,-34,40,33],
[513,-15,-23,42,30],[515,-14,-17,39,33],[517,-11,-11,39,34],[518,-8,-9,39,38],[521,2,-9,36,36],
[525,12,-17,33,32],[529,18,-26,26,40],[533,29,-22,18,51],[537,36,-11,16,44],[541,35,3,18,41],
[545,28,12,24,41],[549,15,15,27,34],[553,9,24,20,38],[555,13,27,15,34],[557,17,30,11,34],
[561,24,34,-1,40],[565,24,41,-11,33],[569,20,47,-22,42],[573,12,43,-32,35],[577,2,42,-40,44],
[581,-9,40,-33,43],[585,-10,41,-19,43],[589,-9,40,-7,60],[593,-11,40,5,43],[597,-16,39,17,48],
[601,-23,36,29,52],[605,-29,29,38,44],
]
SITES = {471: "S471", 509: "T509", 513: "S513", 515: "S515",
         517: "H517", 518: "K518", 555: "S555"}


def plddt_color(v):
    """AlphaFold pLDDT colour bands."""
    if v >= 90:  return "#0053D6"   # very high
    if v >= 70:  return "#65CBF3"   # confident
    if v >= 50:  return "#FFDB13"   # low
    return "#FF7D45"                 # very low (disordered)


# ---------------------------------------------------------------------------
def fig1_matrix(out_dir):
    mat = np.array([CURATED[s] for s in PHOSPHOSITES])
    fig, ax = plt.subplots(figsize=(5.4, 4.2))
    ax.imshow(mat, cmap=matplotlib.colors.ListedColormap(["#f2f2f2", "#c0392b"]),
              vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(CONDITIONS)))
    ax.set_xticklabels(CONDITIONS, rotation=20, ha="right", fontsize=9)
    ax.set_yticks(range(len(PHOSPHOSITES)))
    ax.set_yticklabels(PHOSPHOSITES, fontsize=10)
    for i in range(len(PHOSPHOSITES)):
        for j in range(len(CONDITIONS)):
            if mat[i, j]:
                ax.text(j, i, "+", ha="center", va="center", color="white",
                        fontsize=13, fontweight="bold")
    ax.set_xticks(np.arange(-.5, len(CONDITIONS), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(PHOSPHOSITES), 1), minor=True)
    ax.grid(which="minor", color="white", lw=2)
    ax.tick_params(which="minor", length=0)
    ax.set_title("Curated Spindly phosphosites per condition", fontsize=11)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"Fig1_phosphosite_matrix.{ext}"), dpi=300)
    plt.close(fig)


def fig2_domain_map(out_dir):
    fig, ax = plt.subplots(figsize=(10, 2.6))
    ax.plot([1, 605], [0, 0], color="0.55", lw=7, solid_capstyle="round", zorder=1)
    regions = [(20, 250, "#9ecae1", "Coiled-coil 1"),
               (255, 440, "#9ecae1", "Coiled-coil 2"),
               (276, 306, "#fdae6b", "Spindly box"),
               (450, 605, "#f2c1c1", "Disordered C-terminus")]
    for a, b, col, name in regions:
        ax.add_patch(plt.Rectangle((a, -0.13), b - a, 0.26, color=col,
                                   alpha=0.85, zorder=2, ec="0.4", lw=0.5))
        ax.text((a + b) / 2, -0.34, name, ha="center", va="top", fontsize=7.5)
    # stagger stem heights so the 509-518 cluster does not overlap
    heights = [0.55, 0.95, 1.35, 0.75, 1.15]
    for i, (pos, lab) in enumerate(sorted(SITES.items())):
        h = heights[i % len(heights)]
        ax.plot([pos, pos], [0.08, h], color="#c0392b", lw=1.1, zorder=3)
        ax.plot(pos, h, "o", color="#c0392b", ms=5, zorder=3)
        ax.text(pos, h + 0.08, lab, ha="center", va="bottom", fontsize=8.5)
    ax.set_xlim(-15, 620)
    ax.set_ylim(-0.7, 1.9)
    ax.set_yticks([])
    ax.set_xlabel("Residue number (Q96EA4, 605 aa)")
    ax.set_title("Phosphosites mapped on the domain architecture of human Spindly",
                 fontsize=11)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"Fig2_domain_map.{ext}"), dpi=300)
    plt.close(fig)


def fig3_alphafold(out_dir):
    arr = np.array(CA, dtype=float)
    resi, xyz, plddt = arr[:, 0], arr[:, 1:4], arr[:, 4]

    fig = plt.figure(figsize=(11, 4.8))
    # -- Panel A: 3D backbone coloured by pLDDT --
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    pts = xyz.reshape(-1, 1, 3)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    seg_cols = [plddt_color((plddt[i] + plddt[i + 1]) / 2) for i in range(len(plddt) - 1)]
    lc = Line3DCollection(segs, colors=seg_cols, linewidths=3.2)
    ax.add_collection3d(lc)
    # highlight phosphosites as spheres (residue identities are given in
    # Fig. 2 and panel B; labels are omitted here because the seven sites
    # co-localise in the compact disordered C-terminus)
    hx = []
    for pos in SITES:
        k = np.where(resi == pos)[0]
        if len(k):
            hx.append(xyz[k[0]])
    hx = np.array(hx)
    ax.scatter(hx[:, 0], hx[:, 1], hx[:, 2], color="#c0392b", edgecolor="black",
               s=70, depthshade=False, zorder=6, label="Phosphosites (7)")
    # single leader annotation for the cluster
    cx = hx.mean(axis=0)
    ax.text(cx[0] + 6, cx[1] + 6, cx[2] + 20,
            "7 phosphosites\n(471–555)", fontsize=8, fontweight="bold",
            color="#c0392b", ha="left")
    _set_axes_equal(ax, xyz)
    ax.set_title("AlphaFold model (Q96EA4, v6)\nCα backbone coloured by pLDDT",
                 fontsize=10)
    ax.set_axis_off()
    # pLDDT colour legend
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], color="#0053D6", lw=4, label="Very high (≥90)"),
           Line2D([0], [0], color="#65CBF3", lw=4, label="Confident (70–90)"),
           Line2D([0], [0], color="#FFDB13", lw=4, label="Low (50–70)"),
           Line2D([0], [0], color="#FF7D45", lw=4, label="Very low (<50)"),
           Line2D([0], [0], marker="o", color="w", markerfacecolor="#c0392b",
                  markeredgecolor="black", markersize=8, label="Phosphosite")]
    ax.legend(handles=leg, loc="upper left", fontsize=6.5, frameon=False,
              bbox_to_anchor=(-0.05, 1.0))

    # -- Panel B: pLDDT profile --
    ax2 = fig.add_subplot(1, 2, 2)
    # reconstruct a per-residue-ish line from the sampled points
    order = np.argsort(resi)
    ax2.axhspan(90, 100, color="#0053D6", alpha=0.08)
    ax2.axhspan(70, 90, color="#65CBF3", alpha=0.10)
    ax2.axhspan(50, 70, color="#FFDB13", alpha=0.12)
    ax2.axhspan(0, 50, color="#FF7D45", alpha=0.12)
    ax2.plot(resi[order], plddt[order], color="0.2", lw=1.4)
    # staggered leader annotations (the 509-518 sites are nearly coincident)
    offsets = [(-18, 20), (-6, 34), (2, 20), (8, 46), (14, 32), (20, 12), (26, 26)]
    for i, (pos, lab) in enumerate(sorted(SITES.items())):
        k = np.where(resi == pos)[0]
        if len(k):
            yv = plddt[k[0]]
            dx, dy = offsets[i % len(offsets)]
            ax2.plot(pos, yv, "o", color="#c0392b", ms=6, zorder=5)
            ax2.annotate(lab, (pos, yv), textcoords="offset points",
                         xytext=(dx, dy), ha="center", fontsize=7.5, color="#c0392b",
                         arrowprops=dict(arrowstyle="-", color="#c0392b", lw=0.5))
    ax2.set_xlim(0, 610)
    ax2.set_ylim(0, 100)
    ax2.set_xlabel("Residue number")
    ax2.set_ylabel("pLDDT")
    ax2.set_title("Per-residue model confidence", fontsize=10)
    ax2.text(0.02, 0.06, "All seven phosphosites lie in the low-confidence\n"
             "(pLDDT < 40) disordered C-terminal region",
             transform=ax2.transAxes, fontsize=7.5, style="italic", color="0.25")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"Fig3_alphafold_structure.{ext}"), dpi=300)
    plt.close(fig)


def _set_axes_equal(ax, xyz):
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    mid = [(x.max() + x.min()) / 2, (y.max() + y.min()) / 2, (z.max() + z.min()) / 2]
    r = max(x.max() - x.min(), y.max() - y.min(), z.max() - z.min()) / 2
    ax.set_xlim(mid[0] - r, mid[0] + r)
    ax.set_ylim(mid[1] - r, mid[1] + r)
    ax.set_zlim(mid[2] - r, mid[2] + r)
    ax.view_init(elev=18, azim=-60)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    fig1_matrix(args.out_dir)
    fig2_domain_map(args.out_dir)
    fig3_alphafold(args.out_dir)
    print("figures written to", args.out_dir)


if __name__ == "__main__":
    main()
