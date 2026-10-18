# spindly_alphafold.pml
# ---------------------------------------------------------------------------
# Publication-grade rendering of the AlphaFold model of human Spindly
# (UniProt Q96EA4) with the seven experimentally identified phosphosites.
#
# Usage:
#   1) Download the model (internet required):
#        curl -O https://alphafold.ebi.ac.uk/files/AF-Q96EA4-F1-model_v6.pdb
#      (If v6 is superseded, get the current pdbUrl from
#       https://alphafold.ebi.ac.uk/api/prediction/Q96EA4 )
#   2) Run:  pymol spindly_alphafold.pml
#
# Produces spindly_phosphosites.png (300 dpi, ray-traced).
# ---------------------------------------------------------------------------

load AF-Q96EA4-F1-model_v6.pdb, spindly
bg_color white
hide everything
show cartoon, spindly
set cartoon_transparency, 0.0

# ---- colour the cartoon by AlphaFold pLDDT (stored in the B-factor column) --
spectrum b, blue_white_orange, spindly, minimum=30, maximum=100

# ---- highlight the seven phosphosites -------------------------------------
select phosphosites, resi 471+509+513+515+517+518+555
show spheres, phosphosites and name CA
color red, phosphosites
set sphere_scale, 0.8, phosphosites

# residue labels
label phosphosites and name CA, "%s%s" % (resn, resi)
set label_size, 16
set label_color, black

# ---- presentation ----------------------------------------------------------
set ray_shadows, 0
set antialias, 2
orient spindly
zoom spindly, 5
set ray_opaque_background, 0
ray 2400, 1600
png spindly_phosphosites.png, dpi=300
