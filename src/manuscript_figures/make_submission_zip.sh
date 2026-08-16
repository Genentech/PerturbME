#!/usr/bin/env bash
# Stage 3: assemble and zip the RESIS package for upload to 1Press.
#
# Per "HOW to upload materials for Resis review":
#   - one single zipped file, labelled "yyyy-mm-dd <unix>"
#   - Images:   TIFF (preferred), AI, BMP, JPG, PNG, PSD, RAW
#   - Data:     Excel, TXT, CSV, PRISM
#   - Main text: PDF or DOCX
#   - 2 GB limit
# Anything outside those formats is pruned here so the archive contains only
# material RESIS can actually open.
set -euo pipefail

REPO=/gnet/is1/p01/shares/regevlab/hanchen/Pert_PG/perturb-me
SUB="$REPO/submission"
PKG="$SUB/20260721_Resis_package"
ZIP="$SUB/$(date +%F) $(whoami).zip"

# ---- supplementary files ship verbatim from doc/Supplementary Files -------------
rm -rf "$PKG/04_Supplementary"
cp -r "$REPO/doc/Supplementary Files" "$PKG/04_Supplementary"

# MAGeCK writes its Sweave report scaffolding next to the results: .Rnw (LaTeX +
# R chunks), the .R it extracts, and .log. Those are build artifacts, not data, and
# none of the three is an accepted format. .npy is a NumPy binary, likewise.
find "$PKG/04_Supplementary" -type f \( -name '*.Rnw' -o -name '*.R' -o -name '*.log' -o -name '*.npy' \) -delete

# PDF is accepted for the main text, not as an image format; both Extended Data
# figures ship as 600 dpi TIFF in 02_Figures.
rm -f "$PKG/04_Supplementary/Extended Data Figure "*.pdf
rm -rf "$PKG/02_Figures/vector_originals"

# Supplementary Code 1 ships as PDF; .ipynb is not an accepted format.
conda run -n perturbme python "$REPO/src/manuscript_figures/nb_to_pdf.py"
cp "$SUB/build/Supplementary Code 1 - ElasticNet beta module analysis.pdf" "$PKG/04_Supplementary/"
rm -f "$PKG/04_Supplementary/Supplementary Code 1 - ElasticNet beta module analysis.ipynb"

# [PerturbME_CoSci].pptx is an internal meeting deck ("Goal of the Meeting", "Next steps
# and Logistics"), not a supplementary file - it is not in the manuscript's supplementary
# list and should not go to the journal.
rm -f "$PKG/04_Supplementary/[PerturbME_CoSci].pptx"

# README as .txt (markdown is not an accepted format).
cp "$PKG/README.md" "$PKG/README.txt"

if [ -f "$PKG/01_Main_text/PUT_MANUSCRIPT_HERE.txt" ]; then
  echo "WARNING: 01_Main_text/ still holds the placeholder - the main text has not been added." >&2
fi

# ---- zip ------------------------------------------------------------------------
cd "$PKG"
rm -f "$ZIP" "$SUB/PerturbME_source_data.zip"
zip -r -q "$ZIP" 01_Main_text 02_Figures 03_Source_Data 04_Supplementary README.txt
zip -r -q "$SUB/PerturbME_source_data.zip" 03_Source_Data 04_Supplementary README.txt

echo "1Press upload:  $ZIP  ($(du -h "$ZIP" | cut -f1))"
echo "source data only: $SUB/PerturbME_source_data.zip  ($(du -h "$SUB/PerturbME_source_data.zip" | cut -f1))"
echo
echo "extensions in the archive:"
unzip -l "$ZIP" | awk '{print $NF}' | grep -oE '\.[A-Za-z0-9]+$' | sort | uniq -c | sort -rn
