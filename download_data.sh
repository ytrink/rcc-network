#!/usr/bin/env bash
# Download the processed data from Zenodo and unpack it into data/ (run from the repository root).
# Each zip holds its files under data/..., so unzipping at the repository root puts them in place.
set -euo pipefail

RECORD="23648073"   # Zenodo record of version v1 (concept DOI 10.5281/zenodo.23648072 always points to the latest version)
BASE="https://zenodo.org/records/${RECORD}/files"

ZIPS="scplus_grn.zip li_metadata.zip li_mac_mono_object.zip zvirblyte_endothelial.zip monocle_umaps.zip
      perturbation_mac_mono.zip perturbation_CD8.zip"

for z in $ZIPS; do
  echo "Downloading $z"
  curl -fL --retry 3 -o "$z" "${BASE}/${z}?download=1"
  unzip -o -q "$z"
  rm -f "$z"
done
rm -f $ZIPS   # v1 of perturbation_mac_mono.zip also contains copies of two other zips; remove them
echo "Done. Files are in data/."
