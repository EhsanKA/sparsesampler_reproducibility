#!/bin/bash

# Download script for sparsesampler reproducibility datasets

set -e  # Exit on any error

echo "Downloading datasets for sparsesampler reproducibility..."

# Create directories if they don't exist
mkdir -p lcmv
mkdir -p mcc

# Download LCMV data
echo "Downloading LCMV data..."
cd lcmv
wget -O 2024-02-27_LCMV_all_cells.csv "https://zenodo.org/records/10694407/files/2024-02-27_LCMV_all_cells.csv?download=1"
cd ..

# Download MCC data
echo "Downloading MCC data..."
cd mcc
wget -O adata.h5ad "https://datasets.cellxgene.cziscience.com/b8cfe635-1cf6-4a6d-8bc0-5059477b9a8c.h5ad"
cd ..

echo "Download completed successfully!"
echo "LCMV data saved to: lcmv/2024-02-27_LCMV_all_cells.csv"
echo "MCC data saved to: mcc/adata.h5ad"
