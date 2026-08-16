#!/bin/bash
# PerturbME Environment Setup Script
# Run this when conda network connectivity is stable

echo "🚀 Setting up PerturbME conda environment..."

# Load Miniforge3 module
echo "Loading Miniforge3 module..."
ml load Miniforge3/24.1.2-0

# Create conda environment
echo "Creating conda environment 'perturbme' with Python 3.10..."
conda create -n perturbme python=3.10 -y

# Activate environment
echo "Activating environment..."
source activate perturbme

# Install core scientific packages via conda
echo "Installing core packages via conda..."
conda install -c conda-forge numpy pandas scipy matplotlib seaborn scikit-learn statsmodels h5py tqdm -y

# Install single-cell packages via pip
echo "Installing single-cell analysis packages..."
pip install scanpy[leiden]>=1.8.0
pip install anndata>=0.8.0
pip install scvi-tools>=0.14.0

# Install specialized packages
echo "Installing specialized packages..."
pip install harmony-pytorch
pip install adjustText
pip install gseapy
pip install fastcluster
# pip install spams

echo "✅ Environment setup complete!"
echo ""
echo "To activate the environment, run:"
echo "  ml load Miniforge3/24.1.2-0"
echo "  conda activate perturbme"
echo ""
echo "To test the environment, run:"
echo "  python -c \"import scanpy, anndata, scvi; print('Environment ready!')\""
