#!/bin/bash
# install_deerfold.sh

# Install Mamba
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh" -O mamba_installer.sh
bash mamba_installer.sh -b -p $HOME/mambaforge
rm mamba_installer.sh
source $HOME/mambaforge/etc/profile.d/conda.sh

# Clone DEERFold
git clone https://github.com/CAAIPD/DEERFold.git
cd DEERFold

# Create and activate environment
mamba env create -f deerfold_environment.yml -p ./env/deerfold_venv
mamba activate ./env/deerfold_venv

# Install PyTorch with CUDA support (adjust based on system)
pip install torch==1.12.1+cu113 torchvision==0.13.1+cu113 torchaudio==0.12.1 --extra-index-url https://download.pytorch.org/whl/cu113

# Install MDAnalysis and chiLife
pip install setuptools==60.3.0 gsd==2.4.2 mdanalysis==2.0.0
git clone https://github.com/StollLab/chiLife.git
cd chiLife
pip install -e .
cd ..

echo "Installation complete! Activate the environment with: mamba activate ./env/deerfold_venv"