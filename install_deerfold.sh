#!/bin/bash
# install_deerfold.sh

# Check if Miniforge3 is already installed
if [ -d "$HOME/miniforge3" ]; then
    echo "Miniforge3 already installed at $HOME/miniforge3. Skipping installation."
else
    # Download and install Miniforge3
    echo "Downloading Miniforge3 installer..."
    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O miniforge_installer.sh
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download Miniforge3 installer. Check your internet connection or the URL."
        exit 1
    fi
    bash miniforge_installer.sh -b -p "$HOME/miniforge3"
    rm miniforge_installer.sh
fi

# Source the conda setup script
source "$HOME/miniforge3/etc/profile.d/conda.sh"
export PATH="$HOME/miniforge3/bin:$PATH"

# Create and activate environment
echo "Creating Deerfold environment..."
mamba env create -f deerfold_environment.yml -p ./env/deerfold_venv
if [ $? -ne 0 ]; then
    echo "Error: Failed to create environment. Check deerfold_environment.yml."
    exit 1
fi

# Activate the environment
conda activate ./env/deerfold_venv
export PATH="$PWD/env/deerfold_venv/bin:$PATH"

# Install MDAnalysis and chiLife
pip install setuptools==60.3.0 gsd==2.4.2 mdanalysis==2.0.0
pip install chilife==1.0.2
pip install numpy==1.22.4 --force-reinstall 2>/dev/null

# Install TMscore
wget https://zhanggroup.org/TM-score/TMscoreLinux.zip
unzip TMscoreLinux.zip
rm TMscoreLinux.zip

echo "Installation complete! Activate the environment with: conda activate ./env/deerfold_venv"
