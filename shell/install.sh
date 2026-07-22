#!/bin/bash

# install.sh - Installs REVEAL and its dependencies into the active conda environment

set -euo pipefail

WD="$(dirname "$(dirname "$(realpath "$0")")")"

# Check conda is available
if ! command -v conda &> /dev/null; then
    echo "Error: conda not found. Please install conda first."
    exit 1
fi

# Check we're inside a conda environment
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: No conda environment is active. Please activate one first:"
    echo "  conda activate <your-env>"
    exit 1
fi

echo "Installing REVEAL into: $CONDA_PREFIX"
echo ""

# Install dependencies
conda install -y \
    -c bioconda \
    -c conda-forge \
    --no-channel-priority \
    "conda-forge::python>=3.10" \
    bioconda::pysam \
    conda-forge::numpy \
    conda-forge::pandas \
    "conda-forge::r-base>=4.0" \
    "conda-forge::r-tidyverse>=2.0.0"

# Remove existing installation if present
if [[ -d "$CONDA_PREFIX/lib/reveal" ]]; then
    echo "Existing REVEAL installation found, replacing..."
    rm -rf "$CONDA_PREFIX/lib/reveal"
fi
if [[ -f "$CONDA_PREFIX/bin/REVEAL" ]]; then
    rm -f "$CONDA_PREFIX/bin/REVEAL"
fi

# Copy scripts into the conda environment
mkdir -p "$CONDA_PREFIX/lib/reveal"
cp "$WD/src/"*.py "$CONDA_PREFIX/lib/reveal/"
cp "$WD/src/"*.R  "$CONDA_PREFIX/lib/reveal/"

# Create the REVEAL entry point
cat > "$CONDA_PREFIX/bin/REVEAL" << 'EOF'
#!/bin/bash
exec python "$(dirname "$(dirname "$(realpath "$0")")")/lib/reveal/reveal.py" "$@"
EOF
chmod +x "$CONDA_PREFIX/bin/REVEAL"

echo ""
echo "Done! You can now run:"
echo "  REVEAL --help"
