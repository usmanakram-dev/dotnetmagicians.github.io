# Publication-Ready Circular Phylogenetic Tree Generator

## Overview

This repository contains a Python script that generates publication-quality circular phylogenetic trees from Mycobacterium tuberculosis drug resistance mutation data.

## Features

- **Hamming Distance Calculation**: Computes pairwise distances from binary mutation matrices
- **UPGMA Tree Construction**: Builds phylogenetic trees using the Unweighted Pair Group Method with Arithmetic Mean
- **Circular Layout**: Beautiful circular tree visualization suitable for publications
- **Color-Coded Labels**: Mutations groups are distinguished by color:
  - **S450L** (Serine → Leucine at position 450): Blue (#0066CC)
  - **S531L** (Serine → Leucine at position 531): Red (#CC0000)
  - **H526Y** (Histidine → Tyrosine at position 526): Green (#009933)
  - **D516V** (Aspartic acid → Valine at position 516): Purple (#9933CC)
- **Branch Lengths**: Displays evolutionary distances
- **High-Resolution Output**: 2000×2000 pixels at 300 DPI

## Dataset

The script includes a hard-coded dataset of 15 Mycobacterium tuberculosis samples with 20 genomic positions representing drug resistance mutations in the rpoB gene (rifampicin resistance).

## Installation

### Prerequisites

- Python 3.7 or higher
- pip (Python package installer)

### Install Dependencies

```bash
pip install -r requirements.txt
```

Or install packages individually:

```bash
pip install numpy ete3 PyQt5
```

## Usage

### Basic Usage

```bash
python3 publication_phylogenetic_tree.py
```

### Headless Environments

For systems without a display (servers, CI/CD pipelines):

```bash
QT_QPA_PLATFORM=offscreen python3 publication_phylogenetic_tree.py
```

Or using xvfb:

```bash
xvfb-run -a python3 publication_phylogenetic_tree.py
```

### Making the Script Executable

```bash
chmod +x publication_phylogenetic_tree.py
./publication_phylogenetic_tree.py
```

## Output

The script generates a file named `circular_phylogenetic_tree.png` in the current directory.

**Output specifications:**
- Format: PNG
- Dimensions: 2000×2000 pixels
- DPI: 300 (publication quality)
- Color depth: 8-bit RGBA

## Script Structure

The script consists of several key functions:

1. **`compute_hamming_distance()`**: Calculates Hamming distance between binary sequences
2. **`build_distance_matrix()`**: Builds pairwise distance matrix
3. **`upgma()`**: Implements UPGMA clustering algorithm
4. **`style_tree()`**: Applies publication-quality styling
5. **`main()`**: Orchestrates the entire pipeline

## Example Output

The generated tree includes:
- 15 samples arranged in a circular layout
- Color-coded labels by mutation group
- Branch lengths indicating evolutionary distances
- Title and legend
- Colored spherical markers for each sample

## Customization

### Modifying the Dataset

Edit the `MUTATION_MATRIX` array in the script to use your own binary mutation data:

```python
MUTATION_MATRIX = np.array([
    [1, 0, 0, 1, ...],  # Sample 1
    [0, 1, 0, 0, ...],  # Sample 2
    # ... more samples
])
```

### Changing Colors

Modify the `GROUP_COLORS` dictionary:

```python
GROUP_COLORS = {
    "S450L": "#0066CC",  # Blue
    "S531L": "#CC0000",  # Red
    # ... other groups
}
```

### Adjusting Tree Style

Edit parameters in the `style_tree()` function:
- Font sizes: `fsize` parameter in `TextFace()`
- Circle sizes: `radius` parameter in `CircleFace()`
- Branch spacing: `branch_vertical_margin` in `TreeStyle()`
- Arc angles: `arc_start` and `arc_span` for partial circles

### Output Resolution

Modify the render parameters in `main()`:

```python
tree.render(output_file, w=3000, h=3000, units="px", tree_style=tree_style, dpi=600)
```

## Scientific Background

### UPGMA Algorithm

UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is a simple agglomerative hierarchical clustering method used to construct phylogenetic trees. It assumes a constant evolutionary rate (molecular clock hypothesis).

### Hamming Distance

The Hamming distance between two binary sequences is the number of positions at which the corresponding symbols differ. It's appropriate for binary mutation matrices where 1 indicates presence and 0 indicates absence of a mutation.

### rpoB Gene

The rpoB gene encodes the β-subunit of bacterial RNA polymerase. Mutations in this gene are the primary mechanism of rifampicin resistance in Mycobacterium tuberculosis.

## Troubleshooting

### Qt Platform Plugin Error

If you see: `Could not load the Qt platform plugin "xcb"`

**Solution 1**: Use offscreen rendering
```bash
QT_QPA_PLATFORM=offscreen python3 publication_phylogenetic_tree.py
```

**Solution 2**: Install Qt dependencies
```bash
sudo apt-get install libxcb-xinerama0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1
```

**Solution 3**: Use xvfb
```bash
sudo apt-get install xvfb
xvfb-run -a python3 publication_phylogenetic_tree.py
```

### Import Errors

Ensure all dependencies are installed:
```bash
pip install --upgrade numpy ete3 PyQt5
```

## References

- ETE Toolkit: http://etetoolkit.org/
- UPGMA: Sokal, R. R., & Michener, C. D. (1958)
- Rifampicin resistance mutations: Telenti et al. (1993), Lancet

## License

This script is provided as-is for academic and research purposes.

## Author

Generated for publication-ready phylogenetic tree visualization of Mycobacterium tuberculosis drug resistance mutations.
