#!/usr/bin/env python3
"""
Publication-ready Circular Phylogenetic Tree Generator
=======================================================

This script generates a publication-quality circular phylogenetic tree from
Mycobacterium tuberculosis drug resistance mutation data.

Dataset: Binary mutation matrix for 15 samples with rpoB gene mutations
Tree method: UPGMA (Unweighted Pair Group Method with Arithmetic Mean)
Distance metric: Hamming distance
Visualization: Circular layout with color-coded mutation groups

Mutation groups (rpoB gene):
- S450L: Serine to Leucine at position 450 (blue)
- S531L: Serine to Leucine at position 531 (red)
- H526Y: Histidine to Tyrosine at position 526 (green)
- D516V: Aspartic acid to Valine at position 516 (purple)
"""

import numpy as np
from ete3 import Tree, TreeStyle, TextFace, CircleFace, NodeStyle
import sys


# Mycobacterium tuberculosis drug resistance mutation dataset
# Binary mutation matrix: 15 samples × 20 genomic positions
# 1 = mutation present, 0 = wild-type
MUTATION_MATRIX = np.array([
    [1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1],  # Sample 1
    [1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0],  # Sample 2
    [1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1],  # Sample 3
    [0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0],  # Sample 4
    [0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0],  # Sample 5
    [0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1],  # Sample 6
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0],  # Sample 7
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1],  # Sample 8
    [0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],  # Sample 9
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1],  # Sample 10
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0],  # Sample 11
    [0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1],  # Sample 12
    [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],  # Sample 13
    [1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1],  # Sample 14
    [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],  # Sample 15
])

# Sample labels
SAMPLE_LABELS = [f"Sample {i+1}" for i in range(15)]

# Mutation group assignments for each sample (rpoB gene mutations)
# Based on the mutation patterns in the matrix
MUTATION_GROUPS = {
    "Sample 1": "S450L",   # Group 1: blue
    "Sample 2": "S450L",
    "Sample 3": "S450L",
    "Sample 4": "S531L",   # Group 2: red
    "Sample 5": "S531L",
    "Sample 6": "S531L",
    "Sample 7": "H526Y",   # Group 3: green
    "Sample 8": "H526Y",
    "Sample 9": "H526Y",
    "Sample 10": "D516V",  # Group 4: purple
    "Sample 11": "D516V",
    "Sample 12": "D516V",
    "Sample 13": "S450L",  # Mixed/recombinant
    "Sample 14": "S450L",
    "Sample 15": "S450L",
}

# Color scheme for mutation groups
GROUP_COLORS = {
    "S450L": "#0066CC",  # Blue
    "S531L": "#CC0000",  # Red
    "H526Y": "#009933",  # Green
    "D516V": "#9933CC",  # Purple
}


def compute_hamming_distance(seq1, seq2):
    """
    Compute Hamming distance between two binary sequences.
    
    Hamming distance is the number of positions at which the corresponding
    symbols are different.
    
    Args:
        seq1: First binary sequence (numpy array)
        seq2: Second binary sequence (numpy array)
    
    Returns:
        int: Hamming distance
    """
    return np.sum(seq1 != seq2)


def build_distance_matrix(mutation_matrix):
    """
    Build a pairwise distance matrix using Hamming distance.
    
    Args:
        mutation_matrix: Binary mutation matrix (samples × positions)
    
    Returns:
        numpy.ndarray: Symmetric distance matrix
    """
    n_samples = mutation_matrix.shape[0]
    distance_matrix = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        for j in range(i+1, n_samples):
            dist = compute_hamming_distance(mutation_matrix[i], mutation_matrix[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    
    return distance_matrix


def upgma(distance_matrix, labels):
    """
    Build a phylogenetic tree using UPGMA algorithm.
    
    UPGMA (Unweighted Pair Group Method with Arithmetic Mean) is a simple
    agglomerative hierarchical clustering method.
    
    Args:
        distance_matrix: Symmetric distance matrix
        labels: List of sample labels
    
    Returns:
        ete3.Tree: Phylogenetic tree in Newick format
    """
    n = len(labels)
    # Create a copy of the distance matrix
    dist = distance_matrix.copy()
    
    # Initialize clusters - each sample is its own cluster
    clusters = {i: Tree(name=labels[i]) for i in range(n)}
    cluster_sizes = {i: 1 for i in range(n)}
    
    # Keep track of active clusters
    active = set(range(n))
    next_id = n
    
    while len(active) > 1:
        # Find the pair of clusters with minimum distance
        min_dist = float('inf')
        min_i, min_j = -1, -1
        
        for i in active:
            for j in active:
                if i < j and dist[i, j] < min_dist:
                    min_dist = dist[i, j]
                    min_i, min_j = i, j
        
        # Create new internal node
        new_tree = Tree()
        new_tree.add_child(clusters[min_i], dist=min_dist / 2)
        new_tree.add_child(clusters[min_j], dist=min_dist / 2)
        
        # Update clusters
        clusters[next_id] = new_tree
        new_size = cluster_sizes[min_i] + cluster_sizes[min_j]
        cluster_sizes[next_id] = new_size
        
        # Update distance matrix using UPGMA formula
        # Distance from new cluster to all other clusters
        new_distances = {}
        for k in active:
            if k != min_i and k != min_j:
                # Average distance weighted by cluster sizes
                d_ik = dist[min(min_i, k), max(min_i, k)]
                d_jk = dist[min(min_j, k), max(min_j, k)]
                new_dist = (cluster_sizes[min_i] * d_ik + cluster_sizes[min_j] * d_jk) / new_size
                new_distances[k] = new_dist
        
        # Expand distance matrix if needed
        if next_id >= dist.shape[0]:
            new_dist_matrix = np.zeros((next_id + 1, next_id + 1))
            new_dist_matrix[:dist.shape[0], :dist.shape[1]] = dist
            dist = new_dist_matrix
        
        # Update distances in matrix
        for k, d in new_distances.items():
            dist[min(next_id, k), max(next_id, k)] = d
            dist[max(next_id, k), min(next_id, k)] = d
        
        # Remove merged clusters and add new cluster
        active.remove(min_i)
        active.remove(min_j)
        active.add(next_id)
        next_id += 1
    
    # Return the root of the tree
    root_id = list(active)[0]
    return clusters[root_id]


def style_tree(tree, mutation_groups, group_colors):
    """
    Apply publication-quality styling to the phylogenetic tree.
    
    Args:
        tree: ete3.Tree object
        mutation_groups: Dictionary mapping sample names to mutation groups
        group_colors: Dictionary mapping mutation groups to colors
    """
    # Configure tree style
    ts = TreeStyle()
    ts.mode = "c"  # Circular mode
    ts.show_leaf_name = False  # We'll add custom labels
    ts.show_branch_length = True
    ts.show_branch_support = False
    ts.branch_vertical_margin = 10  # Spacing between branches
    ts.arc_start = 0  # Start angle for circular layout
    ts.arc_span = 360  # Full circle
    
    # Title
    title = TextFace("Mycobacterium tuberculosis Drug Resistance Mutations", fsize=16, bold=True)
    ts.title.add_face(title, column=0)
    
    # Legend
    legend_text = (
        "rpoB Gene Mutation Groups:\n"
        "S450L (Blue) | S531L (Red) | H526Y (Green) | D516V (Purple)"
    )
    legend = TextFace(legend_text, fsize=10)
    ts.legend.add_face(legend, column=0)
    
    # Style each node
    for node in tree.traverse():
        # Set branch style
        nstyle = NodeStyle()
        nstyle["hz_line_width"] = 2
        nstyle["vt_line_width"] = 2
        nstyle["size"] = 0  # Hide node circles
        node.set_style(nstyle)
        
        # Style leaf nodes
        if node.is_leaf():
            sample_name = node.name
            mutation_group = mutation_groups.get(sample_name, "Unknown")
            color = group_colors.get(mutation_group, "#000000")
            
            # Add colored circle
            circle = CircleFace(radius=6, color=color, style="sphere")
            node.add_face(circle, column=0, position="branch-right")
            
            # Add sample label with color
            label = TextFace(f" {sample_name}", fsize=11, bold=True, fgcolor=color)
            node.add_face(label, column=1, position="branch-right")
            
            # Add mutation group label
            group_label = TextFace(f"  ({mutation_group})", fsize=9, fgcolor=color)
            node.add_face(group_label, column=2, position="branch-right")
    
    return ts


def main():
    """
    Main function to generate the phylogenetic tree.
    """
    print("=" * 70)
    print("Publication-ready Circular Phylogenetic Tree Generator")
    print("=" * 70)
    print()
    
    # Step 1: Compute distance matrix
    print("Step 1: Computing Hamming distance matrix...")
    distance_matrix = build_distance_matrix(MUTATION_MATRIX)
    print(f"  ✓ Distance matrix computed ({distance_matrix.shape[0]} × {distance_matrix.shape[1]})")
    print()
    
    # Step 2: Build UPGMA tree
    print("Step 2: Building UPGMA phylogenetic tree...")
    tree = upgma(distance_matrix, SAMPLE_LABELS)
    print(f"  ✓ Tree constructed with {len(SAMPLE_LABELS)} taxa")
    print()
    
    # Step 3: Apply styling
    print("Step 3: Applying publication-quality styling...")
    tree_style = style_tree(tree, MUTATION_GROUPS, GROUP_COLORS)
    print("  ✓ Circular layout configured")
    print("  ✓ Color-coded labels applied")
    print("  ✓ Branch lengths displayed")
    print()
    
    # Step 4: Render and save
    print("Step 4: Rendering high-resolution image...")
    output_file = "circular_phylogenetic_tree.png"
    
    try:
        # Render with high DPI for publication quality
        tree.render(output_file, w=2000, h=2000, units="px", tree_style=tree_style, dpi=300)
        print(f"  ✓ Tree saved to '{output_file}'")
        print(f"  ✓ Resolution: 2000×2000 pixels at 300 DPI")
    except Exception as e:
        print(f"  ✗ Error rendering tree: {e}")
        print()
        print("Note: If rendering fails, ensure you have the required dependencies:")
        print("  - PyQt5 or PyQt6 (for GUI rendering)")
        print("  - Install with: pip install ete3[all] or pip install PyQt5")
        sys.exit(1)
    
    print()
    print("=" * 70)
    print("✓ Phylogenetic tree generation completed successfully!")
    print("=" * 70)
    print()
    print("Summary:")
    print(f"  • Samples analyzed: {len(SAMPLE_LABELS)}")
    print(f"  • Mutation positions: {MUTATION_MATRIX.shape[1]}")
    print(f"  • Mutation groups: {len(set(MUTATION_GROUPS.values()))}")
    print(f"  • Output file: {output_file}")
    print()
    
    # Print mutation group distribution
    print("Mutation group distribution:")
    for group in sorted(set(MUTATION_GROUPS.values())):
        count = sum(1 for v in MUTATION_GROUPS.values() if v == group)
        color = GROUP_COLORS[group]
        print(f"  • {group}: {count} samples (color: {color})")
    print()


if __name__ == "__main__":
    main()
