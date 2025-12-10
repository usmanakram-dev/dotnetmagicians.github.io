#!/usr/bin/env python3
"""
Phylogenetic Tree Builder for Mycobacterium tuberculosis Drug Resistance Mutations

This script processes drug resistance mutation data to build a circular phylogenetic tree
using UPGMA clustering based on Hamming distances between samples.
"""

import os
# Set Qt platform to offscreen for headless environments
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, to_tree
from ete3 import Tree
from ete3.treeview import TreeStyle, NodeStyle, TextFace


# Dataset definition
MUTATIONS = [
    "rpoB_S450L", "katG_S315T", "inhA_-15C_T", "rpoB_S531L", "inhA_-8T_C",
    "rpoB_H526Y", "rpoB_D516V", "pncA_H57D", "pncA_D49A", "pncA_G132S",
    "katG_R463L", "pncA_A102V", "ahpC_-46A_G", "pncA_L182P", "embB_M306V", "embB_M306I"
]

SAMPLES = {
    "Sample_1": [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_2": [1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_3": [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_4": [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_5": [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_6": [0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_7": [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_8": [0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    "Sample_9": [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    "Sample_10": [0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    "Sample_11": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],
    "Sample_12": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
    "Sample_13": [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    "Sample_14": [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "Sample_15": [0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
}


def compute_hamming_distance_matrix(samples_dict):
    """
    Compute Hamming distance matrix between all samples.
    
    Args:
        samples_dict: Dictionary mapping sample names to binary vectors
        
    Returns:
        tuple: (distance_matrix, sample_names)
    """
    sample_names = list(samples_dict.keys())
    vectors = np.array([samples_dict[name] for name in sample_names])
    
    # Compute pairwise Hamming distances
    distances = pdist(vectors, metric='hamming')
    distance_matrix = squareform(distances)
    
    return distance_matrix, sample_names


def build_upgma_tree(distance_matrix, sample_names):
    """
    Build a UPGMA tree from a distance matrix.
    
    Args:
        distance_matrix: Square distance matrix
        sample_names: List of sample names
        
    Returns:
        ete3.Tree: UPGMA tree
    """
    # Perform UPGMA clustering using condensed distance matrix
    condensed_dist = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_dist, method='average')
    
    # Convert scipy linkage to ete3 Tree
    tree = linkage_to_ete3(linkage_matrix, sample_names)
    
    return tree


def linkage_to_ete3(linkage_matrix, labels):
    """
    Convert scipy linkage matrix to ete3 Tree object.
    
    Args:
        linkage_matrix: Scipy linkage matrix
        labels: List of leaf labels
        
    Returns:
        ete3.Tree: Converted tree
    """
    # Build tree from scipy linkage matrix
    scipy_tree = to_tree(linkage_matrix, rd=False)
    
    # Convert to ete3 tree
    def build_ete3_tree(node, node_id):
        if node.is_leaf():
            return Tree(name=labels[node_id])
        else:
            ete_node = Tree()
            ete_node.add_child(build_ete3_tree(node.left, node.left.id))
            ete_node.add_child(build_ete3_tree(node.right, node.right.id))
            return ete_node
    
    return build_ete3_tree(scipy_tree, scipy_tree.id)


def get_leaf_color(sample_name, samples_dict):
    """
    Determine leaf color based on rpoB mutation group.
    
    Color scheme:
    - Blue: S450L (mutation at index 0)
    - Red: S531L (mutation at index 3)
    - Green: H526Y (mutation at index 5)
    - Purple: D516V (mutation at index 6)
    - Black: None of the above
    
    Args:
        sample_name: Name of the sample
        samples_dict: Dictionary mapping sample names to binary vectors
        
    Returns:
        str: Color name
    """
    vector = samples_dict[sample_name]
    
    # Check for rpoB mutations
    if vector[0] == 1:  # rpoB_S450L
        return "blue"
    elif vector[3] == 1:  # rpoB_S531L
        return "red"
    elif vector[5] == 1:  # rpoB_H526Y
        return "green"
    elif vector[6] == 1:  # rpoB_D516V
        return "purple"
    else:
        return "black"


def render_circular_tree(tree, samples_dict, output_file="phylogenetic_tree_circular.png"):
    """
    Render the tree in circular layout with color-coded leaves.
    
    Args:
        tree: ete3.Tree object
        samples_dict: Dictionary mapping sample names to binary vectors
        output_file: Output PNG file path
    """
    # Set up tree style
    ts = TreeStyle()
    ts.mode = "c"  # Circular mode
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.arc_start = 0
    ts.arc_span = 360
    
    # Style the leaves based on rpoB mutation groups
    for leaf in tree.iter_leaves():
        # Get color for this leaf
        color = get_leaf_color(leaf.name, samples_dict)
        
        # Create node style with colored text
        nstyle = NodeStyle()
        nstyle["size"] = 0  # Hide node circle
        
        # Add colored text face
        text_face = TextFace(leaf.name, fgcolor=color, fsize=12)
        leaf.add_face(text_face, column=0, position="branch-right")
        
        leaf.set_style(nstyle)
    
    # Render tree to file
    tree.render(output_file, w=800, h=800, tree_style=ts)
    print(f"Circular phylogenetic tree saved to {output_file}")


def main():
    """Main execution function."""
    print("Building circular phylogenetic tree for M. tuberculosis drug resistance mutations...")
    print(f"Processing {len(SAMPLES)} samples with {len(MUTATIONS)} mutations")
    
    # Compute distance matrix
    print("\nComputing Hamming distance matrix...")
    distance_matrix, sample_names = compute_hamming_distance_matrix(SAMPLES)
    
    # Build UPGMA tree
    print("Building UPGMA tree...")
    tree = build_upgma_tree(distance_matrix, sample_names)
    
    # Render circular tree
    print("Rendering circular tree with color-coded leaves...")
    render_circular_tree(tree, SAMPLES)
    
    print("\nDone!")


if __name__ == "__main__":
    main()
