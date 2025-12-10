#!/usr/bin/env python3
"""
Circos Plot for KEGG Pathways Visualization
============================================

This script generates a Circos plot visualizing closely connected KEGG pathways
based on the Mycobacterium tuberculosis drug resistance mutations dataset.

Requirements:
    - matplotlib: Install via `pip install matplotlib`
    - numpy: Install via `pip install numpy`

Note: While pycircos is available, this script uses matplotlib directly for better
      compatibility and control over the visualization.

Usage:
    python circos_kegg_pathways.py

Output:
    circos_kegg_pathways.png
"""

import sys

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import numpy as np
except ImportError as e:
    print(f"Error: {e}")
    print("\nPlease install the required packages:")
    print("  pip install matplotlib numpy")
    sys.exit(1)


# Configuration constants
GAP_DEGREES = 5  # Gap between sectors in degrees
TEXT_ROTATION_ADJUSTMENT = 90  # Text rotation adjustment in degrees
OUTPUT_DPI = 300  # Output image DPI


def draw_bezier_curve(ax, start_angle, end_angle, radius, color, alpha=0.5):
    """Draw a Bezier curve connecting two points on the circle."""
    # Calculate start and end points
    start_x = radius * np.cos(start_angle)
    start_y = radius * np.sin(start_angle)
    end_x = radius * np.cos(end_angle)
    end_y = radius * np.sin(end_angle)
    
    # Control points at the center for a nice curve
    ctrl_x1, ctrl_y1 = 0, 0
    ctrl_x2, ctrl_y2 = 0, 0
    
    # Create the bezier curve
    from matplotlib.path import Path
    verts = [
        (start_x, start_y),
        (ctrl_x1, ctrl_y1),
        (ctrl_x2, ctrl_y2),
        (end_x, end_y),
    ]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor=color, 
                               alpha=alpha, linewidth=2)
    ax.add_patch(patch)


def main():
    """Main function to generate the Circos plot."""
    
    # Define KEGG pathways associated with mutated genes in M. tuberculosis
    pathways = [
        ('Transcription', 'rpoB', '#E74C3C'),
        ('DNA Repair', 'katG', '#3498DB'),
        ('Fatty Acid Biosynthesis', 'inhA', '#2ECC71'),
        ('Nicotinate Metabolism', 'pncA', '#F39C12'),
        ('Glycan Biosynthesis', 'embB', '#9B59B6'),
        ('Oxidative Stress', 'ahpC', '#1ABC9C')
    ]
    
    # Define connections between closely related pathways
    # Format: (source_pathway_idx, target_pathway_idx, color)
    connections = [
        (0, 1, '#E8B4B4'),  # Transcription to DNA Repair - Light red
        (2, 4, '#A8D5BA'),  # Fatty Acid Biosynthesis to Glycan Biosynthesis - Light green
        (3, 2, '#F5D6A8'),  # Nicotinate Metabolism to Fatty Acid Biosynthesis - Light orange
        (4, 5, '#C5B4D5')   # Glycan Biosynthesis to Oxidative Stress - Light purple
    ]
    
    print("Preparing pathway data...")
    n_pathways = len(pathways)
    
    # Calculate angles for each pathway sector
    gap_angle = np.radians(GAP_DEGREES)
    sector_angle = (2 * np.pi - n_pathways * gap_angle) / n_pathways
    
    # Create figure and axis
    print("Initializing Circos plot...")
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    
    # Outer and inner radius for sectors
    outer_radius = 1.0
    inner_radius = 0.85
    label_radius = 1.15
    
    # Draw each pathway sector
    print("Adding pathway sectors...")
    sector_angles = []
    
    for idx, (pathway_name, gene, color) in enumerate(pathways):
        # Calculate start and end angles for this sector
        start_angle = idx * (sector_angle + gap_angle)
        end_angle = start_angle + sector_angle
        sector_angles.append((start_angle, end_angle))
        
        # Create sector (wedge)
        theta = np.linspace(start_angle, end_angle, 100)
        
        # Outer arc
        outer_x = outer_radius * np.cos(theta)
        outer_y = outer_radius * np.sin(theta)
        
        # Inner arc
        inner_x = inner_radius * np.cos(theta)
        inner_y = inner_radius * np.sin(theta)
        
        # Create polygon for the sector
        verts = list(zip(outer_x, outer_y)) + list(zip(inner_x[::-1], inner_y[::-1]))
        sector_patch = patches.Polygon(verts, facecolor=color, edgecolor='#333333',
                                       linewidth=2, alpha=0.85)
        ax.add_patch(sector_patch)
        
        # Add pathway label
        mid_angle = (start_angle + end_angle) / 2
        label_x = label_radius * np.cos(mid_angle)
        label_y = label_radius * np.sin(mid_angle)
        
        # Calculate rotation for text
        rotation = np.degrees(mid_angle)
        if 90 < rotation < 270:
            rotation = rotation - 180
            ha = 'right'
        else:
            ha = 'left'
        
        ax.text(label_x, label_y, f'{pathway_name}\n({gene})',
                rotation=rotation - TEXT_ROTATION_ADJUSTMENT,
                ha=ha, va='center',
                fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                         edgecolor='none', alpha=0.7))
    
    # Draw connections between pathways
    print("Adding pathway connections...")
    for source_idx, target_idx, link_color in connections:
        # Get midpoint angles for source and target
        source_mid = (sector_angles[source_idx][0] + sector_angles[source_idx][1]) / 2
        target_mid = (sector_angles[target_idx][0] + sector_angles[target_idx][1]) / 2
        
        # Draw bezier curve from source to target
        draw_bezier_curve(ax, source_mid, target_mid, inner_radius - 0.05, 
                         link_color, alpha=0.6)
    
    # Add title
    plt.title('KEGG Pathways in M. tuberculosis Drug Resistance\n'
              'Circular Visualization of Pathway Connections',
              fontsize=16, fontweight='bold', pad=20)
    
    # Save the plot
    output_file = 'circos_kegg_pathways.png'
    print(f"Saving plot to {output_file}...")
    plt.savefig(output_file, dpi=OUTPUT_DPI, bbox_inches='tight', facecolor='white')
    print(f"✓ Circos plot successfully saved to {output_file}")
    
    # Close the plot
    plt.close()
    
    # Optionally display the plot
    # plt.show()


if __name__ == '__main__':
    main()
