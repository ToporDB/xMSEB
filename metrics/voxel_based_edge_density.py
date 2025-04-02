import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as mcolors


def read_vtk_ascii(filename):

    if os.path.exists(filename):
        print("File exists")
    else:
        print("File does not exist")

    """Reads a VTK (ASCII) file and extracts points & edges"""
    with open(filename, "r") as file:
        lines = file.readlines()

    header = lines[:4]  # Preserve original header

    points_start = next(i for i, l in enumerate(lines) if l.startswith("POINTS"))
    num_points = int(lines[points_start].split()[1])

    # Read points until we hit another section (like LINES)
    points = []
    i = points_start + 1
    while i < len(lines) and not lines[i].startswith(("LINES", "POLYGONS", "CELLS")):
        points.extend(map(float, lines[i].split()))
        i += 1
    points = np.array(points).reshape(-1, 3)

    # Read edges
    edges_start = next((i for i, l in enumerate(lines) if l.startswith("LINES")), None)
    edges = []
    if edges_start:
        num_edges = int(lines[edges_start].split()[1])
        i = edges_start + 1
        while len(edges) < num_edges:
            parts = list(map(int, lines[i].split()))
            edges.append(parts[1:])  # Skip first value (number of indices in edge)
            i += 1
        edges = np.array(edges, dtype=int)

    return header, points, edges, lines[edges_start:] if edges_start else []


def get_bounding_box(nodes):
    """Finds the bounding box that contains all the nodes."""
    min_coords = np.min(nodes, axis=0)
    max_coords = np.max(nodes, axis=0)
    return min_coords, max_coords


def get_voxel_size(min_coords, max_coords, num_voxels=10):
    """Calculates the voxel size by dividing the bounding box into smaller voxels."""
    # Calculate the size of the bounding box
    box_size = max_coords - min_coords

    # Divide the bounding box size by the desired number of voxels (e.g., 10x10x10)
    voxel_size = box_size / num_voxels
    return voxel_size


def map_nodes_to_voxels(nodes, min_coords, voxel_size):
    """Assigns each node to a voxel index."""
    return np.floor((nodes - min_coords) / voxel_size).astype(int)


def map_edges_to_voxels(nodes, edges, voxel_indices):
    """Maps edges to voxels and returns a dictionary of voxel-edge mappings."""
    voxel_edges = {}

    for edge in edges:
        for i in range(1, len(edge)):  # Iterate over each edge segment
            start_voxel = tuple(voxel_indices[edge[i - 1]])
            end_voxel = tuple(voxel_indices[edge[i]])

            # Track both start and end voxel
            voxel_edges.setdefault(start_voxel, []).append(edge)
            voxel_edges.setdefault(end_voxel, []).append(edge)

    return voxel_edges


def compute_overplotted_percentage(voxel_edges):
    """Computes the percentage of overplotted voxels (voxels with more than one edge)."""
    used_voxels = set(voxel_edges.keys())
    overplotted_voxels = [v for v, es in voxel_edges.items() if len(es) > 1]

    return 100 * len(overplotted_voxels) / len(used_voxels) if used_voxels else 0


def compute_overcrowded_percentage(voxel_edges, total_edges):
    """Computes the percentage of edges that are bundled (multiple edges in the same voxel)."""
    bundled_edges = sum(len(es) - 1 for es in voxel_edges.values() if len(es) > 1)
    return 100 * bundled_edges / total_edges if total_edges else 0


def compute_ink_paper_ratio(voxel_edges, num_voxels):
    """Computes the Ink-Paper Ratio: used voxels divided by available voxels."""
    used_voxels = len(voxel_edges)
    available_voxels = num_voxels**3
    return used_voxels / available_voxels


def compute_metrics(
    nodes, edges, min_coords, max_coords, voxel_size, voxel_indices, num_voxels=10
):
    """Main function to compute all three metrics."""

    # Step 1: Map edges to voxels
    voxel_edges = map_edges_to_voxels(nodes, edges, voxel_indices)

    # Step 2: Compute metrics
    overplotted_percent = compute_overplotted_percentage(voxel_edges)
    overcrowded_percent = compute_overcrowded_percentage(voxel_edges, len(edges))
    ink_paper_ratio = compute_ink_paper_ratio(voxel_edges, num_voxels)

    return {
        "Overplotted%": overplotted_percent,
        "Overcrowded%": overcrowded_percent,
        "Ink-Paper Ratio": ink_paper_ratio,
    }


def get_voxel_faces(x, y, z, size):
    """Returns the 6 faces of a voxel as polygons for 3D plotting."""
    return [
        [
            [x, y, z],
            [x + size, y, z],
            [x + size, y + size, z],
            [x, y + size, z],
        ],  # Bottom face
        [
            [x, y, z + size],
            [x + size, y, z + size],
            [x + size, y + size, z + size],
            [x, y + size, z + size],
        ],  # Top face
        [
            [x, y, z],
            [x, y, z + size],
            [x + size, y, z + size],
            [x + size, y, z],
        ],  # Front face
        [
            [x, y + size, z],
            [x, y + size, z + size],
            [x + size, y + size, z + size],
            [x + size, y + size, z],
        ],  # Back face
        [
            [x, y, z],
            [x, y, z + size],
            [x, y + size, z + size],
            [x, y + size, z],
        ],  # Left face
        [
            [x + size, y, z],
            [x + size, y, z + size],
            [x + size, y + size, z + size],
            [x + size, y + size, z],
        ],  # Right face
    ]


def visualize_voxel_grid(nodes, edges, voxel_indices, voxel_size, min_coords):
    """Plots the voxel grid and nodes in 3D with correct voxel placement."""
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Plot nodes
    ax.scatter(
        nodes[:, 0], nodes[:, 1], nodes[:, 2], c="blue", marker="o", label="Nodes"
    )

    # Generate colors for voxels
    unique_voxels = np.unique(voxel_indices, axis=0)
    colors = list(mcolors.TABLEAU_COLORS.values())  # Get distinct colors
    color_map = {
        tuple(voxel): colors[i % len(colors)] for i, voxel in enumerate(unique_voxels)
    }

    # Draw voxel grid with correct placement
    for voxel in unique_voxels:
        voxel_corner = (
            min_coords + voxel * voxel_size
        )  # Calculate exact voxel corner position
        x, y, z = voxel_corner
        size = voxel_size[0]  # Assuming cubic voxels

        faces = get_voxel_faces(x, y, z, size)
        color = color_map[tuple(voxel)]
        voxel_poly = Poly3DCollection(
            faces, alpha=0.3, linewidths=0.3, edgecolors="black"
        )
        voxel_poly.set_facecolor(color)
        ax.add_collection3d(voxel_poly)

    # Plot edges
    for edge in edges:
        edge_nodes = nodes[edge]
        ax.plot(
            edge_nodes[:, 0], edge_nodes[:, 1], edge_nodes[:, 2], c="black", alpha=0.5
        )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # Load VTK file
    VTK_FILE = "test_cleaned.vtk"
    _, nodes, edges, _ = read_vtk_ascii(VTK_FILE)

    NUM_VOXELS = 50

    # Get voxel grid parameters
    min_coords, max_coords = get_bounding_box(nodes)

    voxel_size = get_voxel_size(min_coords, max_coords, num_voxels=NUM_VOXELS)

    # Assign nodes to voxels
    voxel_indices = map_nodes_to_voxels(nodes, min_coords, voxel_size)

    # Visualize in 3D
    visualize_voxel_grid(nodes, edges, voxel_indices, voxel_size, min_coords)
    metrics = compute_metrics(
        nodes,
        edges,
        min_coords,
        max_coords,
        voxel_size,
        voxel_indices,
        num_voxels=NUM_VOXELS,
    )

    print(metrics)
