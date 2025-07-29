import os
import sys
import math
import numpy as np
from collections import defaultdict
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


def bresenham_3d(start, end):
    """Generates all voxels along a 3D line from start to end using Bresenham's algorithm."""
    x1, y1, z1 = map(int, start)
    x2, y2, z2 = map(int, end)

    voxels = []

    dx = abs(x2 - x1)
    dy = abs(y2 - y1)
    dz = abs(z2 - z1)

    xs = 1 if x2 > x1 else -1
    ys = 1 if y2 > y1 else -1
    zs = 1 if z2 > z1 else -1

    x, y, z = x1, y1, z1

    # Driving axis is X-axis
    if dx >= dy and dx >= dz:
        p1 = 2 * dy - dx
        p2 = 2 * dz - dx
        for _ in range(dx + 1):
            voxels.append((x, y, z))
            if p1 >= 0:
                y += ys
                p1 -= 2 * dx
            if p2 >= 0:
                z += zs
                p2 -= 2 * dx
            p1 += 2 * dy
            p2 += 2 * dz
            x += xs

    # Driving axis is Y-axis
    elif dy >= dx and dy >= dz:
        p1 = 2 * dx - dy
        p2 = 2 * dz - dy
        for _ in range(dy + 1):
            voxels.append((x, y, z))
            if p1 >= 0:
                x += xs
                p1 -= 2 * dy
            if p2 >= 0:
                z += zs
                p2 -= 2 * dy
            p1 += 2 * dx
            p2 += 2 * dz
            y += ys

    # Driving axis is Z-axis
    else:
        p1 = 2 * dy - dz
        p2 = 2 * dx - dz
        for _ in range(dz + 1):
            voxels.append((x, y, z))
            if p1 >= 0:
                y += ys
                p1 -= 2 * dz
            if p2 >= 0:
                x += xs
                p2 -= 2 * dz
            p1 += 2 * dy
            p2 += 2 * dx
            z += zs

    return voxels


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
    """Maps edges to voxels, including all intermediate voxels along the edge path."""
    voxel_edges = {}

    for j, edge in enumerate(edges):
        visited_voxels = set()
        for i in range(1, len(edge)):
            start_voxel = tuple(voxel_indices[edge[i - 1]])
            end_voxel = tuple(voxel_indices[edge[i]])
            traversed_voxels = bresenham_3d(start_voxel, end_voxel)

            for voxel in traversed_voxels:
                if (voxel, j) not in visited_voxels:
                    voxel_edges.setdefault(voxel, set()).add(j)
                    visited_voxels.add((voxel, j))

    return voxel_edges


def compute_overplotted_percentage(voxel_edges):
    """Computes the percentage of overplotted voxels (voxels with more than one edge)."""
    used_voxels = set(voxel_edges.keys())
    overplotted_voxels = [v for v, es in voxel_edges.items() if len(es) > 1]

    return 100 * len(overplotted_voxels) / len(used_voxels) if used_voxels else 0


def quantize_point(point, precision=3):
    """Rounds the point to the given decimal precision (default: 3 = 0.001 resolution)."""
    return tuple(np.round(point, decimals=precision))

def compute_overcrowded_percentage(edges, total_edges, nodes: list[list[float, float, float]], precision=2):
    """
    Computes the percentage of edges that are bundled.
    Bundling means sharing at least one intermediate node (rounded to 0.001 by default) with another edge.
    """
    point_to_edges = defaultdict(set)

    # Step 1: Build map of quantized interior points → edge indices
    for edge_idx, edge in enumerate(edges):
        if len(edge) <= 2:
            continue

        for node_idx in edge[1:-1]:  # skip first and last
            point = quantize_point(nodes[node_idx], precision)
            point_to_edges[point].add(edge_idx)

    # Step 2: Identify bundled edges
    bundled_edges = set()

    for edge_idx, edge in enumerate(edges):
        if len(edge) <= 2:
            continue

        for node_idx in edge[1:-1]:
            point = quantize_point(nodes[node_idx], precision)
            if len(point_to_edges[point]) > 1:
                bundled_edges.add(edge_idx)
                break

    return 100 * len(bundled_edges) / total_edges if total_edges else 0


def compute_ink_paper_ratio(voxel_edges, num_voxels):
    """Computes the Ink-Paper Ratio: used voxels divided by available voxels."""
    used_voxels = len(voxel_edges)
    available_voxels = num_voxels**3
    return used_voxels / available_voxels


def compute_edge_smoothness(nodes, edges):
    total_avg_bend = 0.0
    total_length_ratio = 0.0
    total_composite_score = 0.0
    count = 0

    for edge in edges:
        if len(edge) < 3:
            continue  # Skip straight lines

        # Edge points
        pts = nodes[edge]

        # Average bend
        bend_sum = 0.0
        bend_count = 0
        for i in range(1, len(pts) - 1):
            v1 = pts[i] - pts[i - 1]
            v2 = pts[i + 1] - pts[i]
            norm1 = np.linalg.norm(v1)
            norm2 = np.linalg.norm(v2)

            if norm1 < 1e-6 or norm2 < 1e-6:
                continue

            v1 /= norm1
            v2 /= norm2
            dot = np.clip(np.dot(v1, v2), -1.0, 1.0)
            angle = math.acos(dot)
            bend_sum += angle
            bend_count += 1

        avg_bend = bend_sum / bend_count if bend_count > 0 else 0.0

        # Edge lengths
        bundled_len = np.sum(np.linalg.norm(np.diff(pts, axis=0), axis=1))
        straight_len = np.linalg.norm(pts[-1] - pts[0])
        length_ratio = bundled_len / straight_len if straight_len > 0 else 1.0

        composite = avg_bend * length_ratio

        total_avg_bend += avg_bend
        total_length_ratio += length_ratio
        total_composite_score += composite
        count += 1

    if count == 0:
        return {
            "Avg Bend (rad)": 0.0,
            "Avg Length Ratio": 1.0,
            "Smoothness Score": 0.0,
        }

    return {
        "Avg Bend (rad)": total_avg_bend / count,
        "Avg Length Ratio": (total_length_ratio / count).item(),
        "Smoothness Score": (total_composite_score / count).item(),
    }


def compute_metrics(
    nodes, edges, min_coords, max_coords, voxel_size, voxel_indices, num_voxels=10
):
    """Main function to compute all three metrics."""

    # Step 1: Map edges to voxels
    voxel_edges = map_edges_to_voxels(nodes, edges, voxel_indices)

    # Step 2: Compute metrics
    overplotted_percent = compute_overplotted_percentage(voxel_edges)
    overcrowded_percent = compute_overcrowded_percentage(edges, len(edges), nodes)
    ink_paper_ratio = compute_ink_paper_ratio(voxel_edges, num_voxels)
    smoothness = compute_edge_smoothness(nodes, edges)

    return {
        "Overplotted%": overplotted_percent,
        "Overcrowded%": overcrowded_percent,
        "Ink-Paper Ratio": ink_paper_ratio,
        **smoothness,  # Merges the smoothness metrics
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


def plot_voxels(ax, voxel_indices, voxel_size, min_coords, color_map=None, alpha=0.3):
    """Plots voxels from a list of indices."""
    unique_voxels = np.unique(voxel_indices, axis=0)

    if color_map is None:
        colors = list(mcolors.TABLEAU_COLORS.values())
        color_map = {
            tuple(voxel): colors[i % len(colors)]
            for i, voxel in enumerate(unique_voxels)
        }

    for voxel in unique_voxels:
        voxel_corner = min_coords + voxel * voxel_size
        x, y, z = voxel_corner
        size = voxel_size[0]  # Assuming cubic voxels

        faces = get_voxel_faces(x, y, z, size)
        voxel_poly = Poly3DCollection(
            faces, alpha=alpha, linewidths=0.3, edgecolors="black"
        )
        voxel_poly.set_facecolor(color_map[tuple(voxel)])
        ax.add_collection3d(voxel_poly)


def plot_nodes(ax, nodes, color="blue", label="Nodes"):
    """Plots nodes as 3D scatter."""
    ax.scatter(nodes[:, 0], nodes[:, 1], nodes[:, 2], c=color, marker="o", label=label)


def plot_edges(ax, nodes, edges, color="black", alpha=0.5):
    """Plots edges between nodes."""
    for edge in edges:
        edge_nodes = nodes[edge]
        ax.plot(
            edge_nodes[:, 0], edge_nodes[:, 1], edge_nodes[:, 2], c=color, alpha=alpha
        )


def visualize_voxel_grid(
    nodes, edges, voxel_indices, voxel_size, min_coords, show_edges=True
):
    """Main visualization wrapper."""
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection="3d")

    plot_nodes(ax, nodes)
    plot_voxels(ax, voxel_indices, voxel_size, min_coords)

    if show_edges:
        plot_edges(ax, nodes, edges)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # Check if a filename was provided
    if len(sys.argv) < 2:
        print("Usage: python script.py <filename>")
        sys.exit(1)

    # Load VTK file
    VTK_FILE = sys.argv[1]
    _, nodes, edges, _ = read_vtk_ascii(VTK_FILE)

    NUM_VOXELS = 50

    # Get voxel grid parameters
    min_coords, max_coords = get_bounding_box(nodes)

    voxel_size = get_voxel_size(min_coords, max_coords, num_voxels=NUM_VOXELS)

    # Assign nodes to voxels
    voxel_indices = map_nodes_to_voxels(nodes, min_coords, voxel_size)

    # Visualize in 3D
    # visualize_voxel_grid(nodes, edges, voxel_indices, voxel_size, min_coords)
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
