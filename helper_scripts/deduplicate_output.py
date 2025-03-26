import numpy as np


def read_vtk_ascii(filename):
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


def deduplicate_nodes(points, edges):
    """Deduplicate nodes and update edge indices"""
    unique_points, inverse_indices = np.unique(points, axis=0, return_inverse=True)

    # Create a mapping from old indices to new unique indices
    index_mapping = {
        old_idx: new_idx for old_idx, new_idx in enumerate(inverse_indices)
    }

    # Update edges with new indices
    new_edges = np.array([[index_mapping[idx] for idx in edge] for edge in edges])

    return unique_points, new_edges


def write_vtk_ascii(output_filename, header, points, edges):
    """Writes a cleaned VTK file in ASCII format"""
    with open(output_filename, "w") as file:
        file.writelines(header)

        # Write unique points
        file.write(f"POINTS {len(points)} float\n")
        for p in points:
            file.write(f"{p[0]} {p[1]} {p[2]}\n")

        # Write edges
        if len(edges) > 0:
            file.write(
                f"LINES {len(edges)} {len(edges) + sum(len(e) for e in edges)}\n"
            )
            for edge in edges:
                file.write(f"{len(edge)} {' '.join(map(str, edge))}\n")


def process_vtk_ascii(input_filename, output_filename):
    """Main function to process an ASCII VTK file"""
    header, points, edges, _ = read_vtk_ascii(input_filename)
    unique_points, new_edges = deduplicate_nodes(points, edges)
    write_vtk_ascii(output_filename, header, unique_points, new_edges)
    print(f"Processed VTK saved as: {output_filename}")


# Example usage:
input_vtk_file = "test.vtk"
output_vtk_file = "test_cleaned.vtk"
process_vtk_ascii(input_vtk_file, output_vtk_file)
