import random

# Parameters
num_clusters = 2
nodes_per_cluster = 10  # Nodes per cluster
inter_cluster_edges = 40  # Random edges between clusters
space_size = 100  # Space range for clusters
edge_weight_range = (1, 3)
spaceness = 5

# Data storage
nodes = []
edges = []
node_index = 0  # Global node index tracker
cluster_centers = []

# Generate cluster centers
for _ in range(num_clusters):
    cluster_centers.append(
        (
            random.uniform(-space_size, space_size),
            random.uniform(-space_size, space_size),
            random.uniform(-space_size, space_size),
        )
    )

# Generate nodes
cluster_nodes = []  # Store node indices for each cluster
for center in cluster_centers:
    cluster = []
    for _ in range(nodes_per_cluster):
        x = center[0] + random.uniform(-spaceness, spaceness)
        y = center[1] + random.uniform(-spaceness, spaceness)
        z = center[2] + random.uniform(-spaceness, spaceness)
        nodes.append((x, y, z))
        cluster.append(node_index)
        node_index += 1
    cluster_nodes.append(cluster)

# Generate intra-cluster edges
for cluster in cluster_nodes:
    for i in range(len(cluster)):
        for j in range(i + 1, len(cluster)):  # Connect nodes within the cluster
            if random.random() < 0.4:  # 40% chance of an edge
                weight = round(random.uniform(*edge_weight_range), 2)
                edges.append((cluster[i], cluster[j], weight))

# Generate inter-cluster edges
for _ in range(inter_cluster_edges * num_clusters):
    c1, c2 = random.sample(range(num_clusters), 2)  # Pick two clusters
    n1 = random.choice(cluster_nodes[c1])
    n2 = random.choice(cluster_nodes[c2])
    weight = round(random.uniform(*edge_weight_range), 2)
    edges.append((n1, n2, weight))

# Write nodes to file
with open(
    f"nodes_nrclts={num_clusters}_npcltrs{nodes_per_cluster}_edgesprct={inter_cluster_edges}.txt",
    "w",
) as f:
    for node in nodes:
        f.write(f"{node[0]} {node[1]} {node[2]}\n")

# Write edges to file
with open(
    f"edges_nrclts={num_clusters}_npcltrs{nodes_per_cluster}_edgesprct={inter_cluster_edges}.txt",
    "w",
) as f:
    for edge in edges:
        f.write(f"{edge[0]} {edge[1]} {edge[2]}\n")

print("Graph generation complete! Files saved as nodes.txt and edges.txt.")
