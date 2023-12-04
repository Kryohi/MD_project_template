import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# Function to parse the PDB file and extract alpha carbon atoms
def parse_pdb_for_ca_atoms(pdb_file):
    ca_atoms = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and " CA " in line:
                residue_index = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                ca_atoms.append((residue_index, np.array([x, y, z])))
    return ca_atoms

# Function to calculate Euclidean distance
def euclidean_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

# Function to calculate the path centrality for a pair of nodes
def path_centrality(G, node1, node2, centrality_scores):
    try:
        shortest_path = nx.shortest_path(G, source=node1, target=node2)
        return sum(centrality_scores[node] for node in shortest_path)
    except nx.NetworkXNoPath:
        return None

# Parse the PDB file
pdb_file_path = '/path/to/your/pdb/file.pdb'  # Change this to your PDB file path
ca_atoms = parse_pdb_for_ca_atoms(pdb_file_path)

# Create the graph
G = nx.Graph()
for i in range(len(ca_atoms)):
    for j in range(i + 1, len(ca_atoms)):
        dist = euclidean_distance(ca_atoms[i][1], ca_atoms[j][1])
        if dist < 5.0:  # 5 Angstroms threshold
            G.add_edge(ca_atoms[i][0], ca_atoms[j][0])

# Visualization of the graph
plt.figure(figsize=(12, 12))
pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, node_size=100, node_color='blue')
nx.draw_networkx_edges(G, pos, alpha=0.5)
nx.draw_networkx_labels(G, pos, font_size=8, font_color="black")
plt.title("Protein Structure Network Graph")
plt.axis('off')
plt.show()

# Centrality measures for mutated residues
mutations = {"Ile143Val": 143, "Lys178Arg": 178, "Leu84Phe": 84}
residues = list(mutations.values())

# Compute betweenness centrality
betweenness_centrality = nx.betweenness_centrality(G)

# Compute other centrality measures
degree_centrality = nx.degree_centrality(G)
closeness_centrality = nx.closeness_centrality(G)
eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
katz_centrality = nx.katz_centrality(G, max_iter=1000)

# Collect centrality measures for mutated residues
centrality_measures = {}
for residue in residues:
    if residue in G.nodes:
        centrality_measures[residue] = {
            "Betweenness Centrality": betweenness_centrality[residue],
            "Degree Centrality": degree_centrality[residue],
            "Closeness Centrality": closeness_centrality[residue],
            "Eigenvector Centrality": eigenvector_centrality[residue],
            "Katz Centrality": katz_centrality[residue]
        }

# Calculate path centrality for all pairs of mutated residues
path_centrality_results = {}
for i in range(len(residues)):
    for j in range(i + 1, len(residues)):
        pair = (residues[i], residues[j])
        path_cent = path_centrality(G, residues[i], residues[j], betweenness_centrality)
        path_centrality_results[pair] = path_cent

print(centrality_measures)
print(path_centrality_results)

