import sys
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg

def compute_energy_matrix(sequence, alpha):
    length = len(sequence)
    energy_matrix = [[0] * length for _ in range(length)] 

    for gap in range(1, length):
        for i in range(length - gap):
            j = i + gap
            min_energy = float('inf')
            min_energy = min(min_energy, energy_matrix[i + 1][j])
            min_energy = min(min_energy, energy_matrix[i][j - 1])
            if sequence[i] + sequence[j] in alpha:
                min_energy = min(min_energy, energy_matrix[i + 1][j - 1] + alpha[sequence[i] + sequence[j]])
            for k in range(i + 1, j):
                min_energy = min(min_energy, energy_matrix[i][k - 1] + energy_matrix[k][j])

            energy_matrix[i][j] = min_energy

    return energy_matrix

def predict_secondary_structure(sequence, alpha):
    length = len(sequence)
    energy_matrix = compute_energy_matrix(sequence, alpha)
    score = energy_matrix[0][length - 1]

    return energy_matrix, score

def traceback(sequence, energy_matrix, alpha):
    length = len(sequence)
    traceback_pairs = ["."] * length 
    paired_positions = []
    structures = []

    def traceback_helper(i, j):
        nonlocal traceback_pairs, paired_positions, structures

        if i >= j:
            return
        
        if sequence[i] + sequence[j] in alpha and energy_matrix[i][j] == energy_matrix[i + 1][j - 1] + alpha[sequence[i] + sequence[j]]:
            traceback_pairs[i] = "("
            traceback_pairs[j] = ")"
            paired_positions.append((i + 1, j + 1))  
            traceback_helper(i + 1, j - 1)
            return

        if energy_matrix[i][j] == energy_matrix[i + 1][j]:
            traceback_helper(i + 1, j)
            return
    
        if energy_matrix[i][j] == energy_matrix[i][j - 1]:
            traceback_helper(i, j - 1)
            return
        
        for k in range(i + 1, j):
            if energy_matrix[i][j] == energy_matrix[i][k - 1] + energy_matrix[k][j]:
                traceback_helper(i, k - 1)
                traceback_helper(k, j)
                return
    
    traceback_helper(0, length - 1)

    visited = [False] * length
    for i, pair in enumerate(traceback_pairs):
        if pair == "(" and not visited[i]:
            j = traceback_pairs.index(")", i + 1)
            visited[i] = visited[j] = True
            if j - i > 1:
                if all(traceback_pairs[x] == "." for x in range(i + 1, j)):
                    structures.append((i + 1, j + 1, "Lazo: Bucle interno"))
                else:
                    structures.append((i + 1, j + 1, "Tallo"))
            else:
                if i == 0 or j == length - 1:
                    structures.append((i + 1, j + 1, "Lazo: Base no emparejada"))
                else:
                    structures.append((i + 1, j + 1, "Bulbo"))

    return traceback_pairs, paired_positions, structures

def plot_structure(sequence, pairs, structures, filename, node_color="#808080", edge_color="#808080", highlight_color="#FF5733"): 
    plt.figure(figsize=(max(4, len(sequence) / 8), max(4, len(sequence) / 8)))

    G = nx.Graph()
    for idx, base in enumerate(sequence):
        G.add_node(idx + 1, base=base, color=node_color)

    for idx in range(1, len(sequence)):
        G.add_edge(idx, idx + 1, color=edge_color, style="-")

    for pair in pairs:
        G.add_edge(pair[0], pair[1], color=highlight_color, style="--")

    pos = nx.kamada_kawai_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=200, edgecolors=node_color, linewidths=1.5)
    nx.draw_networkx_labels(G, pos, labels=nx.get_node_attributes(G, "base"), font_color="black")

    edges = G.edges()
    edge_colors = [G[u][v]['color'] for u, v in edges]
    edge_styles = [G[u][v]['style'] for u, v in edges]
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=edge_colors, style=edge_styles, width=3)

    # Adding annotations for structures
    for (i, j, struct) in structures:
        x = (pos[i][0] + pos[j][0]) / 2
        y = (pos[i][1] + pos[j][1]) / 2
        plt.text(x, y, struct, fontsize=10, color='red', ha='center')

    plt.axis('off')
    plt.savefig(filename)
    plt.close()

def read_sequence(file_path):
    with open(file_path, 'r') as file:
        sequence = file.readline().strip().upper()
    return sequence

def format_energy_matrix(energy_matrix, sequence):
    formatted_matrix = "     " + " ".join(f"{base:3}" for base in sequence) + "\n"
    for i, row in enumerate(energy_matrix):
        formatted_row = " ".join(f"{cell:3}" for cell in row)
        formatted_matrix += f"{sequence[i]:<3}" + formatted_row + "\n"
    return formatted_matrix

sequence = "TCAAGCGTTAGAGAAGTCATTATGTGATAAAAAAATTCAACTTGGTATCAACTTAA"
alpha_dict = {"CG": -1, "GC": -1, "AU": -1, "UA": -1, "GU": -1, "UG": -1}

energy_matrix, min_energy_score = predict_secondary_structure(sequence, alpha_dict)
traceback_structure, paired_positions, structures = traceback(sequence, energy_matrix, alpha_dict)

plot_structure(sequence, paired_positions, structures, 'structure.png', node_color="#00BFFF", edge_color="#D3D3D3", highlight_color="#FF4500")

with open('output.txt', 'w') as f:
    f.write(f"Puntaje mínimo de energía: {min_energy_score}\n\n")
    f.write("Matriz de energía:\n")
    formatted_matrix = format_energy_matrix(energy_matrix, sequence)
    f.write(formatted_matrix + "\n")
