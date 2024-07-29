import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

def single_linkage_clustering(distance_matrix):
    n = len(distance_matrix)
    clusters = [{i} for i in range(n)]
    steps = []

    while len(clusters) > 1:
        min_dist = float('inf')
        pair_to_merge = (0, 0)

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                min_pair_dist = min(distance_matrix[list(clusters[i])[:, None], list(clusters[j])])
                if min_pair_dist < min_dist:
                    min_dist = min_pair_dist
                    pair_to_merge = (i, j)

        cluster1, cluster2 = pair_to_merge
        clusters[cluster1].update(clusters[cluster2])
        del clusters[cluster2]
        steps.append((cluster1, cluster2, min_dist))

    return steps

def create_linkage_matrix(distance_matrix, steps):
    linkage_matrix = []
    current_cluster = len(distance_matrix)
    cluster_map = {i: i for i in range(len(distance_matrix))}

    for i, j, dist in steps:
        linkage_matrix.append([cluster_map[i], cluster_map[j], dist, len(distance_matrix)])
        cluster_map[i] = current_cluster
        current_cluster += 1

    return np.array(linkage_matrix)

def plot_dendrogram(linkage_matrix, labels):
    plt.figure(figsize=(10, 7))
    dendro = dendrogram(
        linkage_matrix,
        labels=labels,
        leaf_rotation=90,
        leaf_font_size=12,
        show_contracted=True
    )

    for i, d, c in zip(dendro['icoord'], dendro['dcoord'], dendro['color_list']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'o', c=c)
        plt.annotate(f'{y:.2f}', (x, y), xytext=(0, -5), textcoords='offset points',
                     va='top', ha='center', fontsize=10, color=c)

    plt.title("Dendrograma de enlace único")
    plt.xlabel("Índice de muestra")
    plt.ylabel("Distancia")
    plt.show()
