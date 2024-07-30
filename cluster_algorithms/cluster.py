import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import io
import base64

def min_linkage_clustering(distance_matrix):
    # Asegurarse de que la matriz de distancias es de tipo float
    dist_matrix = distance_matrix.astype(float)

    # Número de elementos iniciales
    n = len(dist_matrix)

    # Mantener un seguimiento de las uniones
    clusters = {i: [i] for i in range(n)}
    steps = []
    level = 0

    
    while len(clusters) > 1:
        # Encuentra los clústeres más cercanos
        min_dist = np.inf
        min_pair = None
        for i in range(len(dist_matrix)):
            for j in range(i + 1, len(dist_matrix)):
                if dist_matrix[i][j] < min_dist:
                    min_dist = dist_matrix[i][j]
                    min_pair = (i, j)

        # Fusión de los clústeres más cercanos
        i, j = min_pair
        steps.append((i, j, min_dist))

        # Unir los clústeres
        clusters[i].extend(clusters[j])
        del clusters[j]

        # Actualizar la matriz de distancias
        for k in range(len(dist_matrix)):
            if k != i and k in clusters:
                dist_matrix[i][k] = dist_matrix[k][i] = min(dist_matrix[i][k], dist_matrix[j][k])

        dist_matrix[:, j] = np.inf
        dist_matrix[j, :] = np.inf

        # Incrementar el nivel
        level += 1

        # Imprimir estado actual
        current_clusters = []
        for cluster in clusters.values():
            merged_cluster = ''.join([chr(65 + x) for x in sorted(cluster)])
            current_clusters.append(merged_cluster)
        
        # Crear una versión achicada de la matriz de distancias
        active_indices = [k for k in range(len(dist_matrix)) if k in clusters]
        reduced_matrix = dist_matrix[np.ix_(active_indices, active_indices)]
        
    return steps

def max_linkage_clustering(distance_matrix):
    # Asegurarse de que la matriz de distancias es de tipo float
    dist_matrix = distance_matrix.astype(float)

    # Número de elementos iniciales
    n = len(dist_matrix)

    # Mantener un seguimiento de las uniones
    clusters = {i: [i] for i in range(n)}
    steps = []
    level = 0

    while len(clusters) > 1:
        # Encuentra los clústeres más cercanos
        min_dist = np.inf
        min_pair = None
        for i in range(len(dist_matrix)):
            for j in range(i + 1, len(dist_matrix)):
                if dist_matrix[i][j] < min_dist:
                    min_dist = dist_matrix[i][j]
                    min_pair = (i, j)

        # Fusión de los clústeres más cercanos
        i, j = min_pair
        steps.append((i, j, min_dist))

        # Unir los clústeres
        clusters[i].extend(clusters[j])
        del clusters[j]

        # Actualizar la matriz de distancias
        for k in range(len(dist_matrix)):
            if k != i and k in clusters:
                dist_matrix[i][k] = dist_matrix[k][i] = max(dist_matrix[i][k], dist_matrix[j][k])

        dist_matrix[:, j] = np.inf
        dist_matrix[j, :] = np.inf

        # Incrementar el nivel
        level += 1

        # Imprimir estado actual
        current_clusters = []
        for cluster in clusters.values():
            merged_cluster = ''.join([chr(65 + x) for x in sorted(cluster)])
            current_clusters.append(merged_cluster)
  
        # Crear una versión achicada de la matriz de distancias
        active_indices = [k for k in range(len(dist_matrix)) if k in clusters]
        reduced_matrix = dist_matrix[np.ix_(active_indices, active_indices)]
 

    return steps

def promedio_linkage_clustering(distance_matrix):
    # Asegurarse de que la matriz de distancias es de tipo float
    dist_matrix = distance_matrix.astype(float)

    # Número de elementos iniciales
    n = len(dist_matrix)

    # Mantener un seguimiento de las uniones
    clusters = {i: [i] for i in range(n)}
    steps = []
    level = 0

   
    while len(clusters) > 1:
        # Encuentra los clústeres más cercanos
        min_dist = np.inf
        min_pair = None
        for i in range(len(dist_matrix)):
            for j in range(i + 1, len(dist_matrix)):
                if dist_matrix[i][j] < min_dist:
                    min_dist = dist_matrix[i][j]
                    min_pair = (i, j)

        # Fusión de los clústeres más cercanos
        i, j = min_pair
        steps.append((i, j, min_dist))

        # Unir los clústeres
        clusters[i].extend(clusters[j])
        del clusters[j]

        # Actualizar la matriz de distancias
        for k in range(len(dist_matrix)):
            if k != i and k in clusters:
                dist_matrix[i][k] = dist_matrix[k][i] = (dist_matrix[i][k]+dist_matrix[j][k])/2

        dist_matrix[:, j] = np.inf
        dist_matrix[j, :] = np.inf

        # Incrementar el nivel
        level += 1

        # Imprimir estado actual
        current_clusters = []
        for cluster in clusters.values():
            merged_cluster = ''.join([chr(65 + x) for x in sorted(cluster)])
            current_clusters.append(merged_cluster)
       
        # Crear una versión achicada de la matriz de distancias
        active_indices = [k for k in range(len(dist_matrix)) if k in clusters]
        reduced_matrix = dist_matrix[np.ix_(active_indices, active_indices)]

    return steps

def plot_dendrogram_to_base64(distance_matrix, steps, title):
    # Convertir los pasos a formato de linkage para dendrograma
    linkage_matrix = []
    current_cluster = len(distance_matrix)
    cluster_map = {i: i for i in range(len(distance_matrix))}

    for i, j, dist in steps:
        linkage_matrix.append([cluster_map[i], cluster_map[j], dist, len(distance_matrix)])
        cluster_map[i] = current_cluster
        current_cluster += 1

    linkage_matrix = np.array(linkage_matrix)

    # Crear una imagen del dendrograma en memoria
    fig, ax = plt.subplots(figsize=(10, 7))
    dendro = dendrogram(
        linkage_matrix,
        labels=[chr(65 + i) for i in range(len(distance_matrix))],
        leaf_rotation=90,
        leaf_font_size=12,
        show_contracted=True,
        ax=ax
    )

    # Añadir etiquetas de distancia en el dendrograma
    for i, d, c in zip(dendro['icoord'], dendro['dcoord'], dendro['color_list']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'o', c=c)
        plt.annotate(f'{y:.2f}', (x, y), xytext=(0, -5), textcoords='offset points',
                     va='top', ha='center', fontsize=10, color=c)

    plt.title(title)
    plt.xlabel("Índice de muestra")
    plt.ylabel("Distancia")

    # Guardar la imagen en un buffer de memoria
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)

    return img_base64

def plot_dendrogram(distance_matrix, steps, title):
    # Convertir los pasos a formato de linkage para dendrograma
    linkage_matrix = []
    current_cluster = len(distance_matrix)
    cluster_map = {i: i for i in range(len(distance_matrix))}

    for i, j, dist in steps:
        linkage_matrix.append([cluster_map[i], cluster_map[j], dist, len(distance_matrix)])
        cluster_map[i] = current_cluster
        current_cluster += 1

    linkage_matrix = np.array(linkage_matrix)

    # Dibujar el dendrograma con las distancias
    plt.figure(figsize=(10, 7))
    dendro = dendrogram(
        linkage_matrix,
        labels=[chr(65 + i) for i in range(len(distance_matrix))],
        leaf_rotation=90,
        leaf_font_size=12,
        show_contracted=True
    )

    # Añadir etiquetas de distancia en el dendrograma
    for i, d, c in zip(dendro['icoord'], dendro['dcoord'], dendro['color_list']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'o', c=c)
        plt.annotate(f'{y:.2f}', (x, y), xytext=(0, -5), textcoords='offset points',
                     va='top', ha='center', fontsize=10, color=c)

    plt.title(title)
    plt.xlabel("Índice de muestra")
    plt.ylabel("Distancia")
    plt.show()
