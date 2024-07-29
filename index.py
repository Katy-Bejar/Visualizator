from flask import Flask, request, jsonify, render_template
from alignment_algorithms.needleman_wunsch import needleman_wunsch_all
from alignment_algorithms.smith_waterman import smith_waterman_all
import io
import base64
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from cluster_algorithms.cluster import single_linkage_clustering
from cluster_algorithms.cluster import create_linkage_matrix
from cluster_algorithms.cluster import plot_dendrogram
from estructura.estructura_secundaria import predict_secondary_structure, traceback, plot_structure
import networkx as nx
from io import BytesIO
from scipy.cluster.hierarchy import dendrogram, linkage

app = Flask(__name__)

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/La_Secuencia')
def la_secuencia():
    return render_template('Secuencia.html')

@app.route('/El_Alineamiento_2_sec')
def el_alineamiento_2_sec():
    return render_template('Alineamiento_2_sec.html')

@app.route('/El_Alineamiento_multiple')
def el_alineamiento_multiple():
    return render_template('Alineamiento_multiple.html')

@app.route('/El_Cluster')
def el_cluster():
    return render_template('Clusterizacion.html')
@app.route('/analyze_clustering', methods=['POST'])
def analyze_clustering():
    data = request.get_json()
    matrix = np.array(data.get('matrix'))
    method = data.get('method')

    if matrix is None or method not in ['single', 'complete', 'average']:
        return jsonify({'error': 'Datos inválidos.'})

    # Realizar el clustering con el método especificado
    Z = linkage(matrix, method=method)

    # Crear el dendrograma
    plt.figure(figsize=(10, 7))
    dendro = dendrogram(Z, labels=[chr(65 + i) for i in range(len(matrix))], leaf_rotation=90, leaf_font_size=12)

    # Añadir etiquetas de distancia en el dendrograma
    for i, d, c in zip(dendro['icoord'], dendro['dcoord'], dendro['color_list']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'o', c=c)
        plt.annotate(f'{y:.2f}', (x, y), xytext=(0, -5), textcoords='offset points',
                     va='top', ha='center', fontsize=10, color=c)

    img_io = io.BytesIO()
    plt.savefig(img_io, format='png')
    img_io.seek(0)
    img_base64 = base64.b64encode(img_io.getvalue()).decode('utf8')
    plt.close()

    return jsonify({'dendrogram_image': img_base64})



#---------------------------FILOGENIA--------------------------------#
@app.route('/La_Filogenia')
def la_filogenia():
    return render_template('Filogenia.html')


#--------------------------------------------------------------------#



@app.route('/La_Estructura_secundaria')
def la_estructura_secundaria():
    return render_template('Estructura_secundaria.html')
@app.route('/analyze_secondary_structure', methods=['POST'])
def analyze_secondary_structure():
    data = request.get_json()
    sequence = data['sequence']

    # Validar la secuencia aquí (si es necesario)
    if not sequence or any(c not in "ACGT" for c in sequence):
        return jsonify({'error': 'Secuencia no válida. Asegúrate de que solo contiene A, C, G, T.'})

    alpha_dict = {"CG": -1, "GC": -1, "AU": -1, "UA": -1, "GU": -1, "UG": -1}

    energy_matrix, min_energy_score = predict_secondary_structure(sequence, alpha_dict)
    traceback_structure, paired_positions, structures = traceback(sequence, energy_matrix, alpha_dict)

    # Crear imagen y guardarla en un buffer de bytes
    img = BytesIO()
    plot_structure(sequence, paired_positions, structures, img, node_color="#00BFFF", edge_color="#D3D3D3", highlight_color="#FF4500")
    
    # Codificar la imagen en base64
    img.seek(0)
    img_base64 = base64.b64encode(img.getvalue()).decode('utf-8')

    return jsonify({'structureImage': img_base64})



@app.route('/analyze_alignment', methods=['POST'])
def analyze_alignment():
    data = request.get_json()
    sequence1 = data.get('sequence1')
    sequence2 = data.get('sequence2')
    gap = data.get('gap', -2)
    algorithm = data.get('algorithm')

    if not sequence1 or not sequence2:
        return jsonify({'error': 'Las secuencias no pueden estar vacías.'})

    if algorithm == 'needleman_wunsch':
        alignments, score_matrix, match, mismatch, gap = needleman_wunsch_all(sequence1, sequence2, gap=gap)
    elif algorithm == 'smith_waterman':
        alignments, score_matrix = smith_waterman_all(sequence1, sequence2, gap=gap)
        match = 1  # Cambiar según tu configuración
        mismatch = -1  # Cambiar según tu configuración
    else:
        return jsonify({'error': 'Algoritmo no válido seleccionado.'})
    
    result = []
    for align1, align2, colors in alignments:
        alignment = []
        for base1, base2, color in zip(align1, align2, colors):
            alignment.append({'base1': base1, 'base2': base2, 'color': color})
        result.append(alignment)

    # Generate score matrix image as base64
    matrix_image_base64 = plot_score_matrix(sequence1, sequence2, score_matrix, match, mismatch, gap, algorithm)

    return jsonify({'alignments': result, 'matrix_image': matrix_image_base64})



def plot_score_matrix(seq1, seq2, score_matrix, match, mismatch, gap, algorithm):
    fig, ax = plt.subplots(figsize=(10, 10))
    cax = ax.matshow(score_matrix, cmap='Blues')
    fig.colorbar(cax)

    ax.set_xticks(np.arange(len(seq2) + 1))
    ax.set_yticks(np.arange(len(seq1) + 1))
    ax.set_xticklabels([''] + list(seq2))
    ax.set_yticklabels([''] + list(seq1))

    for i in range(len(seq1) + 1):
        for j in range(len(seq2) + 1):
            c = score_matrix[i, j]
            ax.text(j, i, str(int(c)), va='center', ha='center', color='black')

    if algorithm == 'needleman_wunsch':
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                if score_matrix[i, j] == score_matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                    ax.add_patch(patches.FancyArrowPatch((j-1, i-1), (j, i), color='black', arrowstyle='->'))
                if score_matrix[i, j] == score_matrix[i-1, j] + gap:
                    ax.add_patch(patches.FancyArrowPatch((j, i-1), (j, i), color='black', arrowstyle='->'))
                if score_matrix[i, j] == score_matrix[i, j-1] + gap:
                    ax.add_patch(patches.FancyArrowPatch((j-1, i), (j, i), color='black', arrowstyle='->'))

    if algorithm == 'smith_waterman':
        for i in range(1, len(seq1) + 1):
            for j in range(1, len(seq2) + 1):
                if score_matrix[i, j] > 0:
                    if score_matrix[i, j] == score_matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                        ax.add_patch(patches.FancyArrowPatch((j-1, i-1), (j, i), color='red', arrowstyle='->'))
                    if score_matrix[i, j] == score_matrix[i-1, j] + gap:
                        ax.add_patch(patches.FancyArrowPatch((j, i-1), (j, i), color='black', arrowstyle='->'))
                    if score_matrix[i, j] == score_matrix[i, j-1] + gap:
                        ax.add_patch(patches.FancyArrowPatch((j-1, i), (j, i), color='black', arrowstyle='->'))
                    if score_matrix[i, j] == 0:
                        ax.add_patch(patches.FancyArrowPatch((j, i), (j, i), color='black', arrowstyle='->'))


    plt.title(f'Matriz de Puntuación de {algorithm.replace("_", " ").title()}')
    plt.xlabel('Secuencia 2')
    plt.ylabel('Secuencia 1')

    img_io = io.BytesIO()
    plt.savefig(img_io, format='png')
    img_io.seek(0)
    img_base64 = base64.b64encode(img_io.getvalue()).decode('utf8')
    plt.close(fig)
    return img_base64

if __name__ == '__main__':
    app.run(debug=True)
