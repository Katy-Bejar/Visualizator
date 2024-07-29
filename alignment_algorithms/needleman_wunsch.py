import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import io
import base64

def needleman_wunsch_all(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)

    # Inicialización de la matriz de puntuación
    score_matrix = np.zeros((n+1, m+1))
    for i in range(n+1):
        score_matrix[i][0] = i * gap
    for j in range(m+1):
        score_matrix[0][j] = j * gap

    # Llenado de la matriz de puntuación
    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete_score = score_matrix[i-1][j] + gap
            insert_score = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(match_score, delete_score, insert_score)

    # Traza hacia atrás para encontrar todos los alineamientos
    def traceback(i, j, align1, align2, colors):
        if i == 0 and j == 0:
            alignments.append((align1, align2, colors))
            return
        if i > 0 and score_matrix[i][j] == score_matrix[i-1][j] + gap:
            traceback(i-1, j, seq1[i-1] + align1, "-" + align2, ['gap'] + colors)
        if j > 0 and score_matrix[i][j] == score_matrix[i][j-1] + gap:
            traceback(i, j-1, "-" + align1, seq2[j-1] + align2, ['gap'] + colors)
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            color = 'match' if seq1[i-1] == seq2[j-1] else 'mismatch'
            traceback(i-1, j-1, seq1[i-1] + align1, seq2[j-1] + align2, [color] + colors)

    alignments = []
    traceback(n, m, "", "", [])

    return alignments, score_matrix, match, mismatch, gap

''' 
def plot_score_matrix(seq1, seq2, score_matrix, match, mismatch, gap):
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

    # Draw arrows
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            if score_matrix[i, j] == score_matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                ax.add_patch(patches.FancyArrowPatch((j-1, i-1), (j, i), color='black', arrowstyle='->'))
            if score_matrix[i, j] == score_matrix[i-1, j] + gap:
                ax.add_patch(patches.FancyArrowPatch((j, i-1), (j, i), color='black', arrowstyle='->'))
            if score_matrix[i, j] == score_matrix[i, j-1] + gap:
                ax.add_patch(patches.FancyArrowPatch((j-1, i), (j, i), color='black', arrowstyle='->'))

    plt.title('Matriz de Puntuación de Needleman-Wunsch')
    plt.xlabel('Secuencia 2')
    plt.ylabel('Secuencia 1')

    # Convert the plot to a PNG image and return as base64
    img_io = io.BytesIO()
    plt.savefig(img_io, format='png')
    img_io.seek(0)
    img_base64 = base64.b64encode(img_io.getvalue()).decode('utf8')
    plt.close(fig)
    return img_base64
'''