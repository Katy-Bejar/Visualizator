import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def smith_waterman_all(seq1, seq2, match=1, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)

    # Inicialización de la matriz de puntuación
    score_matrix = np.zeros((n+1, m+1))
    max_score = 0
    max_positions = []

    # Llenado de la matriz de puntuación y búsqueda de los máximos puntajes
    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete_score = score_matrix[i-1][j] + gap
            insert_score = score_matrix[i][j-1] + gap
            score_matrix[i][j] = max(0, match_score, delete_score, insert_score)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_positions = [(i, j)]
            elif score_matrix[i][j] == max_score:
                max_positions.append((i, j))

    # Función para obtener todos los alineamientos locales posibles
    def traceback(i, j, align1, align2):
        alignments = []
        while score_matrix[i][j] > 0:
            current_score = score_matrix[i][j]
            if current_score == 0:
                break
            if current_score == score_matrix[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif current_score == score_matrix[i-1][j] + gap:
                align1 = seq1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            elif current_score == score_matrix[i][j-1] + gap:
                align1 = "-" + align1
                align2 = seq2[j-1] + align2
                j -= 1
        alignments.append((align1, align2))
        return alignments

    all_alignments = []
    for pos in max_positions:
        i, j = pos
        alignments = traceback(i, j, "", "")
        all_alignments.extend(alignments)

    return all_alignments, score_matrix