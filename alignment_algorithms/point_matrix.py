def point_matrix(seq1, seq2):
    matrix = []
    for i in range(len(seq1)):
        row = []
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                row.append(1)
            else:
                row.append(0)
        matrix.append(row)
    return matrix
