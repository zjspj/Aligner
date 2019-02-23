#!/usr/bin/env python
### This program takes in two command line sequence that contains ATCG in capital form and align both sequences


import numpy as np
import sys

transi_list = ["AG","GA","CT","TC"]

#penalty box
match  = 3
transi_mismatch   = 2
transver_mismatch   = -3
gap        = -2

#sequences. seq2 should always be larger than seq1
seq1 =  ''
seq2 = ''

#scoring matrix size
rows = 0
cols = 0

def create_score_matrix(rows, cols):
    '''
    Create a matrix of scores representing trial alignments of the two sequences.
    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    '''
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    # Fill the scoring matrix.
    max_score = 0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, rows):
        for j in range(1, cols):
            score = calc_score(score_matrix, i, j)
            if score > max_score:
                max_score = score
                max_pos   = (i, j)
            score_matrix[i][j] = score
    assert max_pos is not None, 'the x, y position with the highest score was not found'
    return score_matrix, max_pos

def calc_score(matrix, x, y):
    #Calculate score for a given x, y position in the scoring matrix.
    #The score is based on the up, left, and upper-left diagnol neighbors
    if seq1[x - 1] == seq2[y - 1]:
        similarity = match
    elif seq1[x - 1] + seq2[y - 1] in transi_list:
        similarity = transi_mismatch
    else:
        similarity = transver_mismatch
    diag_score = matrix[x - 1][y - 1] + similarity
    up_score   = matrix[x - 1][y] + gap
    left_score = matrix[x][y - 1] + gap
    return max(0, diag_score, up_score, left_score)

def print_matrix(matrix):
    '''
    Print the scoring matrix.
    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    print(np.matrix(matrix).T)


def path_str(matrix, start_pos):
    """
    This is the main path finding program. It will also add remaining sequence after the string if needed
    It takes the scoring matrix and starting position
    It returns a list that contains three string, first string is the first string entered
    second string is the second string entered
    third string is the symbol for all of the match and mismatch
    """
    print(seq1, seq2)
    path_str_list = ["","",""]
    pos = start_pos
    print(pos)
    path_str_list = append_path_str(pos, path_str_list, "m")
    pos_pair = path_find(matrix, pos)
    while matrix[pos_pair[0][0]][pos_pair[0][1]] != 0 and pos_pair[0][0] > 0 and pos_pair[0][1] > 0:
        pos, pos_status = pos_pair[0], pos_pair[1]
        path_str_list = append_path_str(pos, path_str_list, pos_status)
        pos_pair = path_find(matrix, pos)
    if pos_pair[0][0] > 1:
        for extra_base in range(pos_pair[0][0] - 1):
            path_str_list[0] += seq1[extra_base]
            path_str_list[2] += "/"
            path_str_list[1] += " "
    elif pos_pair[0][1] > 1:
        for extra_base in pos_pair[0][1] - 1:
            path_str_list[0] += " "
            path_str_list[2] += "/"
            path_str_list[1] += seq2[extra_base]
    return path_str_list

def append_path_str(position, str_list, status):
    """
    This function takes in the position of the cell, string list, and the status of the previous cell as argument
    then append the information of the current cell into the string list which contains all three lines
    match is indicated by "*"
    transition is indicated by ":"
    transversion is indicated by "."
    gap is indicated by a " "
    """
    x = seq1[position[0] - 1]
    y = seq2[position[1] - 1]
    print(position, status)
    print(x,y)
    print(str_list)
    print("")
    if status == "le":
        str_list[0] += x
        str_list[2] += " "
        str_list[1] += "-"

    elif status == "up":
        str_list[0] += "-"
        str_list[2] += " "
        str_list[1] += y

    elif status == "m":
        str_list[0] += x
        if x == y:
            str_list[2] += "*"
        elif x + y in transi_list:
            str_list[2] += ":"
        else:
            str_list[2] += "."
        str_list[1] += y
    return str_list

def path_find(matrix, position):
    """
    This function takes in the scoring matrix and the position of the current cell
    It output a list that contains the coordinate of the next cell and the direction of the movement
    """
    left = (position[0]-1, position[1])
    up = (position[0], position[1]-1)
    diag = (position[0]-1, position[1]-1)
    largest = max(matrix[left[0]][left[1]], matrix[up[0]][up[1]], matrix[diag[0]][diag[1]])
    if matrix[left[0]][left[1]] == largest:
        return [left, "le"]
    elif matrix[up[0]][up[1]] == largest:
        return [up, "up"]
    elif matrix[diag[0]][diag[1]] == largest:
        return [diag, "m"]

if __name__ == '__main__':
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    for base in seq1 + seq2:
        if base != "A" and base != "T" and base != "C" and base != "G":
            raise Exception("sequence must only contain ATCG upper case letters")
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    print(seq1,seq2)
    score_matrix, start_pos = create_score_matrix(rows, cols)
    path_str_list = path_str(score_matrix, start_pos)
    for path_str in path_str_list:
        print(path_str[::-1])
    print_matrix(score_matrix)
