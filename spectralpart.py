"""
Spectral Partition Algorithm implementation

"""

import numpy as np

# Creates an adjacency matrix from a file with one edge per line
# with the following format : "node1 node2"
def adjacency_matrix(filename):
    with open(filename) as f:
        lines = f.readlines()
        for i in range(len(lines)):
            lines[i] = [int(j) for j in lines[i].strip().split()]
        # The size is the maximum index of a node + 1
        size = max(max(lines, key=lambda x : max(x))) + 1
        A = [[0 for i in range(size)] for j in range(size)]
        for l in lines:
            # Keeping A symmetric (the graph is not directed)
            A[l[0]][l[1]] = 1
            A[l[1]][l[0]] = 1
        return A

# Creates a diagonal degree matrix from an adjacency matrix
def degree_matrix(adj_matrix):
    size = len(adj_matrix)
    D = [[0 for i in range(size)] for j in range(size)]
    for i in range(size):
        D[i][i] = sum(adj_matrix[i])
    return D

# Computes the laplacian matrix from the adjacency matrix and the degree matrix
def laplacian_matrix(adj_matrix, deg_matrix):
    size = len(adj_matrix)
    return [[deg_matrix[i][j] - adj_matrix[i][j] for j in range(size)] for i in range(size)]

# Generates a grouping from the laplacian matrix
def grouping_vector(lap_matrix):
    # The laplacian matrix is symmetric so we use the eigh function
    eigenvalues, eigenvectors = np.linalg.eigh(np.array(lap_matrix))
    grouping = [e for e in eigenvectors[:,1]]
    return grouping

# Returns two groups making partition of the initial matrix
def make_groups(lap_matrix, grouping_vector):
    group1, group2 = [], []
    for i in range(len(grouping_vector)):
        if grouping_vector[i] > 0:
            group1.append(i)
        else:
            group2.append(i)
    return group1, group2

# Same function, but with a hint for the size of the first group
def make_groups_sizehint(adj_matrix, lap_matrix, grouping_vector, size_1):
    l = [i for i in range(len(g))]
    l.sort(key=lambda x : grouping_vector[x])
    g11, g12 = l[:size_1], l[size_1:]    # We try to put the smallest or the biggest components of g in the first group
    g21, g22 = l[len(l)-size_1:], l[:len(l)-size_1]
    if groups_connectivity(adj_matrix, g11, g12) <= groups_connectivity(adj_matrix, g21, g22):
        return g11, g12     # And take the solution with the smallest connectivity between the groups
    else:
        return g21, g22

def groups_connectivity(adjacency_matrix, group1, group2):
    return sum(A[i][j] for j in group2 for i in group1)

if __name__ == '__main__':
    A = adjacency_matrix('p2p-Gnutella25.txt')
    print("Matrix size:", len(A))
    D = degree_matrix(A)
    L = laplacian_matrix(A, D)
    g = grouping_vector(L)

    group1, group2 = make_groups(L, g)
    print(len(group1), len(group2), groups_connectivity(A, group1, group2))
    # Outputs: 754 3285 40

    group1, group2 = make_groups_sizehint(A, L, g, int(len(A) / 2))
    print(len(group1), len(group2), groups_connectivity(A, group1, group2))
    # Outputs: 2019 2020 606

    # We test with 2 random groups of the same size to see if our algorithm gives a smaller connectivity
    group1, group2 = [i for i in range(int(len(A) / 2))], [i for i in range(int(len(A) / 2), len(A))]
    print(len(group1), len(group2), groups_connectivity(A, group1, group2))
    # Outputs: 2019 2020 8264

    # Same test with different group size
    group1, group2 = make_groups_sizehint(A, L, g, int(len(A) / 3))
    print(len(group1), len(group2), groups_connectivity(A, group1, group2))
    # Outputs: 1346 2693 2367

    group1, group2 = [i for i in range(int(len(A) / 3))], [i for i in range(int(len(A) / 3), len(A))]
    print(len(group1), len(group2), groups_connectivity(A, group1, group2))
    # Outputs: 1346 2693 13774




