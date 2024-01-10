import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import bmat, coo_matrix
from collision import *

def A_matrix(theta, m=1, l=1):
    N = len(theta)
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i == j:
                A[i, j] = m*l**2
            else:
                A[i, j] = m*l**2*(N-i)/(N-i+1)*np.cos(theta[j] - theta[i])

    return A

def M(index, THETA):
    return A_matrix(THETA[index])

def HcSolid(t, n, theta, l):
    """
    t and n : respectively tangeant and normal of collision
    theta : array of angular positions of rods for a solid (theta0, ... )
    l : length of rods
    """
    # Calculate the matrix Pc containing coordinates of t and n
    Pc = np.array([t, n]).T  # Transpose to make it a 2x2 matrix    
    Jacobian = np.vstack((l * np.cos(theta), l * np.sin(theta)))
    return Pc @ Jacobian

def H_c(nc, ndl, nsc, collisions, index_mapping, THETA):
    """
    nc : nombre de point de contact
    nd : nombre de degré de liberté (theta)
    nsc : nombre de solide en contact
    collisions : array of solids in collision, example : [(i1, j1), (i2, j2), ...]
    """
    block_size = (2, ndl)
    num_blocks = (nc, nsc)
    blocks = [[coo_matrix(np.zeros((block_size[0], block_size[1]))) for _ in range(num_blocks[1])] for _ in range(num_blocks[0])]

    for s, c in enumerate(collisions): # collisions is a list of object collision each object contain all important informations
        i, j = c.solids
        blocks[s][index_mapping[i]] = - coo_matrix(HcSolid(c.t, c.n, THETA[i], 1))
        blocks[s][index_mapping[j]] =   coo_matrix(HcSolid(c.t, c.n, THETA[j], 1))

    Hc = bmat(blocks)

    return Hc

def M_c(ndl, nsc, index_mapping, THETA):
    block_size = (ndl, ndl)
    num_blocks = (nsc, nsc)
    blocks = [[coo_matrix(np.zeros((block_size[0], block_size[1]))) for _ in range(num_blocks[1])] for _ in range(num_blocks[0])]

    solidsInContact = list(index_mapping)
    for s in range(len(solidsInContact)):
        blocks[s][s] = coo_matrix(M(index=solidsInContact[s], THETA=THETA))

    Mc = bmat(blocks)

    return Mc

def test_Hc():
    # Generate some sample collision objects
    collision1 = Collision(np.array([1.0, 0.0]), np.array([0.0, 1.0]), (0, 5), (2, 3))
    collision2 = Collision(np.array([0.0, 1.0]), np.array([1.0, 0.0]), (3, 2), (3, 4))
    collision3 = Collision(np.array([0.0, 1.0]), np.array([1.0, 0.0]), (1, 6), (3, 4))
    collisions = [collision1, collision2, collision3]

    # Create a sample index_mapping dictionary
    index_mapping = {0: 0, 1: 1, 2: 2, 3: 3, 5: 4, 6: 5}

    # Call Hc function with sample collisions and index_mapping
    return H_c(nc=3, ndl=10, nsc=6, collisions=collisions, index_mapping=index_mapping, THETA=[np.array([0.1]*10)]*7)

def test_Mc():
    # Create a sample index_mapping dictionary
    index_mapping = {0: 0, 1: 1, 2: 2, 3: 3, 5: 4, 6: 5}

    # Call Hc function with sample collisions and index_mapping
    return M_c(ndl=10, nsc=6, index_mapping=index_mapping, THETA=[np.array([0.1]*10)]*7)

def test_sparse_matrix():
    # Assuming you have a large matrix H with known dimensions
    # For example, let's say H is 4x4 and created by 2x2 blocks
    block_size = 2
    num_blocks = 10
    H_dim = block_size * num_blocks

    # Create a random large matrix H
    H = np.random.rand(H_dim, H_dim)

    # Create a block matrix using scipy.sparse.bmat
    blocks = [[coo_matrix(np.zeros((block_size, block_size))) for _ in range(num_blocks)] for _ in range(num_blocks)]
    # Create a block of ones with dimensions block_size x block_size
    ones_block = np.ones((block_size, block_size))

    # Replace the upper-left block with the block of ones
    blocks[0][0] = blocks[-1][-1] = coo_matrix(ones_block)


    # Construct the block matrix
    H_sparse = bmat(blocks)

    print(H_sparse)
    # Plot the entire matrix
    plt.figure(figsize=(8, 8))
    plt.spy(H_sparse, markersize=10)
    plt.title('Sparse Matrix Visualization')
    plt.show()

    # Example usage:
    t = np.array([1.0, 0.0])  # Replace with your actual tangent vector
    n = np.array([0.0, 1.0])  # Replace with your actual normal vector
    theta = np.array([0.1, 0.5, 1.0])  # Replace with your actual array of angular positions
    l = 1.0  # Replace with the actual length of the rods

    result_matrix = HcSolid(t, n, theta, l)
    print(result_matrix)

if __name__ == "__main__":
    # Tests
    H_sparse = test_Hc()
    print(H_sparse)

    plt.figure(figsize=(8, 8))
    plt.spy(H_sparse, markersize=10)
    plt.title('Hc')
    plt.show()