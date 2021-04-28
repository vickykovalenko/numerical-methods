import numpy as np


def pagerank(M, num_iterations: int = 100, d: float = 0.85):

    N = M.shape[1]
    v = np.random.rand(N, 1)
    M_hat = (d * M + (1 - d) / N)
    for i in range(num_iterations):
        v = M_hat @ v
        v_norm = np.linalg.norm(v, 1)
        v = v/v_norm
    return v


M = np.loadtxt('D:\\input.txt')
v = pagerank(M, 100, 0.85)
print(v)