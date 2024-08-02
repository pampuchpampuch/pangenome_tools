import numpy as np

def read_phylip_matrix(file):
    res_matrix = None
    n_genomes = None
    with open(file) as f:
        # next(f)
        n_genomes = int(next(f))
        res_matrix = np.array([])
        for line in f:
            res_matrix = np.append(res_matrix, list(map(lambda x: float(x), line.split()[1:])), axis=0)
    return res_matrix.reshape((n_genomes, n_genomes))

    # return np.tril(res_matrix)

def mean_dist(phylip_M):
    n_genomes = phylip_M.shape[0]
    return(np.sum(phylip_M)/(n_genomes*n_genomes - n_genomes)/2)

phy_M = read_phylip_matrix("../phylonium_out/phylonium.out")
print(mean_dist(phy_M))