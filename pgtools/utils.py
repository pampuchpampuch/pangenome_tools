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

def reverse_coords(start: int, end: int, chr_len: int, strand: int):
    """
    ### !prob should a method of sequence class all other sequence class inherit after
    """
    rev_start = chr_len - end
    rev_end = chr_len - start + 1
    return (rev_start, rev_end, -1 * strand)

### 