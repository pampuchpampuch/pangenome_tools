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

def intersection_len(s1_coords, s2_coords):
    """
    seq2 is assumed to be from MAF and seq1 from gfa (checking how well
    seq1 fits into seq2)
    """
    s1_start, s1_end = s1_coords
    s2_start, s2_end = s2_coords
    # if (seq1_coords[0] >= seq2_coords[0]) and (seq1_coords[1] <= seq2_coords[1])
    if s1_end < s2_start or s1_start > s2_end:
        return 0

    if s1_start < s2_start:
        inter_start = s2_start
    else:
        inter_start = s1_start    

    if s1_end > s2_end:
        inter_end = s2_end
    else:
        inter_end = s1_end
        
    return inter_end - inter_start + 1

def contains(s1_coords, s2_coords, threshold = 0.7):
    """
    Assumes s1 is from gfa and s2 from maf. The goal is to check if
    gfa vertex is contained in maf block. Threshold is a fraction of 
    gfa seq that is contained in maf block sequence
    """
    inter_len = intersection_len(s1_coords, s2_coords)
    return inter_len >= threshold