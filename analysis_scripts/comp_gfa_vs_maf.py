import argparse
from tqdm import tqdm
from itertools import combinations
from pgtools import gfa_parser
from pgtools import maf_parser
from pgtools.utils import intersection_len, contains

'''
Jeśli chcesz sprawdzić jak to wygląda dokładniej to dla pary contigów csv:
id bloku z gfa, id bloku maf, coverage1, coverage2
Z tego też by wynikało, czy bardziej V gfa zawierają się w blokach MAF czy na odwrót

1. Find longest contigs for all genomes


1. Change coords for both models to gff
2. 
    avg_ocverages = []
    vert_lens = []
    For V in gfa verts:
        cont_lens = [[],[]]
        for block in maf blocks:
            if gfa.strands == maf.strands or gfa_strands == maf.strands[::-1]:
                inter_1 = intersection(V_1, block_1)
                inter_2 = intersection(V_2, block_2)

        Averaging for this block:
        
            #both sequences in gfa have to have equal lens
            V_len = V[0][1]-V[0][0]+1
            avg1 = sum(cont_lens[0]) / V_len
            avg2 = sum(cont_lens[1]) / V_len
            avg_cov = avg1 + avg2 / 2

            avg_coverages.append(avg_cov)
            vert_lens.append(V_len)

            # -> now we have avereged coverage of this vertex - that much of it is containe
            # -> could be further thresholded if is contained or not - to 0/1
            # it is not that important if small vertex has small coverage (but for very big ones
            # it is not worrying either - panaro has short blocks) - should we weight the mean?
    res_avg = sum(avg_coverages * vert_lens) / sum(vert_lens)


    # it is more important to weight according to contig lens probably!
'''