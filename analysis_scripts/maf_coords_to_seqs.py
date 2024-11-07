from pgtools import panaroo_parser
# from pgtools.gff_parser import parse_GFFs_dir, Gff, parse_gff
from pgtools.utils import intersection_len
from pgtools import maf_parser, gff_parser
from pgtools.pangenome import Pangenome
from itertools import combinations
from copy import deepcopy
from pgtools.pangenome import Pangenome
from pgtools.maf_parser import SyntenyBlock, MAFseq
import os
import argparse
from tqdm import tqdm
import numpy as np 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf")
    parser.add_argument("gff_dir")
    parser.add_argument("maf_out")
    args = parser.parse_args()

    # panaroo_obj = panaroo_parser.parse_panaroo_output(args.panaroo_dir, args.gff_dir)
    # new_blocks = maf_new_blocks(panaroo_obj, args.maf_out)
    # edges = panaroo_obj.get_cluster_edges()
    # print(edges[:4])
    # clusters_dict = panaroo_obj.get_seq_collections_dict()

    maf = maf_parser.parse_maf(args.maf)
    # for seq_coll in maf.seq_collections:
    #     print(seq_coll.to_MAF_block())
    #     break
    scaffolds = gff_parser.scaffolds_from_GFFs_dir(args.gff_dir)
    print(scaffolds)

    

if __name__ == "__main__":
    main()