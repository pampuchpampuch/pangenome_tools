"""
As a baseline use Panaroo and do not map - parse alignment files / gml:
for each full cluster all annotations (optionally with contig names)

Mapping for other models

for each cluster:
    for each seuence:
        get list of mapped annotations

then for each full cluster: list of annotations

then different approaches could be used:
1. Simplest one: every annotation that was in a full cluster is core
2. Less naive one: contents have to be similiar - if multiple annotations mapped then a
    corresponding ponaroo cluser should be found with a significant annotation overlap

How to compare - more clearly defined clusters could get higher weight or be
considered separately - more clearly defined ie proteins with common name or consistenmt annotation id

stats:
common / panaroo and common / the other model

Get csv such that:

seq_id, cluster_id, [mapped annotations], core status


"""

from pgtools import panaroo_parser
from pgtools.gff_parser import parse_gff, parse_GFFs_dir, Pangenome_Gffs
from pgtools.gfa_parser import parse_gfa1
from pgtools.maf_parser import parse_maf
from pgtools.utils import intersection_len, contains
from pgtools.pangenome import Pangenome
import os
import argparse
from Bio.Seq import Seq

def main():
    parser = argparse.ArgumentParser(description="Writes annotation info into csv")
    parser.add_argument("maf",
                        help="path to the maf file")
    # parser.add_argument("gff_dir",
    #                     help="dir with gff files")
    # parser.add_argument("seqs_summary_out",
    #                     help="csv output file")
    # parser.add_argument("annots_out",
    #                     help="csv output file")
    # parser.add_argument("--overlap_threshold", default=0.7, dest="overlap_threshold", type=float,
    #                     help="section of sequences that have to overlap to consider them mapped to annotation")
    args = parser.parse_args()

    maf = parse_maf(args.maf)
    # maf.annotations_to_csv(args.gff_dir, csv_seqs_name=args.seqs_summary_out, csv_annots_name=args.annots_out, overlap_threshold=args.overlap_threshold)
    triplets = maf.get_panaroo_cluster_triplets()
    print(triplets)

if __name__ == "__main__":
    main()