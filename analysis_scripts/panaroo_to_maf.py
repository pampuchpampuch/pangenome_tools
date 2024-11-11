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
    parser.add_argument("panaroo_dir",
                        help="path to panaroo output")
    parser.add_argument("gff_dir",
                        help="dir with gff files")
    parser.add_argument("maf_out",
                        help="maf output file")

    args = parser.parse_args()

    panaroo_obj = panaroo_parser.parse_panaroo_output(args.panaroo_dir, args.gff_dir)
    # panaroo_obj.assigned_annotations_to_csv(args.gff_dir, csv_seqs_name=args.csv_out_seqs, csv_annots_name=args.csv_out_annots)
    panaroo_obj.to_MAF(args.maf_out)
if __name__ == "__main__":
    main()