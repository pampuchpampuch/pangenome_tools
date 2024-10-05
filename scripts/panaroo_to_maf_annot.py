import argparse
import os
from pgtools.gff_parser import GffCDS, Scaffold
from pgtools import panaroo_parser


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("panaroo_dir",
                        help="path to dir with panaroo output")
    parser.add_argument("gff_dir",
                        help="path to dir with gff files used to generate panaroo output")
    parser.add_argument("maf_out",
                    help="path to output maf file")
    parser.add_argument("--include_refound", action="store_true", default=False,
                        dest="include_refound",
                        help="consider genes that do not have corresponding cds in gff files")
    args = parser.parse_args()

    # gff_dict = {}
    # for filename in os.listdir(args.gff_dir): gff_dict.update(parse_gff(os.path.join(args.gff_dir,filename))[1])

    # gene_data_dir = os.path.join(args.panaroo_dir, "gene_data.csv")

    # genes_dict = parse_gene_data(gene_data_dir)

    # aln_files_dir = os.path.join(args.panaroo_dir, "aligned_gene_sequences")

    # maf_out = open(args.maf_out,"w")
    # maf_out.write("##maf version=1 scoring=N/A\n\n")

    # for aln_file in os.listdir(aln_files_dir):
    #     aln_file = os.path.join(aln_files_dir,aln_file)
    #     panaroo_aln_to_maf(aln_file, gff_dict, genes_dict, maf_out)

    # maf_out.close()
    panaroo_obj = panaroo_parser.parse_panaroo_output(args.panaroo_dir, args.gff_dir, include_refound=args.include_refound, store_annotation=True)
    panaroo_obj.to_MAF(args.maf_out)

if __name__ == "__main__":
    main()
            