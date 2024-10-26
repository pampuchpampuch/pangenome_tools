from pgtools import maf_parser
import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculates how frequently\
                                     gfa verticles are present in maf blocks (for two genomes)")
    parser.add_argument("maf",
                        help="path to the maf file")
    parser.add_argument("csv_pref",
                        help="prefix to add to csv files")
    args = parser.parse_args()
    maf = args.maf
    maf = maf_parser.parse_maf(maf, store_seqs=True)

    # maf = maf_parser.MAF(list(maf.seq_collections)[:500])
    maf = maf_parser.MAF(list(maf.seq_collections))

    maf.block_lens_csv(csv_out=args.csv_pref + "_lens.csv")
    maf.block_sizes_csv(csv_out=args.csv_pref + "_sizes.csv")


if __name__ == "__main__":
    main()