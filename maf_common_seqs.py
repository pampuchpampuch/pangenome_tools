import argparse
import os
from pangenome_models.maf_parser import parse_maf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf1",
                        help="path to first maf file")
    parser.add_argument("maf2",
                        help="path to second maf file")
    parser.add_argument("outdir",
                        help="path to directory to write output to")
    
    args = parser.parse_args()

    if not os.path.exists(args.outdir): os.makedirs(args.outdir)

    maf1 = parse_maf(args.maf1)
    maf2 = parse_maf(args.maf2)

    common_chr_names = maf1.common_chr_names(maf2)

    maf1.write_to_file(os.path.join(args.outdir, args.maf1), chr_names=common_chr_names)
    maf2.write_to_file(os.path.join(args.outdir, args.maf2), chr_names=common_chr_names)
    print(maf1.chr_names)
    print(maf2.chr_names)
    print(common_chr_names)

if __name__ == "__main__":
    main()