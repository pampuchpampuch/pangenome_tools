import argparse
import os
from pgtools.maf_parser import parse_maf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf",
                        help="path to maf file")
    
    args = parser.parse_args()

    maf = parse_maf(args.maf)

    identity = maf.identity_pct()
    print(f"{args.maf}: {identity}")

if __name__ == "__main__":
    main()