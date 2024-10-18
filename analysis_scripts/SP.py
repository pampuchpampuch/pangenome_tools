import argparse
from Bio import AlignIO
import itertools

def score(pair):
    if "-" in pair:
        return -1
    elif pair[0]==pair[1]:
        return 1
    else:
        return -2

def evaluate_maf(maf):
    gaps_pcts=[]
    SP=[]
    for MSA in AlignIO.parse(maf, "maf"):
        MSA_SP=0
        n_rows=len(MSA)
        for i in range(MSA.get_alignment_length()):
            col=MSA[:,i]
            n_gaps=col.count("-")
            gaps_pcts.append(n_gaps/n_rows)

            for pair in itertools.combinations(col,2):
                MSA_SP+=score(pair)

        SP.append(MSA_SP)

    return(gaps_pcts,SP)

def main():
    parser = argparse.ArgumentParser(
                description="...")
    parser.add_argument("maf", metavar="maf",
                        help="path to maf file")
    args = parser.parse_args()

    gaps_pcts,SP=evaluate_maf(args.maf)
    print(gaps_pcts)
    print(SP)

if __name__=='__main__':
    main()
