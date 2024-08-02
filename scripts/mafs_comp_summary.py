import argparse
import os
from pgtools import analysis
    
def maf_comp_summary(maf_comp_dir, out_csv):
    """
    Summarizes mafComparator outputs for files in subdirectories
    (dir name is saved in csv)
    
    Parameters:
        maf_comp_dir: str
            directory with directories containing mafComparator outputs
    """
    csv_file = open(out_csv,"w")
    csv_file.write("dataset,recall,precision\n")

    # for dataset in os.listdir(maf_comp_dir):
    for f in os.listdir(maf_comp_dir):
        #file names have to follow strict naming convention
        if not f.endswith(".xml"): continue
        f_path = os.path.join(maf_comp_dir, f)
        print(f_path)
        csv_line=",".join(analysis.parse_maf_comp(f_path))+"\n"
        csv_file.write(csv_line)

def main():
    parser = argparse.ArgumentParser(description="Summarizes maf files statistics.\
                                     Apropriate file structure is necessary. Outputs are\
                                     placed in dir maf_comp, and in subdirectories of stats")
    parser.add_argument("comp_dir",
                        help="path to the dir with maf files and mafTools results")
    # parser.add_argument("gffs_dir",
    #                     help="path to the dir with corresponding gff files")  
    args = parser.parse_args()



    # results_dir = args.results_dir
    # gffs_dir = args.gffs_dir

    # mafComparator summary
    comp_dir = args.comp_dir
    comp_csv = os.path.join(comp_dir, "mafComp.csv")
    maf_comp_summary(comp_dir, comp_csv)    
    
if __name__ == "__main__":
    main()
