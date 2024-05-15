import argparse
import os
from pgtools import analysis
    
def main():
    parser = argparse.ArgumentParser(description="Summarizes maf files statistics.\
                                     Apropriate file structure is necessary. Outputs are\
                                     placed in dir maf_comp, and in subdirectories of stats")
    parser.add_argument("results_dir",
                        help="path to the dir with maf files and mafTools results")
    parser.add_argument("gffs_dir",
                        help="path to the dir with corresponding gff files")  
    args = parser.parse_args()

    COMP_DIR = "maf_comp"
    STATS_DIR = "stats"

    results_dir = args.results_dir
    gffs_dir = args.gffs_dir

    # mafComparator summary
    comp_dir = os.path.join(results_dir, COMP_DIR)
    comp_csv = os.path.join(results_dir, "mafComp.csv")
    analysis.maf_comp_summary(comp_dir, comp_csv)    

    stats_dir = os.path.join(results_dir, STATS_DIR)
    for mafs_dir in os.listdir(stats_dir):
        # mafStats summary (for every subdirectory in stats)
        stats_csv = os.path.join(results_dir, f"{mafs_dir}_mafStats.csv")
        mafs_dir_path = os.path.join(stats_dir, mafs_dir)
        analysis.maf_stats_summary(mafs_dir_path, stats_csv)
        # cds statistics (for every subdirectory in stats)
        cds_csv = os.path.join(results_dir, f"{mafs_dir}_cds_stats.csv")
        analysis.maf_cds_summary(mafs_dir_path, gffs_dir, cds_csv)
    
if __name__ == "__main__":
    main()
