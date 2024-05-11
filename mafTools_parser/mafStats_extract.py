'''
columns:
dir, cact/panaroo, char_n, gap_n, columns_n, blocks_n, ave_area, max_area, total_area, ave_degree, max_degree, ave_seq, max_seq, n_seq
char pct is char_n/total_area
'''
import argparse
import os

def parse_maf_stats(dir,file_name):
    file = os.path.join(dir,file_name)
    lines = open(file).readlines()
    maf_source = lines[0].split('.')[0].split('_')[0]
    csv_line=[dir,maf_source]

    for i in range(12,24):
        if i == 15: continue
        else: csv_line.append(lines[i].split(":")[1].split()[0])

    csv_line.append(lines[25].split()[0])
    return csv_line

def main():
    parser = argparse.ArgumentParser(
                description="...")
    parser.add_argument("dir", metavar="dir",
                        help="path to csv file to append data to")
    parser.add_argument("stats", metavar="stats",
                        help="path mafStats file")

    parser.add_argument("csv", metavar="csv",
                        help="path to csv file to append data to")
    args = parser.parse_args()

    # stats=parse_maf_stats(args.xml)
    # name=args.xml
    # scale=re.search("\d+",name).group(0)
    # scale=args.xml.split("_")[-1].split("/")[0]
    out_file=open(args.csv,"a")
    
    line = ",".join(parse_maf_stats(args.dir,args.stats))
    out_file.write(line+"\n")


if __name__=='__main__':
    main()
