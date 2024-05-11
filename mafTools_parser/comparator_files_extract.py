'''
cactus as truth (recall is intersecion/panaroo)
dir,recall,precision
'''
import argparse
import os

def parse_maf_stats(dir,file_name):
    # file = os.path.join(dir,file_name)
    lines = open(file_name).readlines()
    maf_source = lines[0].split('_')[0]
    # csv_line=[dir[:-len("_out_maf")],maf_source]
    # csv_line =[dir[:-len("_out_maf")]]
    csv_line =[dir]

    for i in range(len(lines)):
        if lines[i].strip().startswith("<homologyTests fileA"):
            csv_line.append(str(round(float(lines[i+2].split()[-1].split('"')[1]),4)))
            print(lines[i+2])
    print(csv_line)
    return csv_line

def main():
    parser = argparse.ArgumentParser(
                description="...")
    parser.add_argument("dir", metavar="dir",
                        help="path to csv file to append data to")
    parser.add_argument("comp_out", metavar="comp_out",
                        help="path to xml file with homologytest from alignathon")

    parser.add_argument("csv", metavar="csv",
                        help="path to csv file to append data to")
    args = parser.parse_args()

    # stats=parse_maf_stats(args.xml)
    '''
    csv t. że
    długość, sensitivity, specificity, precision, F-score
    '''
    # name=args.xml
    # scale=re.search("\d+",name).group(0)
    # scale=args.xml.split("_")[-1].split("/")[0]
    out_file=open(args.csv,"a")
    
    line = ",".join(parse_maf_stats(args.dir,args.comp_out))
    out_file.write(line+"\n")


if __name__=='__main__':
    main()
