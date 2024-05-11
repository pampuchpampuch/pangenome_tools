import argparse
import re
#chce mieć jaki to xml (np MRCA neutral i z jakiego maf (mammals.4500))
#pętla bo wszytskch plikach tak żeby potem w jednym pliku csv były kolumny dla neurtral,
#repetetive itd, i wiersze dla każdego z plików 100, 500 itd.

'''
w jednym pliku są 2 fileA= i fileB=
w tym na górze average to sensitivity, a na dole specifity
'''

def get_recall(xml):
    stats={}
    lines=open(xml, 'r').readlines()

    for i in range(len(lines)):
        if lines[i].strip().startswith("<homologyTests fileA"):
            if "sensitivity" in stats:
                stats["specificity"] = float(lines[i+2].split()[4][9:-3])
                TP = int(lines[i+2].split()[2][11:-1])
                print(lines[i])
                print(lines[i+2])
                print("TP",TP)
                print("spec",stats["specificity"])

            else:
                stats["sensitivity"] = float(lines[i+2].split()[4][9:-3])
                FP = int(lines[i+2].split()[3][12:-1])
                print(lines[i])
                print(lines[i+2])
                print("FP",FP)

    stats["precision"] = TP/(TP+FP)
    prec=stats["precision"]
    sens=stats["sensitivity"]
    stats["Fscore"] =  2*prec*sens/(prec+sens)

    return stats

def main():
    parser = argparse.ArgumentParser(
                description="...")
    parser.add_argument("xml", metavar="xml",
                        help="path to xml file with homologytest from alignathon")
    parser.add_argument("csv", metavar="csv",
                        help="path to csv file to append data to")
    args = parser.parse_args()

    stats=get_recall(args.xml)
    '''
    csv t. że
    długość, sensitivity, specificity, precision, F-score
    '''
    name=args.xml
    scale=re.search("\d+",name).group(0)
    # scale=args.xml.split("_")[-1].split("/")[0]
    out_file=open(args.csv,"a")
    
    line= scale+","+str(stats["sensitivity"])+","+str(stats["specificity"])+","+str(stats["precision"])+","+str(stats["Fscore"])+"\n"
    out_file.write(line)


if __name__=='__main__':
    main()
