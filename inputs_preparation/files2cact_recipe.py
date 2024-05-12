import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dir",
                        help="path to dir with gz fasta files")
    parser.add_argument("out_file",
                        help="path to output file - cactus recipe")
    
    args = parser.parse_args()

    out_file = open(args.out_file, "w")

    for filename in os.listdir(args.dir):
        out_file.write(filename.split(".fa")[0]+"\t"+args.dir+"/"+filename+"\n")
if __name__ == "__main__":
    main()