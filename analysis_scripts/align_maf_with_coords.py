import argparse
from pgtools import maf_parser, gff_parser
from Bio.Seq import Seq
from tqdm import tqdm

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf_in",
                        help="dir to pseudo maf with coords")
    parser.add_argument("gff_dir")
    parser.add_argument("maf_out")
    args = parser.parse_args()

    maf_dir = args.maf_in
    maf_obj = maf_parser.parse_maf(maf_dir)
    all_scaffs = gff_parser.scaffolds_from_GFFs_dir(args.gff_dir)



    new_blocks = []

    for block in maf_obj.seq_collections:
        new_block = set()
        for seq in block.sequences:

            # print(all_scaffs[seq.get_genome_name()][seq.get_contig_name()][:3])
            # seq_scaff = all_scaffs[seq.get_genome_name()][seq.get_contig_name()][seq.start:seq.end]
            # if not seq.get_genome_name() in all_scaffs:
            #     continue
            seq_scaff = all_scaffs[seq.seq_name]
            if seq.strand < 0:
                seq_scaff = Seq(seq_scaff)
                seq_scaff = str(seq_scaff.reverse_complement())
                # seq_scaff = seq_scaff
            seq_scaff = seq_scaff[seq.start:seq.end]
            seq.seq = seq_scaff
            new_block.add(seq)
            # seq_0 = seq.seq.replace("-","")
            # assert len(seq_scaff) == len(seq_0), f"{seq_0}, {seq_scaff}"
        new_blocks.append(maf_parser.SyntenyBlock(new_block))
    # break
    
    new_maf = maf_parser.MAF(new_blocks)

    maf_out = args.maf_out
    res_maf = open(maf_out, "w")
    for block in tqdm(list(new_maf.seq_collections)):
        res_maf.write(block.to_MAF_block(align=True))
        res_maf.flush()
    res_maf.close()

if __name__ == "__main__":
    main()