import argparse
import sys

class MAF:
    '''
    represents a synteny file as a set of synteny blocks
    '''
    def __init__(self, synteny_blocks):
        self.synteny_blocks = synteny_blocks
        chr_names = []
        for syn_block in synteny_blocks:
            chr_names += [maf_seq.seq_name for maf_seq in syn_block.block_seqs]
        self.chr_names = list(set(chr_names))
    
    def add_synteny_block(self, syn_block):
        self.synteny_blocks.append(syn_block)

    def write_to_file(self, out_file, chr_names='all'):
        out_file = open(out_file,"w")
        out_file.write("##maf version=1 scoring=N/A\n\n")

        if chr_names=='all':
            for block in self.synteny_blocks:
                out_file.write(block.MAF_repr())
        
        else:
             for block in self.synteny_blocks:
                out_file.write(block.MAF_repr(chr_names))           

        out_file.close()

    def common_chr_names(self, other_maf):
        return set(self.chr_names).intersection(set(other_maf.chr_names))
    
    def get_sequences(self):
        sequences = {}
        for block in self.synteny_blocks:
            for seq in block.block_seqs:
                if not seq.seq_name in sequences:
                    sequences[seq.seq_name] = [seq]
                else:
                    sequences[seq.seq_name].append(seq)
        return sequences

class SyntenyBlock:
    '''
    represents synteny block as a set of sequences with an id
    from coords file
    '''

    def __init__(self, block_seqs, aligned=False):
        # consider adding id?
        self.block_seqs=block_seqs
        self.aligned = aligned

    def add_block_seq(self,block_seq):
        self.block_seqs.append(block_seq)

    def MAF_repr(self, chr_names = None):
        '''
        return a string representing block in MAF format
        '''
        if self.aligned:

            maf_str = "a\n"

            if chr_names:
                for maf_seq in self.block_seqs:
                    if maf_seq.seq_name in chr_names:
                        maf_str += maf_seq.MAF_repr() 

            else:
                for maf_seq in self.block_seqs:
                    maf_str += maf_seq.MAF_repr() 

            maf_str += '\n'       

        else:
            sys.exit("MAF block alignment to be implemented")
        
        return maf_str

class MAFseq:
    '''
    Represents sequence in maf block
    '''
    def __init__(self, chr_name, start, end, strand, chr_size, seq):
        self.seq_name = chr_name
        self.start = start
        self.end = end
        #self.seq_len = int(end) - int(start)
        self.strand = strand
        self.chr_size = chr_size
        self.seq = seq

    def __len__(self):
        return (self.end - self.start)
    
    def MAF_repr(self):
        '''
        returns string representing seguence in MAF format
        '''
        
        strand_sign = "+" if self.strand > 0 else "-"

        s_line=f"s\t{self.seq_name}\t{self.start}\t{len(self)}\t{strand_sign}\t{self.chr_size}\t{self.seq}\n"

        return(s_line)

    @staticmethod
    def one_based_coords(start, end, chr_size, strand):
        """
        Converts zero based, half-open coords (used in maf) to one-based, fully closed
        """
        if strand>0:
            one_start = start
            one_end = end-1
        else:
            one_start = chr_size - end + 1
            one_end = chr_size - start + 1

        return(one_start, one_end)

def parse_maf(maf_file, store_seqs=True):

    synteny_blocks = []
    block = False
    block_seqs = []
    with open(maf_file) as f:
        for line in f:

            if line.startswith("a"):
                block = True
                continue

            if block:
                if line.startswith("s"):
                    # print(line.split())
                    _, chr_name, start, seq_len, strand_sign, chr_size, seq = line.split()

                    strand = 1 if strand_sign == "=" else -1
                    start = int(start)
                    chr_size = int(chr_size)
                    end = start + int(seq_len)

                    if not store_seqs: seq = None
                    block_seqs.append(MAFseq(chr_name, start, end, strand, chr_size, seq))
                    # print(block_seqs)
                else:
                    if len(block_seqs) > 1:
                        # print(synteny_blocks)
                        synteny_blocks.append(SyntenyBlock(block_seqs, aligned=True))
                    block = False
                    block_seqs = []
        
        if block_seqs:
            if len(block_seqs) > 1:
                synteny_blocks.append(SyntenyBlock(block_seqs, aligned=store_seqs))

    return MAF(synteny_blocks)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("maf",
                        help="path to maf file")
    
    args = parser.parse_args()


if __name__ == "__main__":
    main()