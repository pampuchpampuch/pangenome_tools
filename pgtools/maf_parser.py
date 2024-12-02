import argparse
import sys
from itertools import combinations, product
from pgtools.pangenome import BaseSeq, SeqCollection, Pangenome

from pgtools.utils import reverse_coords

class MAF(Pangenome):
    """
    represents a synteny file as a set of synteny blocks
    """
    def __init__(self, synteny_blocks):
        super().__init__(synteny_blocks)
        # self.synteny_blocks = synteny_blocks
        chr_names = []
        for syn_block in synteny_blocks:
            chr_names += [maf_seq.seq_name for maf_seq in syn_block.sequences]
        self.chr_names = list(set(chr_names))
    
    def add_synteny_block(self, syn_block):
        self.seq_collections[syn_block.id] = syn_block
        # self.synteny_blocks.append(syn_block)

    # def write_to_file(self, out_file, chr_names=None):
    #     """
    #     Writes MAF object into a maf file. Filtering on sequence names
    #     is possible
    #     """
    #     out_file = open(out_file,"w")
    #     out_file.write("##maf version=1 scoring=N/A\n\n")

    #     if not chr_names:
    #         for block in self.synteny_blocks:
    #             out_file.write(block.MAF_repr())
        
    #     else:
    #          for block in self.synteny_blocks:
    #             out_file.write(block.MAF_repr(chr_names))           

    #     out_file.close()

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
            
    def get_filtered_blocks_by_strand(self, contig_names, symmetrical_invert = False):
        """
        filteres blocks by contig names and categorizes them by strandness
        """
        maf_blocks = {"".join(strand_comb): [] for strand_comb in product(["+","-"], repeat=len(contig_names))}
        # maf_blocks = {"++": [], "+-": [], "--": [], "-+": []}
        for block in self.synteny_blocks:
            if not set(contig_names).issubset({seq.seq_name for seq in block.block_seqs}):
                continue
            dict_key = ""
            coords = []
            seqs_dict = {seq.seq_name: seq for seq in block.block_seqs}
            for genome in contig_names:
                seq = seqs_dict[genome]
                coords.append((seq.start, seq.end))
                if seq.strand > 0:
                    dict_key += "+"
                else:
                    dict_key += "-"       

            maf_blocks[dict_key].append(tuple(coords))
        
        if symmetrical_invert:
            assert len(contig_names) == 2, "Symmetrcial invert (considering +- and -+\
            simultaniously is possible only for filtering with two sequences)"
            for strandness, coords in maf_blocks.items():
                opposite_strandness = ""
                for el in strandness:
                    if el == "+":
                        opposite_strandness += "-"
                    else:
                        opposite_strandness += "+"

                rev_coords = [(reverse_coords(coords_[0])[:2], reverse_coords(coords_[1])[:2]) for coords_ in coords]
                for coords_1, coords_2 in coords:
                    rev_coords_1 = reverse_coords(coords_1[0], coords_1[1], 1)[:2]


        return maf_blocks
    
    def mismatch_pct(self):
        """
        Calculates identity for each pair of sequence in each block.
        Does not distinguish between them
        """
        n_gaps = 0
        n_cols = 0

        for block in self.synteny_blocks:
            n_gaps_, n_cols_ = block.count_mismatch_cols()
            n_gaps += n_gaps_
            n_cols += n_cols_
        
        return round((n_cols - n_gaps) * 100 / n_cols, 2)
    
    def filter(self, chr_names):
        chr_names = set(chr_names)
        # blocks_dict = {block.id: set(seq.seq_name for seq in block.block_seqs)
        #                for block in self.synteny_blocks}
        blocks_dict = {block.id: block for block in self.synteny_blocks}
        # for id, genome_names in blocks_dict.items():
        filtered_blocks = [blocks_dict[id].filter(chr_names) for id, block in blocks_dict.items()
                              if chr_names.issubset({seq.seq_name for seq in block.block_seqs})]
        return(MAF(filtered_blocks))

class SyntenyBlock(SeqCollection):

    """
    represents synteny block as a set of sequences with an id
    from coords file
    """

    def __init__(self, block_seqs, aligned=False, id=None):
        super().__init__(id, block_seqs)
        # consider adding id?
        self.id = id
        # self.block_seqs = block_seqs
        self.aligned = aligned
    
    def get_block_seqs(self):
        return list(self.seq_dict.values())

    def add_block_seq(self,block_seq):
        # self.block_seqs.append(block_seq)
        self.seq_dict[block_seq.seq_name] = block_seq
    
    # def get_contig_names(self):
    #     return set(self.seq_dict.keys())
        # return {seq.seq_name for seq in self.block_seqs}

    #         maf_str = "a\n"

    #         if chr_names:
    #             for maf_seq in self.seq_dict.values():
    #                 if maf_seq.seq_name in chr_names:
    #                     maf_str += maf_seq.MAF_repr() 

    #         else:
    #             for maf_seq in self.seq_dict.values():
    #                 maf_str += maf_seq.MAF_repr() 

    #         maf_str += '\n'       

    #     else:
    #         sys.exit("MAF block alignment to be implemented")
        
    #     return maf_str
    
    def filter(self, seq_names):
        seqs = [seq for seq in self.seq_dict.values() if seq.seq_name in seq_names]
        return SyntenyBlock(seqs, id = self.id)
    
    def count_mismatch_cols(self):
        assert self.aligned
        seq_idx = [i for i in range(len(self.seq_dict.values()))]
        seqs = [maf_seq.seq for maf_seq in self.seq_dict.values()]

        n_cols = 0
        n_gaps = 0
        for idx_1, idx_2 in combinations(seq_idx, 2):
            n_gaps_, n_cols_ = MAFseq.count_mismatch_cols(seqs[idx_1], seqs[idx_2])
            n_cols += n_cols_
            n_gaps += n_gaps_
        
        return n_gaps, n_cols
    
    def mean_gap_content(self):
        all_seqs = list(self.sequences)
        n_seqs = len(all_seqs)
        gap_cont = []
        seq_len = len(all_seqs[0].seq)
        for i in range(seq_len):
            i_pos_char = [seq_.seq[i] for seq_ in all_seqs]
            gap_cont.append(i_pos_char.count("-") / n_seqs)
        
        return sum(gap_cont) / len(gap_cont)

class MAFseq(BaseSeq):
    """
    Represents sequence in maf block
    """
    def __init__(self, seq_name: str, start: int, end: int, strand: int, src_size: int, in_format="maf", seq: str = None, coord_system="maf"):
        super().__init__(seq_name, start, end, strand, src_size, in_format, seq, coord_system)
    # def __init__(self, chr_name, start, end, strand, chr_size, seq):
    #     super().__init__(chr_name, start, end, strand)
    #     # self.seq_name = chr_name
    #     # self.start = start
    #     # self.end = end
    #     # #self.seq_len = int(end) - int(start)
    #     # self.strand = strand
    #     self.chr_size = chr_size
    #     self.seq = seq


    def __len__(self):
        return (self.end - self.start + 1)
    
        
    #     strand_sign = "+" if self.strand > 0 else "-"

    #     s_line=f"s\t{self.seq_name}\t{self.start}\t{len(self)}\t{strand_sign}\t{self.chr_size}\t{self.seq}\n"

    #     return(s_line)
            
    
    @staticmethod
    def count_mismatch_cols(seq1, seq2):
        """
        Count number of columns where at list one sequence
        has a gap character. Provided sequences have to be aligned.
        """
        assert len(seq1) == len(seq2)

        n_gaps = 0
        n_cols = 0
        for i in range(len(seq1)):
            if (seq1[i] == "-" or seq2[i] == "-"):
                continue
            elif seq1[i] != seq2[i]:
                n_gaps +=1
                
            n_cols += 1
            
        return n_gaps, n_cols

    @staticmethod
    def gff_coords(start, end, chr_size, strand):
        """
        For now, converts to gff compatible coords
        """
        #TO DO: Converts zero based, half-open coords (used in maf) to one-based, fully closed
        # print(strand)
        if strand>0:
            one_start = start
            one_end = end
        else:
            one_start = chr_size - end
            one_end = chr_size - start

        return(one_start, one_end)
    


def parse_maf(maf_file, store_seqs=True):
    """
    Parses maf file and return MAF object.
    Possible to not store sequences (as they
    are not needed for all applications to conserve memory)
    """
    ### TO DO ##
    # possiblity to input compressed (.gz) maf files
    ############
    synteny_blocks = []
    block = False
    block_seqs = []
    id = 0
    with open(maf_file) as f:
        for line in f:

            if line.startswith("a"):
                block = True
                continue

            if block:
                if line.startswith("s"):
                    # print(line.split())
                    _, chr_name, start, seq_len, strand_sign, chr_size, seq = line.split()
                    
                    # print(strand_sign)
                    strand = 1 if strand_sign == "+" else -1
                    start = int(start)
                    chr_size = int(chr_size)
                    end = start + int(seq_len)

                    if not store_seqs: seq = None
                    block_seqs.append(MAFseq(chr_name, start, end, strand, chr_size, in_format="maf", seq=seq))
                    # print(block_seqs)
                else:
                    if len(block_seqs) > 1:
                        # print(synteny_blocks)
                        synteny_blocks.append(SyntenyBlock(block_seqs, aligned=True, id=id))
                        id += 1
                    block = False
                    block_seqs = []
        
        if block_seqs:
            if len(block_seqs) > 1:
                synteny_blocks.append(SyntenyBlock(block_seqs, aligned=store_seqs, id=id))
                id += 1
    return MAF(synteny_blocks)
