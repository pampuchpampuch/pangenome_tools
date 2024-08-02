def reverse_coords(start: int, end: int, chr_len: int, strand: int):
    """
    ### !prob should a method of sequence class all other sequence class inherit after
    """
    rev_start = chr_len - end
    rev_end = chr_len - start + 1
    return (rev_start, rev_end, -1 * strand)

