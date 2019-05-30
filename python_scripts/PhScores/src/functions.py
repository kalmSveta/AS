#!/usr/bin/env python2

def Seq_to_bin(seq):
    Dict = {'A':0b0,'T':0b10,'G':0b11,'C':0b1}
    s = 0
    for char in seq:
        try:
            s = s<<2|Dict[char]
        except KeyError:
            s = False
    return(s)

def Index_seq(seq, kmer_length = 5):
    seq_bin = Seq_to_bin(seq)
    if seq_bin is False:
        return(False)
    else:
        seq_indxd = []
        mask = 2**(kmer_length * 2) - 1
        for i in range(len(seq) - (kmer_length - 1)):
            seq_indxd.append(mask & seq_bin)
            seq_bin >>= 2
        return(seq_indxd)