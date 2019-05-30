#!/usr/bin/env python2

import pandas as pd
import numpy as np
import glob, re, sys, getopt
from functools import partial
import sys
sys.path.append('../../../tools/')
from pyfaidx import Fasta
from subprocess import call
import binascii
import multiprocessing as mp
import pickle
import itertools
from Bio.Seq import Seq
sys.path.insert(0, './')
from functions import *

stacking_matrix = np.load("../lib/stacking_matrix.npy")

def CalculateStackingEnergy(seq, seq_compl):
    Dic_bp = {'CG': 0, 'GC': 1, 'GT': 2, 'TG': 3, 'AT': 4, 'TA': 5}
    energy = 0
    i = 0
    j = 0
    while i < len(seq) - 1 and j < len(seq) - 1:
        energy_add = stacking_matrix[Dic_bp.get(seq_compl[i + 1] + seq[j + 1],6)][Dic_bp.get(seq[j] + seq_compl[i],6)]
        i += 1
        j += 1
        energy += energy_add
    return(energy / 100)


def CalculatekmerFrequency(path_to_kmer_freq_df, kmer_length = 5, Threads = 1):
    bases = ['A', 'T', 'G', 'C']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=kmer_length)]
    kmerEnergy = {kmer: CalculateStackingEnergy(kmer, str(Seq(kmer).complement())) for kmer in kmers}
    kmerEnergy_df = pd.DataFrame.from_dict(kmerEnergy, orient='index')
    kmerEnergy_df['kmer'] = kmerEnergy_df.index
    kmerEnergy_df.rename({0: 'energy'}, inplace=True, axis=1)
    Index_seq2 = partial(Index_seq, kmer_length = kmer_length)
    kmerEnergy_df['kmer_indxd'] = map(Index_seq2, kmerEnergy_df.kmer)
    kmerEnergy_df['kmer_indxd'] = [item[0] for item in kmerEnergy_df.kmer_indxd]
    kmerEnergy = kmerEnergy_df[['kmer_indxd', 'energy']].set_index('kmer_indxd')['energy'].to_dict()
    kmer_freq_df = pd.read_csv(path_to_kmer_freq_df, sep='\t')
    merged = pd.merge(kmerEnergy_df, kmer_freq_df, on=['kmer', 'kmer_indxd'])
    merged.to_csv('../out/kmerEnergy.tsv', sep='\t', header=True, quoting=False)
    with open('../out/kmerEnergy.pkl', 'wb') as f:
        pickle.dump(kmerEnergy, f, pickle.HIGHEST_PROTOCOL)
    return (kmerEnergy_df)



def main(argv):
    path_to_kmer_freq_df = '../out/kmerFreq.tsv'
    kmer_length = 5
    Threads = 10
    try:
        opts, args = getopt.getopt(argv, "h:k:T:",
                                   ["help=", "kmer=", "Threads="])
    except getopt.GetoptError:
        print('CalculatekmerEnergy.py -k <kmer> -T <Threads>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('CalculatekmerEnergy.py -k <kmer> -T <Threads>')
            sys.exit()
        elif opt in ("-k", "--kmer"):
            kmer = int(arg)
        elif opt in ("-T", "--Threads"):
            Threads = int(arg)
    CalculatekmerFrequency(path_to_kmer_freq_df, kmer_length = 5, Threads = Threads)

if __name__ == '__main__':
    main(sys.argv[1:])