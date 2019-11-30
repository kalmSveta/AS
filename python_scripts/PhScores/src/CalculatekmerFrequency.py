#!/usr/bin/env python2

import pandas as pd
import numpy as np
import glob, re, sys, getopt
from functools import partial
import sys
sys.path.append('../../../tools/')
from pyfaidx import Fasta
from subprocess import call
import pickle
sys.path.insert(0, './')
from functions import *

def GetSequencesForDF(genome, row):
    return (str(genome[row['chr']][row['start']-1:row['end']]))

def readFasta(path_to_fasta):
    kmerFreqDict = {}
    with open(path_to_fasta) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                continue
            sequence = line
            if sequence not in kmerFreqDict:
                kmerFreqDict[sequence] = []
            kmerFreqDict[sequence].append(active_sequence_name)
    return(kmerFreqDict)

def CalculatekmerFrequency(path_to_population, genome_file, kmer_length = 5, Threads = 1):
    call('bedtools sort -i  ' + path_to_population + ' | bedtools merge -i stdin > ../out/population.bed', shell=True)
    df = pd.read_csv('../out/population.bed', sep = '\t')
    genome = Fasta(genome_file)
    GetSequencesForDF2 = partial(GetSequencesForDF, genome)
    df.columns = ['chr', 'start', 'end']
    seqs = df.apply(GetSequencesForDF2, axis=1)
    out_file = open("../out/population.fa", "w")
    for i in range(len(seqs)):
        out_file.write(">" + str(i) + "\n" + seqs[i] + "\n")
    out_file.close()

    call('./jellyfish-linux count -m ' + str(kmer_length) + ' -s 100M -t ' + str(Threads) + ' ../out/population.fa', shell=True)
    call('./jellyfish-linux dump mer_counts.jf > ../out/kmer_counts.fa', shell=True)
    call('mv mer_counts.jf ../out/kmer_counts.jf', shell=True)

    kmerFreqDict = readFasta('../out/kmer_counts.fa')
    kmerFreq_df = pd.DataFrame.from_dict(kmerFreqDict, orient='index')
    kmerFreq_df['kmer'] = kmerFreq_df.index
    kmerFreq_df.rename({0: 'frequency'}, inplace=True, axis=1)
    Index_seq2 = partial(Index_seq, kmer_length = kmer_length)
    kmerFreq_df['kmer_indxd'] = map(Index_seq2, kmerFreq_df.kmer)
    kmerFreq_df['kmer_indxd'] = [item[0] for item in kmerFreq_df.kmer_indxd]
    kmerFreq_df.frequency = kmerFreq_df.frequency.values.astype('int')
    s = sum(kmerFreq_df.frequency.values)
    kmerFreq_df.frequency = kmerFreq_df.frequency / s
    kmerFreq = kmerFreq_df[['kmer_indxd', 'frequency']].set_index('kmer_indxd')['frequency'].to_dict()
    kmerFreq_df.to_csv('../out/kmerFreq.tsv', sep='\t', header=True, quoting=False)
    with open('../out/kmerFreq.pkl', 'wb') as f:
        pickle.dump(kmerFreq, f, pickle.HIGHEST_PROTOCOL)
    return(0)


def main(argv):
    path_to_population = '../../../conservative_features/whole_human_genome/preStructure/conin_grange_prot_coding.bed'
    genome_file = './../../../genomes/GRCh37.p13.genome.fa'
    kmer_length = 5
    Threads = 10

    try:
        opts, args = getopt.getopt(argv, "h:p:g:k:T:",
                                   ["help=", "population=", "genome=", "kmer=", "Threads="])
    except getopt.GetoptError:
        print('CalculatekmerFrequency.py -p <population> -g <genome> -k <kmer> -T <Threads>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('CalculatekmerFrequency.py -p <population> -g <genome> -k <kmer> -T <Threads>')
            sys.exit()
        elif opt in ("-p", "--population"):
            path_to_population = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-k", "--kmer"):
            kmer_length = int(arg)
        elif opt in ("-T", "--Threads"):
            Threads = int(arg)
    CalculatekmerFrequency(path_to_population, genome_file, kmer_length, Threads)


if __name__ == '__main__':
    main(sys.argv[1:])
