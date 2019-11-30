#!/usr/bin/env python2
from numpy import argmin, unravel_index, full, empty, load, set_printoptions, argwhere, argsort, random, digitize
from math import ceil
import numpy as np
import re, sys, getopt, itertools, binascii, time
from functools import partial
sys.path.append('../../../tools/')
from pyfaidx import Fasta
import pandas as pd
import multiprocessing as mp
sys.path.insert(0, './')
from fold_new_no_intersection import * 
from Bio.Seq import Seq
from collections import Counter


inf = float('inf')


def GetSequencesForDF(genome, row):
    return (str(genome[row['chr']][row['start_interval'] - 1:row['end_interval']]))

def MakeComplement(row):
    if row['strand'] == '-':
        return(str(Seq(row['sequences']).complement()))
    else:
        return(row['sequences'])


def Find_panhandles_one_gene(lock, df, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                             need_suboptimal, kmers_stacking_matrix, seed, gene):
    with open('../out/genes_done2.txt', 'a') as done_f:
        with lock:
            done_f.write(str(gene) + '\n')
    print('Working with gene ' + gene)
    results_one_gene = []
    results_one_gene_table = pd.DataFrame(
        {'gene': [], 'energy': [],
         'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
         'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
    interval_indexes = df.index[df['gene_chr_start_end_strand'] == gene].tolist()
    counts_close = 0
    for left_interval_index in range(interval_indexes[0], interval_indexes[len(interval_indexes) - 1]):
        seq1 = df.loc[left_interval_index, 'sequences'][:]
        seq1_indxd = df.loc[left_interval_index, 'sequences_indxd'][:]
        for right_interval_index in range(left_interval_index + 1, interval_indexes[len(interval_indexes) - 1] + 1):
            seq2 = df.loc[right_interval_index, 'sequences'][:]
            seq2_indxd = df.loc[right_interval_index, 'sequences_indxd'][:]
            if abs(df.loc[right_interval_index, 'start_interval'] - df.loc[left_interval_index, 'end_interval']) \
                    <= panhandle_length_threshold:
                counts_close += 1
                align = FindMinEnLocAlkmer(seq1, seq2, seq1_indxd, seq2_indxd, k, energy_threshold,
                                                     handle_length_threshold, need_suboptimal, kmers_stacking_matrix)
                if align != 0:
                    results_one_gene.append([align, df.loc[left_interval_index, 'interval_chr_start_end_strand'],
                                             df.loc[right_interval_index, 'interval_chr_start_end_strand']])

    for result in results_one_gene:
        for alignment in result[0]:
            results_one_gene_table = results_one_gene_table.append(
                {'gene': gene, 'energy': alignment[0], 'interval1': result[1], 'interval2': result[2],
                 'start_al1': alignment[1], 'end_al1': alignment[2],
                 'start_al2': alignment[3], 'end_al2': alignment[4],
                 'alignment1': alignment[5], 'alignment2': alignment[6], 'structure': alignment[7]}, ignore_index=True)
    with open('../out/random_panhandles' + str(seed) + '.tsv', 'a') as f:
        with lock:
            results_one_gene_table.to_csv(f, sep='\t', index=False, header=False)

    with open('../out/counts_close_' + str(seed) + '.txt', 'a') as f:
        with lock:
            f.write(str(counts_close) + '\n')


def Shuffle(df, what, seed):
    print('Start shuffling..')
    if what == 'genes':
        data = (df.end_interval + df.start_interval)//2
        bins = range(min(data), max(data)+1, 10000)[1:]
        df['bins'] = digitize(data, bins, right=True)
        df.gene_chr_start_end_strand = df.groupby(["bins"])['gene_chr_start_end_strand'].transform(random.RandomState(seed=seed).permutation)
        df['chr'], df['start_gene'], df['end_gene'], df['strand'] = df['gene_chr_start_end_strand'].str.split('_', 4).str
        df.sort_values(by=['gene_chr_start_end_strand', 'start_interval', 'end_interval'], inplace=True)
        df.reset_index(inplace = True)
    elif what == 'sequences':
        print('shuffle on sequences')
        df['length'] = (df.end_interval - df.start_interval)
        df.sort_values(by=['length'], inplace=True)
        df.index = range(1,len(df.length)+1)
        df['bins'] = pd.qcut(df.index, 50, labels=False)
        df['sequences'] = df.groupby(["bins"])['sequences'].transform(random.RandomState(seed=seed).permutation)
        #df['sequences'] = df['sequences'].transform(random.RandomState(seed=seed).permutation)
        df.drop(['bins', 'length'], axis=1, inplace = True)
        df.sort_values(by=['gene_chr_start_end_strand', 'start_interval', 'end_interval'], inplace=True)
        df.reset_index(inplace = True)
    elif what == 'sequences_GC':
        print('shuffle on sequences and GC content')
        df['length'] = (df.end_interval - df.start_interval)
        df['GC_content'] = df.apply(lambda row: (Counter(row['sequences'])['G'] + Counter(row['sequences'])['C']) * 1.0 / row['length'], axis = 1)
        df.sort_values(by=['length', 'GC_content'], inplace=True)
        df.index = range(1,len(df.length)+1)
        df['bins'] = pd.qcut(df.index, 50, labels=False)
        df['sequences'] = df.groupby(["bins"])['sequences'].transform(random.RandomState(seed=seed).permutation)
        #df['sequences'] = df['sequences'].transform(random.RandomState(seed=seed).permutation)
        df.drop(['bins', 'length', 'GC_content'], axis=1, inplace = True)
        df.sort_values(by=['gene_chr_start_end_strand', 'start_interval', 'end_interval'], inplace=True)
        df.reset_index(inplace = True)

    print('Shuffled!')
    return(df)

def Find_Random_panhandles(path_to_intervals, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix, N_seeds, strandness, what):
    start_time = time.time()
    df = pd.read_csv(path_to_intervals, sep='\t')
    df["gene_chr_start_end_strand"] = df.chr + "_" + df.start_gene.map(str) + "_" + df.end_gene.map(str) + "_" + df.strand
    df["interval_chr_start_end_strand"] = df["chr"] + "_" + df["start_interval"].map(str) + "_" + df["end_interval"].map(str) + "_" + df["strand"]
    df['start_interval'] = df['start_interval'].astype(int)
    df['end_interval'] = df['end_interval'].astype(int)
    if not ('sequences' in list(df.columns.values)):
        print('Attaching sequences..')
        genome = Fasta(genome_file)
        GetSequencesForDF2 = partial(GetSequencesForDF, genome)
        df.loc[:, 'sequences'] = df.apply(GetSequencesForDF2, axis=1)
        if strandness:
            print("Making complement of minus strand..")
            df.loc[:, 'sequences'] = df.apply(MakeComplement, axis=1)         
        df.to_csv("../out/intervals_with_seqs.tsv", sep='\t', index=False)
    df.sequences = map(lambda x: x.upper(), df['sequences'])

    for seed in range(1, N_seeds + 1):  
        print(seed)  
        ##shuffle
        df_new = df.copy() 
        df_new = Shuffle(df_new, what, seed)
        df_new["sequences_indxd"] = df_new['sequences'].apply(lambda x: Index_seq(x, k))
        df_new = df_new.loc[df_new.sequences_indxd != False]
        print("Creating files..")
        with open('../out/genes_done2.txt', 'w') as done_f:
            done_f.write('Started alignment: \n')
        results_one_gene_table = pd.DataFrame(
            {'gene': [], 'energy': [],
             'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
             'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
        with open('../out/random_panhandles' + str(seed) + '.tsv', 'w') as f:
            results_one_gene_table.to_csv(f, sep='\t', index=False, header=True)
        with open('../out/counts_close_' + str(seed) + '.txt', 'w') as f:
            f.write('')
        print('Start to align..')
        p = mp.Pool(processes=threads)
        print('Created pool')
        m = mp.Manager()
        print('Created manager')
        lock = m.Lock()
        print('Created lock')
        f_shuffled = '../out/intervals_shuffled_' + str(seed) + '.tsv'
        df_new.to_csv(f_shuffled, sep='\t', index=False, header=True)
        Find_panhandles_one_gene2 = partial(Find_panhandles_one_gene, lock, df_new, energy_threshold, handle_length_threshold,
                                            panhandle_length_threshold, k, need_suboptimal, kmers_stacking_matrix, seed)
        print('Created partial')
        genes = df_new["gene_chr_start_end_strand"].unique()
        print('Created genes')
        p.map(Find_panhandles_one_gene2, genes)
        p.close()
        p.join()
    print("all done!")
    print(time.time() - start_time)
    return (0)


def main(argv):
    path_to_intervals = '../out/folding/intervals_with_seqs.tsv'
    genome_file = './../../../genomes/GRCh37.p13.genome.fa'
    k = 5
    panhandle_length_threshold = 10000
    handle_length_threshold = 10
    threads = 5
    energy_threshold = -15
    need_suboptimal = True
    GT_threshold = 2
    N_seeds = 10
    strandness = True
    what = 'sequences_GC'
    try:
        opts, args = getopt.getopt(argv, "h:i:g:k:p:a:t:e:u:d:N:s:w:",
                                   ["help=", "intervals=", "genome=", "k=", "panh_len_max", "handle_len_min", "threads",
                                    "energy_max", "need_subopt", "gt_threshold", "strand", "what"])
    except getopt.GetoptError:
        print(
            'AssessFDR.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max> ' +
            '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -N <N_seeds> -s <strand> -w <what>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'AssessFDR.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max> ' +
                '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -N <N_seeds> -s <strand> -w <what>')
            sys.exit()
        elif opt in ("-i", "--intervals"):
            path_to_intervals = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-k", "--kmer"):
            k = int(arg)
        elif opt in ("-p", "--panh_len_max"):
            panhandle_length_threshold = int(arg)
        elif opt in ("-a", "--handle_len_min"):
            handle_length_threshold = int(arg)
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-e", "--energy_max"):
            energy_threshold = float(arg)
        elif opt in ("-u", "--need_subopt"):
            need_suboptimal = bool(arg)
        elif opt in ("-d", "--gt_threshold"):
            GT_threshold = int(arg)
        elif opt in ("-N", "-N_seeds"):
            N_seeds = int(arg)
        elif opt in ("-s", "--strand"):
            strandness = bool(arg)
        elif opt in ("-w", "--what"):
            what = arg

    kmers_stacking_matrix = load("../lib/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy")
    Find_Random_panhandles(path_to_intervals, energy_threshold, handle_length_threshold,
                    panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix, N_seeds, strandness, what)


if __name__ == '__main__':
    main(sys.argv[1:])

