#!/usr/bin/env python2
from numpy import argmin, unravel_index, full, empty, load, set_printoptions, argwhere, argsort
from math import ceil
import re, sys, getopt, itertools, binascii, time
from functools import partial
sys.path.append('../../../tools/')
from pyfaidx import Fasta
import pandas as pd
import multiprocessing as mp
sys.path.insert(0, './')
import fold_new_no_intersection

inf = float('inf')


def GetSequencesForDF(genome, row):
    return (str(genome[row['chr']][row['start_interval'] - 1:row['end_interval']]))


def Find_panhandles_one_gene(lock, df, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                             need_suboptimal, kmers_stacking_matrix, gene):
    with open('../out/genes_done.txt', 'a') as done_f:
        with lock:
            done_f.write(str(gene) + '\n')
    print('Working with gene ' + gene)
    results_one_gene = []
    results_one_gene_table = pd.DataFrame(
        {'gene': [], 'energy': [],
         'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
         'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
    interval_indexes = df.index[df['gene_chr_start_end_strand'] == gene].tolist()
    for left_interval_index in range(interval_indexes[0], interval_indexes[len(interval_indexes) - 1]):
        seq1 = df.loc[left_interval_index, 'sequences'][:]
        seq1_indxd = df.loc[left_interval_index, 'sequences_indxd'][:]
        #for right_interval_index in range(left_interval_index, interval_indexes[len(interval_indexes) - 1] + 1):
        for right_interval_index in range(left_interval_index + 1, interval_indexes[len(interval_indexes) - 1] + 1):
            seq2 = df.loc[right_interval_index, 'sequences'][:]
            seq2_indxd = df.loc[right_interval_index, 'sequences_indxd'][:]
            if abs(df.loc[right_interval_index, 'start_interval'] - df.loc[left_interval_index, 'end_interval']) \
                    <= panhandle_length_threshold:
                align = fold_new_no_intersection.FindMinEnLocAlkmer(seq1, seq2, seq1_indxd, seq2_indxd, k, energy_threshold,
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
    with open('../out/panhandles.tsv', 'a') as f:
        with lock:
            results_one_gene_table.to_csv(f, sep='\t', index=False, header=False)



def Find_panhandles(path_to_intervals, energy_threshold, handle_length_threshold, panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix):
    start_time = time.time()
    df = pd.read_csv(path_to_intervals, sep='\t')
    df["gene_chr_start_end_strand"] = df.chr + "_" + df.start_gene.map(str) + "_" + df.end_gene.map(str) + "_" + df.strand
    df["interval_chr_start_end_strand"] = df["chr"] + "_" + df["start_interval"].map(str) + "_" + df[
        "end_interval"].map(str) + "_" + df["strand"]
    df['start_interval'] = df['start_interval'].astype(int)
    df['end_interval'] = df['end_interval'].astype(int)
    df.sort_values(by=['gene_chr_start_end_strand', 'start_interval', 'end_interval'], inplace=True)
    df.reset_index(inplace = True)
    if not ('sequences' in list(df.columns.values)):
        print('Attaching sequences..')
        genome = Fasta(genome_file)
        GetSequencesForDF2 = partial(GetSequencesForDF, genome)
        df.loc[:, 'sequences'] = df.apply(GetSequencesForDF2, axis=1)
        df.to_csv("../out/intervals_with_seqs.tsv", sep='\t', index=False)
    df.sequences = map(lambda x: x.upper(), df['sequences'])
    df["sequences_indxd"] = df['sequences'].apply(lambda x: fold_new_no_intersection.Index_seq(x, k))
    df = df.loc[df.sequences_indxd != False]
    
    print("Creating files..")
    with open('../out/genes_done.txt', 'w') as done_f:
        done_f.write('Started alignment: \n')
    results_one_gene_table = pd.DataFrame(
        {'gene': [], 'energy': [],
         'start_al1': [], 'end_al1': [], 'start_al2': [], 'end_al2': [],
         'alignment1': [], 'alignment2': [], 'structure': [], 'interval1': [], 'interval2': []})
    with open('../out/panhandles.tsv', 'w') as f:
        results_one_gene_table.to_csv(f, sep='\t', index=False, header=True)
    print('Start to align..')
    p = mp.Pool(processes=threads)
    m = mp.Manager()
    lock = m.Lock()
    Find_panhandles_one_gene2 = partial(Find_panhandles_one_gene, lock, df, energy_threshold, handle_length_threshold,
                                        panhandle_length_threshold, k, need_suboptimal, kmers_stacking_matrix)
    genes = df["gene_chr_start_end_strand"].unique()
    p.map(Find_panhandles_one_gene2, genes)
    p.close()
    p.join()
    print("all done!")
    print(time.time() - start_time)
    return (0)


def main(argv):
    path_to_intervals = '../out/intervals_with_seqs.tsv'
    genome_file = './../../../genomes/GRCh37.p13.genome.fa'
    k = 5
    panhandle_length_threshold = 10000
    handle_length_threshold = 10
    threads = 5
    energy_threshold = -15
    need_suboptimal = True
    GT_threshold = 2
    try:
        opts, args = getopt.getopt(argv, "h:i:g:k:p:a:t:e:u:d:",
                                   ["help=", "intervals=", "genome=", "k=", "panh_len_max", "handle_len_min", "threads",
                                    "energy_max", "need_subopt", "gt_threshold"])
    except getopt.GetoptError:
        print(
            'FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max> ' +
            '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max> ' +
                '-a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold>')
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
    kmers_stacking_matrix = load("../lib/" + str(k) + str(GT_threshold) + "mers_stacking_energy_binary.npy")
    Find_panhandles(path_to_intervals, energy_threshold, handle_length_threshold,
                    panhandle_length_threshold, k,
                    genome_file, threads, need_suboptimal, kmers_stacking_matrix)


if __name__ == '__main__':
    main(sys.argv[1:])
