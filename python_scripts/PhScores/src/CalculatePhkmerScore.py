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
sys.path.insert(0, './')
import itertools
from functions import *

def CountOnekmerScore(ph, kmer_length, kmerFreq, id):
    summed_freq = sum([kmerFreq[item] for item in Index_seq(ph.loc[ph.id == id, 'alignment1'].values[0], kmer_length)]) + \
    sum([kmerFreq[item] for item in Index_seq(ph.loc[ph.id == id, 'alignment2'].values[0], kmer_length)])
    length = ph.loc[ph.id == id, 'al1_length'].values[0] + ph.loc[ph.id == id, 'al2_length'].values[0]
    Score = summed_freq / length
    return(Score)


def CalculatePhkmerScore(path_to_kmer_freq, path_to_ph, kmer_length, Threads):
    with open(path_to_kmer_freq, 'rb') as f:
        kmerFreq = pickle.load(f)
    ph = pd.read_csv(path_to_ph, sep='\t')
    ph.drop_duplicates(subset = 'id', inplace = True)
    p = mp.Pool(processes=Threads)
    m = mp.Manager()
    CountOnekmerScore2 = partial(CountOnekmerScore, ph, kmer_length, kmerFreq)
    ids = ph["id"].unique()
    ph['kmerScore'] = p.map(CountOnekmerScore2, ids)
    p.close()
    p.join()
    return(ph)


def FindStackedkmers(structure, sequence, kmer_length):
    structure = structure[:structure.find(')')]
    stacks_in_sequence = []
    while 1:
        dot_position = structure.find('.')
        if dot_position < 0:
            stacks_in_sequence.append(sequence)
            break
        else:
            stacks_in_sequence.append(sequence[:dot_position])
            structure = structure[dot_position + 1:]
            sequence = sequence[dot_position + 1:]
    kmers_in_sequence = [Index_seq(stack, kmer_length) for stack in stacks_in_sequence]
    kmers_in_sequence = list(itertools.chain.from_iterable(kmers_in_sequence))
    return(kmers_in_sequence)


def CountOneEnergyScore(ph, kmer_length, kmerEnergy, mean, sigma, id):
    ph_tmp = ph.loc[ph.id == id]
    kmers_indxd = FindStackedkmers(ph_tmp.structure.values[0], ph_tmp.alignment1.values[0], kmer_length)
    mean_energy = sum([kmerEnergy[item] for item in kmers_indxd]) / len(kmers_indxd)
    z_score = (mean_energy - mean) / (sigma / np.sqrt(len(kmers_indxd)))
    return (z_score)


def CalculatePhEnergyScore(path_to_ph, path_to_kmer_energy, path_to_table, kmer_length, Threads, energy_and_freq = True):
    with open(path_to_kmer_energy, 'rb') as f:
        kmerEnergy = pickle.load(f)
    ph = pd.read_csv(path_to_ph, sep='\t')
    ph.drop_duplicates(subset='id', inplace=True)
    table = pd.read_csv(path_to_table, sep='\t')
    if energy_and_freq:
        mean = sum(table.energy * table.frequency)
        sigma = np.sqrt(sum(((table.energy-mean)**2)*table.frequency))
    else:
        mean = np.mean(table.energy)
        sigma = np.std(table.energy)
    p = mp.Pool(processes=Threads)
    CountOneEnergyScore2 = partial(CountOneEnergyScore, ph, kmer_length, kmerEnergy, mean, sigma)
    ids = ph["id"].unique()
    if energy_and_freq:
        ph['Energy_and_freq_Score'] = p.map(CountOneEnergyScore2, ids)
    else:
        ph['Energy_Score'] = p.map(CountOneEnergyScore2, ids)
    p.close()
    p.join()
    return(ph)


def main(argv):
    path_to_kmer_freq = '../out/kmerFreq.pkl'
    path_to_ph = '../../../conservative_features/whole_human_genome/alignments_whole_human_genome_processed_norm_structures.tsv'
    path_to_kmer_energy = '../out/kmerEnergy.pkl'
    path_to_table = '../out/kmerEnergy.tsv'
    mode = 'kmer'
    kmer_length = 5
    Threads = 10
    try:
        opts, args = getopt.getopt(argv, "h:f:p:k:T:e:t:m:",
                                   ["help=", "frequency=", "ph=", "kmer=", "Threads=", "energy=", 'table=', 'mode='])
    except getopt.GetoptError:
        print('CalculatePhkmerScore.py -f <frequency> -p <ph> -k <kmer> -T <Threads> -m kmer \n ' +
              'OR CalculatePhkmerScore.py -e <energy> -t table -p <ph> -k <kmer> -T <Threads> -m energy \n ' +
              'OR CalculatePhkmerScore.py -e <energy> -t table -p <ph> -k <kmer> -T <Threads> -m energy_and_frequency')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('CalculatePhkmerScore.py -f <frequency> -p <ph> -k <kmer> -T <Threads> -m kmer \n ' +
                  'OR CalculatePhkmerScore.py -e <energy> -t table -p <ph> -k <kmer> -T <Threads> -m energy \n ' +
                  'OR CalculatePhkmerScore.py -e <energy> -t table -p <ph> -k <kmer> -T <Threads> -m energy_and_frequency')
            sys.exit()
        elif opt in ("-f", "--frequency"):
            path_to_kmer_freq = arg
        elif opt in ("-p", "--ph"):
            path_to_ph = arg
        elif opt in ("-k", "--kmer"):
            kmer_length = int(arg)
        elif opt in ("-T", "--Threads"):
            Threads = int(arg)
        elif opt in ("-e", "--energy"):
            path_to_kmer_energy = arg
        elif opt in ("-t", "--table"):
            path_to_table = arg
        elif opt in ("-m", "--mode"):
            mode = arg
    if mode == 'kmer':
        ph = CalculatePhkmerScore(path_to_kmer_freq, path_to_ph, kmer_length, Threads)
    else:
        if mode == 'energy':
            ph = CalculatePhEnergyScore(path_to_ph, path_to_kmer_energy, path_to_table, kmer_length, Threads,
                                        energy_and_freq=False)
        else:
            ph = CalculatePhEnergyScore(path_to_ph, path_to_kmer_energy, path_to_table, kmer_length, Threads,
                                   energy_and_freq=True)
    ph.to_csv('../out/ph_w_score_' + mode + '.tsv', sep = '\t', header=True, quoting=False, index = False)
    print(ph.head())


if __name__ == '__main__':
    main(sys.argv[1:])
