#!/usr/bin/env python2
import time
from functools import partial
import sys, getopt
sys.path.insert(0, './')
import RecalculateEnergies
from re import sub
from random import choice
import multiprocessing as mp
from functools import partial
import pandas as pd
import numpy as np
from subprocess import call


def Intersect_with_mutations(path_to_ph, path_to_mut, title):
    ph = pd.read_csv(path_to_ph, sep="\t")
    ph_bed = ph[["chr", "panhandle_start", "panhandle_left_hand", "id", "energy", "strand", "structure", "alignment1",
                 "alignment2", "al1_length", "al2_length", "panhandle_right_hand", "panhandle_end"]]
    ph_bed.panhandle_start = ph.panhandle_start.astype(int)
    ph_bed.panhandle_left_hand = ph.panhandle_left_hand.astype(int)
    ph_bed.panhandle_right_hand = ph.panhandle_right_hand.astype(int)
    ph_bed.panhandle_end = ph.panhandle_end.astype(int)
    ph_bed.to_csv("../out/panhandles.bed13", index=False, header=False, sep="\t")
    call(['./Intersect_with_mutations.sh', "../out/panhandles.bed13", path_to_mut, title])
    ph_left_handles = pd.read_csv('../out/panhandles_with_mutations_left_handles.bed', sep="\t", header=None)
    ph_right_handles = pd.read_csv('../out/panhandles_with_mutations_right_handles.bed', sep="\t", header=None)
    ph_left_handles.drop([13, 14], axis=1, inplace=True)
    ph_left_handles.columns = ["chr", "panhandle_start", "panhandle_left_hand", "panhandle_id", "energy", "strand",
                               "structure", "alignment1", "alignment2", "al1_length", "al2_length", "panhandle_right_hand",
                               "panhandle_end", "mut_coord", "mutation_info"]
    ph_left_handles['panhandle_hand'] = "left"
    ph_right_handles.drop([13, 14], axis=1, inplace=True)
    ph_right_handles.columns = ["chr", "panhandle_right_hand", "panhandle_end", "panhandle_id", "energy", "strand",
                               "structure", "alignment1", "alignment2", "al1_length", "al2_length", "panhandle_start",
                                "panhandle_left_hand", "mut_coord", "mutation_info"]
    ph_right_handles['panhandle_hand'] = "right"
    ph = ph_left_handles.append(ph_right_handles, ignore_index=True)
    ph[['mut_from', "mut_to"]] = ph.mutation_info.str.split("_", 4, expand=True)[[2, 4]]
    ph.to_csv("../out/panhandles_with_mutations.tsv", sep="\t", index=False)


def MakeRelativeCoords(row):
    if row['strand'] == '+': 
        if row['panhandle_hand'] == 'left':
            return (row['mut_coord'] - row['panhandle_start'])
        else:
            return (row['mut_coord'] - row['panhandle_right_hand'])
    else:
        if row['panhandle_hand'] == 'left':
            return (row['panhandle_left_hand'] - row['mut_coord'])
        else:
            return (row['panhandle_end'] - row['mut_coord'])


def IndexesOfSymbol(string, symbol):
    return [i for i, ltr in enumerate(string) if ltr == symbol]

def RandomizeMutations(seq, seq_compl, mut_to, mut_from, seq_long, seq_compl_long, mut_coord_relative, mut_seq, part1_long, part2_long, how = "nt_type"):
    print('Randomizing...')
    if how == "nt_type":
        nts = ['A', 'T', 'G', 'C']
        mut_to_random = choice(list(set(nts) - set([mut_to, mut_from])))
        print('Found random mut_to')
        seq_long_random, seq_compl_long_random, part1_long_random, part2_long_random = RecalculateEnergies.Mutate(
            seq_long, seq_compl_long, mut_coord_relative, mut_seq, mut_to_random, part1_long, part2_long)
        print('Mutated')
        seq_random = sub('-', '', seq_long_random)
        seq_compl_random = sub('-', '', seq_compl_long_random)
        part1_random = sub('-', '', part1_long_random)
        part2_random = sub('-', '', part2_long_random)
        mut_coord_relative_random = mut_coord_relative
        result = [seq_random, seq_compl_random, part1_random, part2_random, mut_to_random, mut_coord_relative_random]
    elif how == "position":
        if(mut_seq == "seq"):
            possible_coords_long = IndexesOfSymbol(seq_long, mut_from)
            possible_coords = IndexesOfSymbol(seq, mut_from)
        else:
            possible_coords_long = IndexesOfSymbol(seq_compl_long, mut_from)
            possible_coords = IndexesOfSymbol(seq_compl, mut_from)
        print(possible_coords)  
        if len(possible_coords_long) == 1: # if there is only one nt of type mut_from in this handle, we don't do randomization
            result = [0, 0, 0, 0, 0, 0]
        else:
            mut_coord_relative_random = choice(list(set(possible_coords) - set([mut_coord_relative])))
            mut_coord_relative_random_long = possible_coords_long[possible_coords.index(mut_coord_relative_random)]
            print('Found random mut_coord_relative')
            seq_long_random, seq_compl_long_random, part1_long_random, part2_long_random = RecalculateEnergies.Mutate(
                seq_long, seq_compl_long, mut_coord_relative_random_long, mut_seq, mut_to, part1_long, part2_long)
            print('Mutated')
            seq_random = sub('-', '', seq_long_random)
            seq_compl_random = sub('-', '', seq_compl_long_random)
            part1_random = sub('-', '', part1_long_random)
            part2_random = sub('-', '', part2_long_random)
            mut_to_random = mut_to
            print(mut_coord_relative_random)
            result = seq_random, seq_compl_random, part1_random, part2_random, mut_to_random, mut_coord_relative_random
    return(result)


def Recalculate_one_panhandle(lock, df, need_random, title, panhandle_id_mut_info):
    Complements = {"A" : "T", "T" : "A", "G" : "C", "C" : "G"}
    with open('../out/panhandles_info_' + title + '.txt', 'a') as error_f:
        with lock:
            error_f.write(str(panhandle_id_mut_info) + '\n')
    print(panhandle_id_mut_info)
    panhandle_id = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'panhandle_id'].values[0]
    seq = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'alignment1'].values[0]
    seq_compl = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'alignment2'].values[0]
    mut_coord_relative = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'mut_coord_relative'].values[0]
    mut_seq = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'panhandle_hand'].values[0]
    mut_seq = 'seq' if mut_seq == 'left' else 'seq_compl'
    mut_coord = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'mut_coord'].values[0]
    strand = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'strand'].values[0]
    mut_from = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'mut_from'].values[0]
    mut_to = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'mut_to'].values[0]
    mut_info = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'mutation_info'].values[0]
    energy = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'energy'].values[0]
    structure = df.loc[df.panhandle_id_mut_info == panhandle_id_mut_info, 'structure'].values[0]
    print('extracted values')
    part1_long, part2_long, seq_long, seq_compl_long, mut_coord_relative_long = RecalculateEnergies.AddDelitions(
        structure, seq, seq_compl, mut_coord_relative, mut_seq)
    print('Added deletions')
    seq_long_mutated, seq_compl_long_mutated, part1_long_mutated, part2_long_mutated = RecalculateEnergies.Mutate(
        seq_long, seq_compl_long, mut_coord_relative_long, mut_seq, mut_to, part1_long, part2_long)
    seq_mutated = sub('-', '', seq_long_mutated)
    seq_compl_mutated = sub('-', '', seq_compl_long_mutated)
    part1_mutated = sub('-', '', part1_long_mutated)
    part2_mutated = sub('-', '', part2_long_mutated)
    print('Mutated')
    align = RecalculateEnergies.FindEnergyFromStructure(seq_mutated, seq_compl_mutated, part1_mutated,
                                                        part2_mutated)
    print('Finded new energy')

    results_one_panhandle_table = pd.DataFrame({'panhandle_id': panhandle_id,
                                                'mut_coord_relative': mut_coord_relative,
                                                'mut_coord': mut_coord,
                                                'new_energy': align[0],
                                                'new_structure': align[1],
                                                'new_alignment1': align[2],
                                                'new_alignment2': align[3],
                                                'mut_to': mut_to,
                                                'old_energy': energy,
                                                'mutation_info': mut_info}, index=[0])
    if need_random:
        seq_random, seq_compl_random, part1_random, part2_random, mut_to_random, mut_coord_relative_random = RandomizeMutations(
            seq, seq_compl, mut_to, mut_from, seq_long, seq_compl_long, mut_coord_relative, mut_seq, part1_long, part2_long, "position")
        if seq_random == 0:
            print("only one nt of these type!")
        else:   
            align = RecalculateEnergies.FindEnergyFromStructure(seq_random, seq_compl_random, part1_random,
                                                                part2_random)
            print('Recalculated energy')
            results_one_panhandle_table['random_energy'] = align[0]
            results_one_panhandle_table['random_structure'] = align[1]
            results_one_panhandle_table['random_alignment1'] = align[2]
            results_one_panhandle_table['random_alignment2'] = align[3]
            results_one_panhandle_table['random_mut_to'] = mut_to_random
            results_one_panhandle_table['random_mut_relative_coord'] = mut_coord_relative_random

    results_one_panhandle_table = results_one_panhandle_table.reindex(sorted(results_one_panhandle_table.columns), axis=1)

    with open('../out/out_mutator_' + title + '.tsv', 'a') as f:
        with lock:
            results_one_panhandle_table.to_csv(f,  sep='\t', index=False, header=False)
    return (0)


def FindComplementNt(column, row):
    Complements = {"A" : "T", "T" : "A", "G" : "C", "C" : "G"}
    if row['strand'] == '-':
        return(Complements[row[column]])
    else:
        return(row[column])

def CorrectStructure(row):
    if row['strand'] == '-':
        structure = row['structure']
        part1, part2 = RecalculateEnergies.SplitStructure(structure)
        part1 = part1[::-1]
        part2 = part2[::-1]
        return(part1+part2)
    else:
        return(row['structure'])

def Recalculate_energy_after_mutations(path_to_ph_with_mut, threads, need_random, title):
    df = pd.read_csv(path_to_ph_with_mut, sep='\t')
    start_time = time.time()
    # Find relative mutation coordinates in handles
    df['mut_coord_relative'] = df.apply(MakeRelativeCoords, axis=1)
    # Correct on minus strand for mut nt
    FindComplementNt_tmp = partial(FindComplementNt, 'mut_from')
    df.mut_from = df.apply(FindComplementNt_tmp, axis=1)
    FindComplementNt_tmp = partial(FindComplementNt, 'mut_to')
    df.mut_to = df.apply(FindComplementNt_tmp, axis=1)
    # Correct in minus strand for structures
    df.structure = df.apply(CorrectStructure, axis=1)
    # Make pairs mut - ph
    df['panhandle_id_mut_info'] = df['panhandle_id'].map(str) + '_' + df['mutation_info'].map(str)
    panhandle_id_mut_info = df['panhandle_id_mut_info'].values
    # Prepare outputfiles
    with open('../out/panhandles_info_' + title + '.txt', 'w') as error_f:
        error_f.write('Started alignment: \n')
        error_f.write(str(len(panhandle_id_mut_info)) + '\n')
    if need_random:
        results_one_panhandle_table = pd.DataFrame({'mut_coord': [],
                                                    'mut_coord_relative': [],
                                                    'mut_to': [],
                                                    'mutation_info': [],
                                                    'new_alignment1': [],
                                                    'new_alignment2': [],
                                                    'new_energy': [],
                                                    'new_structure': [],
                                                    'old_energy': [],
                                                    'panhandle_id': [],
                                                    'random_energy': [],
                                                    'random_structure': [],
                                                    'random_alignment1': [],
                                                    'random_alignment2': [],
                                                    'random_mut_to': [],
                                                    'random_mut_relative_coord' : []})
    else:
        results_one_panhandle_table = pd.DataFrame({'panhandle_id': [],
                                                    'mut_coord_relative': [],
                                                    'mut_coord': [],
                                                    'new_energy': [],
                                                    'new_structure': [],
                                                    'new_alignment1': [],
                                                    'new_alignment2': [],
                                                    'mut_to': [],
                                                    'old_energy': [],
                                                    'mutation_info': []})
    with open('../out/out_mutator_' + title + '.tsv', 'w') as f:
        results_one_panhandle_table.to_csv(f, sep='\t', index=False, header=True)
    p = mp.Pool(processes=threads)
    m = mp.Manager()
    lock = m.Lock()
    func = partial(Recalculate_one_panhandle, lock, df, need_random, title)
    p.map(func, panhandle_id_mut_info)
    p.close()
    p.join()

    print('recalculated energies!')
    print(time.time() - start_time)
    return (0)


def main(argv):
    path_to_ph = ""
    path_to_mut = ""
    path_to_ph_with_mut = ""
    threads = 1
    need_random = True
    title = ''
    try:
        opts, args = getopt.getopt(argv, "h:p:t:r:i:m:w:",
                                   ["help", "path_to_ph", "threads",
                                    "need_random", "title", "path_to_mut"])
    except getopt.GetoptError:
        print('RecalculateAfterMutations.py -p <path_to_ph> -t <threads> -r <need_random> -i <title> -m <path_to_mut.bed> -w <path_to_ph_with_mut>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('RecalculateAfterMutations.py -p <path_to_ph> -t <threads> -r <need_random> -i <title> -m <path_to_mut.bed> -w <path_to_ph_with_mut>')
            sys.exit()
        elif opt in ("-p", "--path_to_ph"):
            path_to_ph = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-r", "--need_random"):
            need_random = bool(arg)
        elif opt in ("-i", "--title"):
            title = arg
        elif opt in ("-m", "--path_to_mut"):
            path_to_mut = arg
        elif opt in ("-w", "--path_to_ph_with_mut"):
            path_to_ph_with_mut = arg
    if path_to_ph_with_mut == "":
        path_to_ph_with_mut = "../out/panhandles_with_mutations.tsv"
        Intersect_with_mutations(path_to_ph, path_to_mut, title)
    Recalculate_energy_after_mutations(path_to_ph_with_mut, threads, need_random, title)


if __name__ == '__main__':
    main(sys.argv[1:])




