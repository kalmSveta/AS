#!/usr/bin/env python2

import sys, getopt
from subprocess import call
import pandas as pd
sys.path.insert(0, './')
import time
from functools import partial
import multiprocessing as mp
import glob
import re
import numpy as np


Dic_bp_plus = {'CG':0, 'GC':1, 'GT':2, 'TG':3,'AT':4,'TA':5}
Dic_bp_min = {'CG':0, 'GC':1, 'CA':2, 'AC':3,'AT':4,'TA':5}

def ReverseCompl(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    rc_seq = ''.join(bases)
    return(rc_seq)


def SplitStructure(structure):
    structure1 = structure[ : structure.find(')')]
    structure2 = structure[structure.find(')') : ]
    return(structure1, structure2)


def AddDelitions(structure, seq, seq_compl, mut_coord = -1, mut_seq = '',  have_mut = False):
    part1, part2 = SplitStructure(structure)
    i = 0
    while i != max(len(part1), len(part2)):
        if part1[len(part1) - i - 1] == '.' and part2[i] != '.':
            if have_mut:
                if mut_seq == 'seq_compl' and mut_coord >= i:
                    mut_coord += 1
            part2 = part2[ : i] + '-' + part2[i : ]
            seq_compl = seq_compl[ : i] + '-' + seq_compl[i : ]
        elif part2[i] == '.' and part1[len(part1) - i - 1] != '.':
            if have_mut:
                if mut_seq == 'seq' and mut_coord > len(part1) - i - 1:
                    mut_coord += 1
            part1 = part1[ : len(part1) - i] + '-' + part1[len(part1) - i : ]
            seq = seq[ : len(seq) - i] + '-' + seq[len(seq) - i : ]
        i += 1
    return(part1, part2, seq, seq_compl, mut_coord)


def AddDeletionsOne(lock, df, outputfile, panhandle_id):
    with open('../out/panhandles_info.txt', 'a') as error_f:
        error_f.write(str(panhandle_id) + '\n')
    seq = df.loc[df.id == panhandle_id,'alignment1'].values[0]
    seq_compl = df.loc[df.id == panhandle_id,'alignment2'].values[0]
    structure = df.loc[df.id == panhandle_id,'structure'].values[0]
    chr = df.loc[df.id == panhandle_id,'chr'].values[0]
    strand = df.loc[df.id == panhandle_id,'strand'].values[0]
    panhandle_start = df.loc[df.id == panhandle_id,'panhandle_start'].values[0]
    panhandle_left_hand = df.loc[df.id == panhandle_id,'panhandle_left_hand'].values[0]
    panhandle_right_hand = df.loc[df.id == panhandle_id,'panhandle_right_hand'].values[0]
    panhandle_end = df.loc[df.id == panhandle_id,'panhandle_end'].values[0]
    part1_long, part2_long, seq_long, seq_compl_long, mut_coord_relative = AddDelitions(structure, seq, seq_compl, have_mut = False)
    results_one_panhandle_table = pd.DataFrame({'panhandle_id':panhandle_id,
                                                'struct_part1':part1_long, 
                                                'struct_part2':part2_long, 
                                                'alignment1':seq_long, 
                                                'alignment2':seq_compl_long, 
                                                'chr':chr, 
                                                'panhandle_start':panhandle_start, 
                                                'panhandle_left_hand':panhandle_left_hand, 
                                                'panhandle_right_hand':panhandle_right_hand, 
                                                'panhandle_end':panhandle_end,
                                                'strand':strand}, index=[0])
    with open(outputfile, 'a') as f:
        with lock:
            results_one_panhandle_table.to_csv(f, sep="\t", index=False)


def MakeNtPairs(lock, df, outputfile, panhandle_id):
    with open('../out/panhandles_info.txt', 'a') as error_f:
        error_f.write(str(panhandle_id) + '\n')
    seq = df.loc[df.panhandle_id == panhandle_id,'alignment1'].values[0]
    seq_compl = df.loc[df.panhandle_id == panhandle_id,'alignment2'].values[0]
    chr_ = df.loc[df.panhandle_id == panhandle_id,'chr'].values[0]
    panhandle_start = float(df.loc[df.panhandle_id == panhandle_id,'panhandle_start'].values[0])
    panhandle_end = float(df.loc[df.panhandle_id == panhandle_id,'panhandle_end'].values[0])
    strand = df.loc[df.panhandle_id == panhandle_id,'strand'].values[0]
    results = list()
    add1 = 0
    add2 = 0
    seq_compl = seq_compl[::-1]
    for i, nt1, nt2 in zip(range(len(seq)), seq, seq_compl):
        print(i,  nt1,  nt2)
        if ((nt1 + nt2 in Dic_bp_plus) & (strand == '+')) | ((nt1 + nt2 in Dic_bp_min) & (strand == '-')):
            name1 = str(panhandle_id) + '_' + seq[i]
            name2 = str(panhandle_id) + '_' + seq_compl[i]
            start1 = panhandle_start + i + add1
            start2 = panhandle_end - i + add2
            strand1 = strand
            strand2 = strand
            results.append(pd.DataFrame([[chr_, start1, start1, name1, '.', strand1, chr_, start2, start2, name2, '.', strand2]], columns = ['chr1', 'start1', 'end1', 'name1', 'score1', 'strand1', 'chr2', 'start2', 'end2', 'name2', 'score2', 'strand2']))
        else:
            if nt1 == '-':
                add1 -= 1
            elif nt2 == '-':
                add2 += 1
    results = pd.concat(results)
    with open(outputfile, 'a') as f:
        with lock:
            results.to_csv(f,  sep='\t',  index=False)


def CountCompensatory(nts,  mut_path,  threshold, save_df = False):
    nts['compensatory'] = False
    # select pairs with mirror mutations:
    nts = nts.loc[(nts.mut_id1.notnull()) & (nts.mut_id2.notnull())]
    nts['mut_to1'] = nts.mut_to1.astype(str)
    nts['mut_to2'] = nts.mut_to2.astype(str)
    nts = nts[~(nts['mut_to1'].str.contains(',')) & ~(nts['mut_to2'].str.contains(','))]
    x = nts.mut_to1 + nts.mut_to2
    nts = nts.loc[((nts.strand1 == '+') & (x.isin(Dic_bp_plus))) | ((nts.strand1 == '-') & (x.isin(Dic_bp_min)))]
    nts = nts.loc[(nts.n_donors1 > threshold) & (nts.n_donors2 > threshold)]
    files_with_mut = glob.glob(mut_path + "*_out.txt")
    for chr_ in nts.chr1.unique():
        print(chr_)
        tmp = nts.loc[nts.chr1 == chr_]
        mut = filter(re.compile(chr_+'\.').findall, files_with_mut)[0]
        mut = pd.read_csv(mut, sep='\t',  header=None)
        for i in tmp.index:
            mut_id1 = tmp.loc[tmp.index == i].mut_id1.values[0]
            mut_tmp = mut.loc[mut[8] == mut_id1]
            if mut_tmp.empty:
                donors_left1 = []
                donors_right1 = []
            else:
                mut_tmp = mut_tmp.iloc[0]
                donors = mut_tmp.squeeze()
                donors = donors.astype(str)
                donors_left1 = list(donors[donors.str.contains('1\|')].index)
                donors_right1 = list(donors[donors.str.contains('\|1')].index)
            mut_id2 = tmp.loc[tmp.index == i].mut_id2.values[0]
            mut_tmp = mut.loc[mut[8] == mut_id2]
            if mut_tmp.empty:
                donors_left2 = []
                donors_right2 = []
            else:        
                mut_tmp = mut_tmp.iloc[0]
                donors = mut_tmp.squeeze()
                donors = donors.astype(str)
                donors_left2 = list(donors[donors.str.contains('1\|')].index)
                donors_right2 = list(donors[donors.str.contains('\|1')].index)
            donors_compensatory = list(set(donors_right1) & set(donors_right2) | set(donors_left1) & set(donors_left2))
            if len(donors_compensatory) > threshold:
                nts.loc[(nts.chr1 == chr_) & (nts.index == i), 'compensatory'] = True
    if save_df:
        nts.to_csv('../out/nt_pairs_with_compensatory.stv', sep = '\t', index = False)
    nts = nts.loc[nts.compensatory == True]
    npairs_with_compensatory_SNP = nts.shape[0]
    return(npairs_with_compensatory_SNP)


def Shuffle(nts, seed, shuffle_on):
    if shuffle_on == 'type':
    # shuffle inside same pair type
        nts_shuffled = nts.groupby(["nt2", "chr2"])['name2'].transform(np.random.RandomState(seed=seed).permutation) 
    elif shuffle_on == 'n_donors':
        # shuffle inside same pair type and same number of donors with mutation (0, 1, 2 ...1672)
        nts_shuffled = nts.groupby(["nt2", "chr2", "n_donors2"])['name2'].transform(
            np.random.RandomState(seed=seed).permutation) 
    elif shuffle_on == 'n_donors_bin':
        # shuffle inside same pair type and same bin of number of donors 
        data = nts.n_donors2.astype(int)    
        bins = range(min(data), max(data)+1, 25)[1:]
        nts['n_donors2_bin'] = np.digitize(data, bins, right=True)
        nts_shuffled = nts.groupby(["nt2", "chr2", "n_donors2_bin"])['name2'].transform(
            np.random.RandomState(seed=seed).permutation) 
    elif shuffle_on == 'ph':
        # shuffle inside same pair type inside one panhandle 
        nts_shuffled = nts.groupby(["nt2", "chr2", "ph2"])['name2'].transform(
            np.random.RandomState(seed=seed).permutation) 
    else:
        print('I do not know such on option')
    nts_unchanged1 = nts[['chr1', 'start1', 'end1', 'name1', 'score1', 'strand1', 'mut_id1', 'mut_from1', 'mut_to1', 'n_donors1']]
    nts_unchanged2 = nts[['start2', 'end2', 'name2', 'score2', 'strand2', 'mut_id2', 'mut_from2', 'mut_to2', 'n_donors2']]
    nts_new = pd.concat([nts_unchanged1.reset_index(drop=True), nts_shuffled], axis=1)
    nts_merged = nts_new.merge(nts_unchanged2, on = 'name2', how = 'inner')
    return(nts_merged)


def ShuffleAnd_Count_par(lock, nts, outputfile_comp, outputfile_pairs, mut_path, threshold, shuffle_on, iterator):
    print('Im here')
    nts_shuffled = Shuffle(nts, seed=iterator, shuffle_on = shuffle_on)
    n_pairs_with_SNP_shuffled = nts_shuffled.loc[(nts_shuffled.n_donors1 > threshold) | (nts_shuffled.n_donors2 > threshold)].shape[0]
    n_pairs_with_compensatory_SNP_shuffled = CountCompensatory(nts_shuffled, mut_path, threshold, save_df=False)
    with open(outputfile_comp, 'a') as f:
        with lock:
            f.write("%s\n" % n_pairs_with_compensatory_SNP_shuffled)
    with open(outputfile_pairs, 'a') as f:
        with lock:
            f.write("%s\n" % n_pairs_with_SNP_shuffled)
    print(iterator)
    return(0)

def ShuffleAnd_Count(nts, outputfile_comp, outputfile_pairs, mut_path, threshold, shuffle_on, iterator):
    print('Im here')
    nts_shuffled = Shuffle(nts, seed=iterator, shuffle_on = shuffle_on)
    n_pairs_with_SNP_shuffled = nts_shuffled.loc[(nts_shuffled.n_donors1 > threshold) | (nts_shuffled.n_donors2 > threshold)].shape[0]
    n_pairs_with_compensatory_SNP_shuffled = CountCompensatory(nts_shuffled, mut_path, threshold, save_df=False)
    with open(outputfile_comp, 'a') as f:
        f.write("%s\n" % n_pairs_with_compensatory_SNP_shuffled)
    with open(outputfile_pairs, 'a') as f:
        f.write("%s\n" % n_pairs_with_SNP_shuffled)
    print(iterator)
    return(0)    


def main(argv):
    threshold = 25
    N_permutations = 100
    threads = 50
    ph_input = '../../folding_pretty_copy/out/hg19_ss_flanks/panhandles_preprocessed_filtered.tsv'
    mut_path = "../../../1000genomes/new/hg19_ss_flanks/phs/"
    shuffle_on = 'n_donors_bin'
    try:
        opts, args = getopt.getopt(argv,"h:p:T:m:t:N:s:", ["help=","ph_input=","threads=","mut_path","threshold=","Npermutations", "shuffle_on"])
    except getopt.GetoptError:
        print('CalculateCompensatory.py -p <ph_input_path> -T <n_threads> -m <mut_path> -t <SNP_threshold> -N <N_perumations> -s shuffle_on')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('CalculateCompensatory.py2 -p <ph_input_path> -T <n_threads> -m <mut_path> -t <SNP_threshold> -N <N_perumations> -s shuffle_on')
            sys.exit()
        elif opt in ("-p", "--ph_input"):
            ph_input = arg
        elif opt in ("-T", "--threads"):
            threads = int(arg)
        elif opt in ("-m", "--mut_path"):
            mut_path = arg
        elif opt in ("-t", "--threshold"):
            threshold = int(arg)
        elif opt in ("-N", "--Npermutations"):
            N_permutations = int(arg)
        elif opt in ("-s", "--shuffle_on"):
            shuffle_on = arg
    ###################### make everything plus stranded 
    ph = pd.read_csv(ph_input, sep = '\t') 
    ph.loc[ph.strand=='-','alignment1'] = map(ReverseCompl, ph.loc[ph.strand == '-'].alignment1)
    ph.loc[ph.strand=='-','alignment2'] = map(ReverseCompl, ph.loc[ph.strand == '-'].alignment2)
    ###################### add deletions
    with open('../out/panhandles_info.txt', 'w') as error_f:
        error_f.write('Started alignment: \n')
        error_f.close()
    with open('../out/ph_with_deletions.tsv', 'w') as f:
        f.write('')
        f.close()
    p = mp.Pool(processes = threads)
    func = partial(AddDeletionsOne, mp.Manager().Lock(), ph, '../out/ph_with_deletions.tsv')
    panhandle_ids = ph['id'].values
    p.map(func, panhandle_ids)
    p.close()
    p.join()
    print('added deletions!')
    ph = pd.read_csv('../out/ph_with_deletions.tsv', sep="\t")
    ph = ph[ph.alignment1 != 'alignment1']
    ph.to_csv('../out/ph_with_deletions.tsv', sep="\t")
    ph["panhandle_id"] = pd.to_numeric(ph["panhandle_id"])
    #################  make pairs of nts
    with open('../out/panhandles_info.txt', 'w') as error_f:
        error_f.write('Started alignment: \n')
    with open('../out/nt_pairs.tsv', 'w') as f:
        f.write('')
    p = mp.Pool(processes = threads)
    func = partial(MakeNtPairs, mp.Manager().Lock(), ph, '../out/nt_pairs.tsv')
    panhandle_ids = ph['panhandle_id'].values
    p.map(func, panhandle_ids)
    p.close()
    p.join()
    results_df = pd.read_csv('../out/nt_pairs.tsv',  sep = '\t')
    results_df = results_df[results_df.start1 != 'start1']
    results_df['own_id'] = range(len(results_df.index))
    results_df['name1'] = results_df.name1.astype(str) + '_' + results_df.own_id.astype(str)
    results_df['name2'] = results_df.name2.astype(str) + '_' + results_df.own_id.astype(str)
    results_df.drop(columns=['own_id'],  inplace = True)
    results_df.to_csv( '../out/nt_pairs.tsv',  sep='\t',  index=False)
    ############################# select mutations with # donors > THR1, count pairs with mut
    call(['./count_donors.sh',  mut_path])
    nts = pd.read_csv( '../out/nt_pairs.tsv',  sep = '\t')
    with open('../out/nt_pairs_with_mut.tsv', 'w') as f:
        f.write('')
    files_with_mut = glob.glob(mut_path + "/*.2")
    for file in files_with_mut:
        chr_ = re.findall("chr[0-9,X,Y,MT,M,mt,x,y]+", file)[0]
        print(chr_)
        merged = nts.loc[nts.chr1 == chr_]
        try:
            mut = pd.read_csv(file,  sep = '\t', header = None)
            mut.drop_duplicates(inplace = True)
            mut.columns = ['chr','mut_coord','mut_id','mut_from','mut_to','n_donors']
            mut['chr'] = 'chr' + mut.chr.astype(str)
            merged = merged.merge(mut, how = 'left', left_on = ['chr1', 'start1'], right_on = ['chr','mut_coord'])
            merged.drop(columns=['chr', 'mut_coord'], inplace = True)
            merged.rename(columns={"mut_id": "mut_id1", "mut_from": "mut_from1", "mut_to":"mut_to1", "n_donors":"n_donors1"}, inplace = True)
            merged = merged.merge(mut, how = 'left', left_on = ['chr2', 'start2'], right_on = ['chr','mut_coord'])
            merged.drop(columns=['chr', 'mut_coord'], inplace = True)
            merged.rename(columns={"mut_id": "mut_id2", "mut_from": "mut_from2", "mut_to":"mut_to2", "n_donors":"n_donors2"}, inplace = True)
            merged.to_csv('../out/nt_pairs_with_mut.tsv',  sep = '\t',   index = False,  na_rep = 'NA',  mode='a')   
        except:
            print('nothing to parse!')
            continue  
    nts = pd.read_csv( '../out/nt_pairs_with_mut.tsv',  sep = '\t') # all nts pairs, some have mut, some have paired mut
    nts = nts.loc[nts.mut_to1 != 'mut_to1']
    nts['n_donors1'] = pd.to_numeric(nts['n_donors1'])
    nts['n_donors2'] = pd.to_numeric(nts['n_donors2'])
    nts.drop_duplicates(inplace = True)
    npairs_with_SNP = nts.loc[(nts.n_donors1 > threshold) | (nts.n_donors2 > threshold)].shape[0]
    print(npairs_with_SNP)
    nts.to_csv('../out/nt_pairs_with_mut.tsv',  sep = '\t',  index = False)
    with open('../out/compensatory_info.txt', 'w') as f:
        f.write('npairs_with_SNP = ' + str(npairs_with_SNP) + '\n')
    ########## select compensatory and count them
    nts = pd.read_csv( '../out/nt_pairs_with_mut.tsv',  sep = '\t')
    npairs_with_compensatory_SNP = CountCompensatory(nts,  mut_path,  threshold, save_df = True)
    print(npairs_with_compensatory_SNP)
    with open('../out/compensatory_info.txt', 'a') as f:
        f.write('npairs_with_compensatory_SNP = ' + str(npairs_with_compensatory_SNP) + '\n')
    #########################  shuffling
    nts = pd.read_csv( '../out/nt_pairs_with_mut.tsv',  sep = '\t')
    nts["nt2"] = nts["name2"].str.split("_", n=2, expand=True)[1]
    nts['n_donors2'].fillna(0, inplace=True)
    if shuffle_on == 'ph':
        nts["ph2"] = nts["name2"].str.split("_", n=2, expand=True)[0]
        nts.sort_values(by=["nt2", "chr2", "ph2"], inplace=True)
    else:
        nts.sort_values(by=["nt2", "chr2", "n_donors2"], inplace=True)
    with open('../out/N_compensatory_shuffled.txt', 'w') as f:
        f.write('')
    with open('../out/N_pairs_with_SNP_shuffled.txt', 'w') as f:
        f.write('')
    parallel = False
    if parallel:
        p = mp.Pool(processes = threads)
        func = partial(ShuffleAnd_Count_par, mp.Manager().Lock(), nts, '../out/N_compensatory_shuffled.txt',
                       '../out/N_pairs_with_SNP_shuffled.txt', mut_path, threshold, shuffle_on)
        p.map(func, range(N_permutations))
        p.close()
        p.join()
    else:
        for i in range(N_permutations):
            ShuffleAnd_Count(nts, '../out/N_compensatory_shuffled.txt',
                       '../out/N_pairs_with_SNP_shuffled.txt', mut_path, threshold, shuffle_on, i)


    
if __name__ == '__main__':
    main(sys.argv[1:])

