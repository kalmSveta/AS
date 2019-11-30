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



def Make_handles(ph_input):
    ph = pd.read_csv(ph_input, sep="\t")
    print(ph.head())
    ph_left = ph[['chr', 'panhandle_start', 'panhandle_left_hand', 'strand', 'id']]
    ph_right = ph[['chr', 'panhandle_right_hand', 'panhandle_end', 'strand', 'id']]
    ph_left.columns = ['chr', 'start', 'stop', 'strand', 'name']
    ph_right.columns = ['chr', 'start', 'stop', 'strand', 'name']
    ph_left['name'] = ph_left.name.apply(str) + "_left"
    ph_right['name'] = ph_right.name.apply(str) + "_right"
    ph_both = ph_left.append(ph_right, ignore_index=True)
    ph_both.chr = ph_both.chr.str[3:]
    ph_both['start'] = ph_both.start.apply(int)
    ph_both['stop'] = ph_both.stop.apply(int)
    ph_both.to_csv('../out/panhandle_handles.tsv', sep="\t", index=False, header=False)
    call(["bedtools sort", "-i", "../out/panhandle_handles.tsv"], stdout = open("../out/panhandle_handles_sorted.tsv", "w")) # NOT WORK
    call("mv ../out/panhandle_handles_sorted.tsv ../out/panhandle_handles.tsv")

def intersect(mut_path, file_):
    call(['./intersect_mut_and_panh.sh', file_, mut_path, '../out/panhandle_handles.tsv'])

def main(argv):
    ph_input = '../../folding_pretty_copy/out/foo.tsv'
    mut_path = "../out/handles_and_mut/"
    threads = 10
    try:
        opts, args = getopt.getopt(argv, "h:p:m:t:",
                                   ["help=", "ph_input=", "mut_path", "threads"])
    except getopt.GetoptError:
        print(
            'Intersect_with_mut.py -p <ph_input_path> -m <mut_path> -t <threads>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print(
                'Intersect_with_mut.py -p <ph_input_path> -m <mut_path> -t <threads>')
            sys.exit()

        elif opt in ("-p", "--ph_input"):
            ph_input = arg
        elif opt in ("-m", "--mut_path"):
            mut_path = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)

    Make_handles(ph_input)
    print("Made handles!")
    print("start to intersect...")
    files_with_mut = glob.glob(mut_path + "*.vcf.gz")
    print(files_with_mut)

    #p = mp.Pool(processes=threads)
    #func = partial(intersect, mut_path)
    #p.map(func, files_with_mut)
    #p.close()
    #p.join()
    call(['./intersect_mut_and_panh2.sh', mut_path, '../out/panhandle_handles.tsv'])



if __name__ == '__main__':
    main(sys.argv[1:])
