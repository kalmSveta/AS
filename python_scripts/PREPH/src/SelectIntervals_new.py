#!/usr/bin/env python2
import pandas as pd
import sys, getopt
from subprocess import call
from functools import partial


def select_coding_genes(path_to_anno):
    headers = 0
    with open(path_to_anno) as f:
        while 1:
            line = f.readline()
            if line.startswith("#"):
                headers += 1
            else:
                break
    df = pd.read_csv(path_to_anno, sep="\t", skiprows=headers, header=None)
    df_CDS = df.loc[df[2] == 'CDS']
    coding_genes = df_CDS[8].str.split("; ", n=2, expand=True)[0]
    gene_id_place = [i for i, s in enumerate(df.iloc[1, 8].split(";")) if 'gene_id' in s][0]
    df['gene_id'] = df[8].str.split("; ", n=2, expand=True)[0]
    df_coding = df.loc[df.gene_id.isin(coding_genes.unique())]
    df_coding = df_coding.loc[df_coding[2] == "gene"]
    gene_name_place = [i for i, s in enumerate(df_coding.iloc[1, 8].split(";")) if 'gene_name' in s][0]
    df_coding['gene_name'] = df_coding[8].str.split("; ", n=5, expand=True)[gene_name_place]
    df_coding['gene_name'] = df_coding["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    df_coding['gene_id'] = df_coding["gene_id"].str.split(' "', n=3, expand=True)[1].str[:-1]
    df_coding['gene_id_name'] = df_coding['gene_id'] + '_' + df_coding['gene_name']
    df_coding_bed = df_coding[[0, 3, 4, 'gene_id_name', 5, 6]]
    df_coding_bed.to_csv("../data/coding_genes.bed", sep="\t", index=False, header=False)
    return(0)


def select_intronic_regions(path_to_anno, flank_length):
    headers = 0
    with open(path_to_anno) as f:
        while 1:
            line = f.readline()
            if line.startswith("#"):
                headers += 1
            else:
                break
    df = pd.read_csv(path_to_anno, sep="\t", skiprows=headers, header=None)
    # select CDSs
    CDS = df[df[2] == 'CDS']

    CDS['gene_id'] = CDS[8].str.split("; ", n=2, expand=True)[0]
    CDS['gene_id'] = CDS.gene_id.replace("gene_id ", "", regex=True)
    CDS['gene_id'] = CDS.gene_id.replace("\"", "", regex=True)  
    gene_name_place = [i for i, s in enumerate(CDS.iloc[1, 8].split(";")) if 'gene_name' in s][0]
    CDS['gene_name'] = CDS[8].str.split("; ", n=5, expand=True)[gene_name_place]
    CDS['gene_name'] = CDS["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    CDS['gene_id_name'] = CDS['gene_id'] + '_' + CDS['gene_name']
    CDS = CDS[[0, 3, 4, 'gene_id_name', 5, 6]]
    CDS.to_csv("../data/CDS.bed", sep="\t", index=False, header=False)

    # select genes
    genes = df[df[2] == 'gene']
    genes['gene_id'] = genes[8].str.split("; ", n=2, expand=True)[0]
    genes['gene_id'] = genes.gene_id.replace("gene_id ", "", regex=True)  
    genes['gene_id'] = genes.gene_id.replace("\"", "", regex=True)    
    gene_name_place = [i for i, s in enumerate(genes.iloc[1, 8].split(";")) if 'gene_name' in s][0]
    genes['gene_name'] = genes[8].str.split("; ", n=5, expand=True)[gene_name_place]
    genes['gene_name'] = genes["gene_name"].str.split(' "', n=3, expand=True)[1].str[:-1]
    genes['gene_id_name'] = genes['gene_id'] + '_' + genes['gene_name']
    genes = genes[[0, 3, 4, 'gene_id_name', 5, 6]]
    genes.to_csv("../data/genes.bed", sep="\t", index=False, header=False)

    # subtract CDS from genes
    x = 'bedtools subtract -a ../data/genes.bed -b ../data/CDS.bed -s > ../data/genes_no_CDS.bed' 
    call([x], shell = True)
    dt = pd.read_csv('../data/genes_no_CDS.bed',  sep="\t", header=None)
    dt = dt.sort_values(by =[3, 1, 2])

    # add flanks
    dt[1] = dt[1] - flank_length
    dt[2] = dt[2] + flank_length
    dt.to_csv('../data/tmp', sep="\t", index=False, header=False)
    call(['bedtools intersect -a ../data/tmp -b ../data/genes.bed -s > ../data/introns_python.bed'], shell = True)
    call(['rm ../data/genes_no_CDS.bed'], shell = True)
    call(['rm ../data/CDS.bed'], shell = True)
    call(['rm ../data/genes.bed'], shell = True)
    call(['rm ../data/tmp'], shell = True)
    return(0)


def main(argv):
    path_to_anno="../../../conservative_features/gencode.v19.annotation.gtf"
    path_to_cons="../../../conservative_features/phastConsElements100way.txt"
    length_min=10

    try:
        opts, args = getopt.getopt(argv, "h:a:c:l:i:r:f:",
                                   ["help=", "annotation=", "cons_regions=", "handle_len_min=", "always_intronic=", "chromosomes=", "flanks="])
    except getopt.GetoptError:
        print('SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -f <flanks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -f <flanks>')
            sys.exit()
        elif opt in ("-a", "--annotation"):
            path_to_anno = arg
        elif opt in ("-c", "--cons_regions"):
            path_to_cons = arg
        elif opt in ("-l", "--handle_len_min"):
            length_min = arg
        elif opt in ("-f", "--flanks"):
            flank_length = int(arg)

    select_coding_genes(path_to_anno)
    print('selected coding genes')
    select_intronic_regions(path_to_anno, flank_length)
    call(['./Select_conins_new.sh', path_to_cons, length_min])
    print('selected intervals')

if __name__ == '__main__':
    main(sys.argv[1:])
