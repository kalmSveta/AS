#!/usr/bin/env python2
import pandas as pd
import sys, getopt
from subprocess import call


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
    df_coding_bed = df_coding[[0, 3, 4, 'gene_id', 'gene_name', 5, 6]]
    df_coding_bed.to_csv("../data/coding_genes.bed", sep="\t", index=False, header=False)
    return(0)


def main(argv):
    path_to_anno = "../conservative_features/gencode.v19.annotation.gtf"
    path_to_cons = "../../phastConsElements100way.txt"
    path_to_CDS = ""
    always_intronic = False
    length_min = 10
    chromosomes = ''

    try:
        opts, args = getopt.getopt(argv, "h:a:c:l:i:r:",
                                   ["help=", "annotation=", "cons_regions=", "handle_len_min=", "always_intronic=", "chromosomes="])
    except getopt.GetoptError:
        print('SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -r <chromosomes>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> -r <chromosomes>')
            sys.exit()
        elif opt in ("-a", "--annotation"):
            path_to_anno = arg
        elif opt in ("-c", "--cons_regions"):
            path_to_cons = arg
        elif opt in ("-l", "--handle_len_min"):
            length_min = arg
        elif opt in ("-i", "--always_intronic"):
            always_intronic = arg
        elif opt in ("-r", "--chromosomes"):
            chromosomes = arg


    select_coding_genes(path_to_anno)
    print('selected coding genes')
    call(['./Select_conins.sh', path_to_anno, path_to_cons, length_min, always_intronic, chromosomes])
    print('selected intervals')

if __name__ == '__main__':
    main(sys.argv[1:])
