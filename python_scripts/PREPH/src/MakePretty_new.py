#!/usr/bin/env python2
import pandas as pd
import sys, getopt


def MakePretty(path_to_ph, path_to_genes):
    df = pd.read_csv(path_to_ph, sep="\t")
    # Divide columns
    df[["start_gene", "end_gene"]] = df.gene.str.split('_', expand=True)[[1,2]].astype(int)
    df[["strand"]] = df.gene.str.split('_', expand=True)[[3]]
    df[["chr"]] = df.gene.str.split('_', expand=True)[[0]]
    df[["interval1_start", "interval1_end"]] = df.interval1.str.split('_', expand=True)[[1,2]].astype(int)
    df[["interval2_start", "interval2_end"]] = df.interval2.str.split('_', expand=True)[[1,2]].astype(int)

    # Make absolute coordinates
    df["panhandle_start"] = df.interval1_start + df.start_al1
    df["panhandle_left_hand"] = df.interval1_start + df.end_al1
    df["panhandle_right_hand"] = df.interval2_start + df.start_al2
    df["panhandle_end"] = df.interval2_start + df.end_al2

    # Calculate handle length
    df["al1_length"] = df.end_al1 - df.start_al1 + 1
    df["al2_length"] = df.end_al2 - df.start_al2 + 1

    # Remove broken phs (ph in one conserved interval that have end lefter than start)
    df = df.loc[df.panhandle_start < df.panhandle_right_hand]
    df = df.loc[df.panhandle_left_hand < df.panhandle_right_hand]
    df.drop(['interval1', 'interval2','start_al1','end_al1','start_al2','end_al2',
             'interval1_start','interval2_start', 'interval1_end','interval2_end', 'gene'], axis=1, inplace=True)

    # Add gene names
    anno = pd.read_csv(path_to_genes, sep='\t', header=None)
    anno[["gene_id", "gene_name"]] = anno[3].str.split('_', expand=True)[[0,1]]
    anno.drop([3], axis = 1, inplace = True)
    anno.columns = ['chr', 'start_gene', 'end_gene', 'score', 'strand', 'gene_id', 'gene_name']
    anno = anno[['chr', 'start_gene', 'end_gene', 'gene_id', 'gene_name', 'score', 'strand']]
    df = pd.merge(df, anno, on=['chr','start_gene','end_gene','strand'], how="inner")
    df.drop('score', axis=1, inplace=True)
    df.drop_duplicates(subset=["chr", "panhandle_start", "panhandle_left_hand",
                               "panhandle_right_hand", "panhandle_end"], inplace=True)

    # Make alignments 5'-3' for all strands
    df['rev1'] = df.alignment1.apply(lambda x: x[::-1])
    df['rev2'] = df.alignment2.apply(lambda x: x[::-1])
    df.loc[df.strand == '-', 'alignment1'] = df.loc[df.strand == '-'].rev1
    df.loc[df.strand == '-', 'alignment2'] = df.loc[df.strand == '-'].rev2
    df.drop(['rev1', 'rev2'], axis=1)

    # Make structures 5'-3' for all strands
    df['rev1'] = df.alignment1.apply(lambda x: x[::-1])
    df['rev2'] = df.alignment2.apply(lambda x: x[::-1])
    df.loc[df.strand == '-', 'alignment1'] = df.loc[df.strand == '-'].rev1
    df.loc[df.strand == '-', 'alignment2'] = df.loc[df.strand == '-'].rev2
    df.drop(['rev1', 'rev2'], axis=1)

    # Add unique ids to phs
    df.sort_values(by=["chr", "panhandle_start", "panhandle_left_hand",
                               "panhandle_right_hand", "panhandle_end"], inplace=True)
    df["id"] = range(1, df.shape[0] + 1)
    df.to_csv("../out/panhandles_preprocessed.tsv", sep="\t", index=False)


def main(argv):
    path_to_ph = "../out/panhandles.tsv"
    path_to_genes = '../data/coding_genes.bed'

    try:
        opts, args = getopt.getopt(argv, "h:p:g:",
                                   ["help=", "panhandles=", "genes="])
    except getopt.GetoptError:
        print('MakePretty.py -p <path_to_ph> -g <path_to_genes>')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('MakePretty.py -p <path_to_ph> -g <path_to_genes>')
            sys.exit()
        elif opt in ("-p", "--path_to_ph"):
            path_to_ph = arg
        elif opt in ("-g", "--need_genes"):
            path_to_genes = arg

    MakePretty(path_to_ph, path_to_genes)


if __name__ == '__main__':
    main(sys.argv[1:])