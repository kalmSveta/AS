# AS

PrecalculateStackingEnergeis.py -k <kmer_length> -g <gt_amount_in_kmer_max>

Need to be run once before usage. Precalculates stacking energies for kmers. 
-k <kmer_length> - is a minimal length of consicutive stacked nt pairs. Must be the same as used for ./FindPanhandles.py Recommended = 5
-g <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2

#######################################################

SelectIntervals.py -a <annotation> -c <cons_regions> -l <handle_len_min> 

Selects conservative intronic regions in coding genes (input for FindPanhandles.py)

-a <annotation> - path to annotation file in gtf format
-c <cons_regions>  - path to phastConsElements100way file
-l <handle_len_min> - min length of ph handles

##########################################################

FindPanhandles.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max>  -a <handle_len_min> -t <threads> -e <energy_max> -s <have_seqs_in_input_df> -u <need_suboptimal>
Predicts panhandles with handles in the intervals. 

where -i <intervals_df> is a tab-separated tsv file with header like this :

chr     start_interval  end_interval    strand  start_gene      end_gene
chr1    136219  		136245  		-       134901  		139379
chr1    136282  		136300  		-       134901  		139379
chr1    136330  		136355  		-       134901  		139379
chr1    136382  		136403  		-       134901  		139379

or like this:
chr     start_interval  end_interval    strand  start_gene      end_gene	sequences
chr1    136219  		136245  		-       134901  		139379		GGCTTTGATAAAAA
chr1    136282  		136300  		-       134901  		139379		TTTTTATAAAGCC
chr1    136330  		136355  		-       134901  		139379		GGCCAGCAGATGG
chr1    136382  		136403  		-       134901  		139379		TGACAAACCACAGGACACTACAC

if column 'sequences' is abscent, -g <genome.fa> must be provided. It is a fasta file with sequences which include interval sequences. E.g. a whole human genome. 'sequences' colunm will be extracted automatically

-k <kmer_length> is a minimal length of consicutive stacked nt pairs. Must be the same as used for PrecalculateStackingEnergeis.py. Recommended = 5
-p <panhandle_len_max> - amximum length of looped out region of panhandle. Recommened = 10000
-a <handle_len_min> - minimal length of handles. Recommended = 10
-t <threads> - number of threads run in parallel. 
-e <energy_max> - maximum energy i kcal/mol. Recommended = -15
-u <need_suboptimal> - if True, will try to find suboptimal structures
-d <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2

###################################################################

AssessFNRate.py -i <intervals_df> -g <genome.fa> -k <kmer_lentgh> -p <panhandle_len_max> -a <handle_len_min> -t <threads> -e <energy_max> -u <need_suboptimal> -d <gt_threshold> -N <N_seeds>

Predicts panhandles with handles in the intervals, but randomizes genes

where -i <intervals_df> is a tab-separated tsv file with header like this :

chr     start_interval  end_interval    strand  start_gene      end_gene
chr1    136219  		136245  		-       134901  		139379
chr1    136282  		136300  		-       134901  		139379
chr1    136330  		136355  		-       134901  		139379
chr1    136382  		136403  		-       134901  		139379

or like this:
chr     start_interval  end_interval    strand  start_gene      end_gene	sequences
chr1    136219  		136245  		-       134901  		139379		GGCTTTGATAAAAA
chr1    136282  		136300  		-       134901  		139379		TTTTTATAAAGCC
chr1    136330  		136355  		-       134901  		139379		GGCCAGCAGATGG
chr1    136382  		136403  		-       134901  		139379		TGACAAACCACAGGACACTACAC

if column 'sequences' is abscent, -g <genome.fa> must be provided. It is a fasta file with sequences which include interval sequences. E.g. a whole human genome. 'sequences' colunm will be extracted automatically

-k <kmer_length> is a minimal length of consicutive stacked nt pairs. Must be the same as used for PrecalculateStackingEnergeis.py. Recommended = 5
-p <panhandle_len_max> - amximum length of looped out region of panhandle. Recommened = 10000
-a <handle_len_min> - minimal length of handles. Recommended = 10
-t <threads> - number of threads run in parallel. 
-e <energy_max> - maximum energy i kcal/mol. Recommended = -15
-u <need_suboptimal> - if True, will try to find suboptimal structures
-d <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2
-N <N_seeds> - repeats randomization N time

###################################################################

fold.py -f <first_seq> -s <second_seq> -k <kmer_lentgh> -a <handle_len_min> -e <energy_max> -u <need_subopt>

Looks for panhandles in 2 sequences

-f <first_seq> -s <second_seq> - two DNA sequences in witch panhandles handles mihjt appear
-a <handle_len_min> - minimal length of handles. Recommended = 10
-k <kmer_length> is a minimal length of consicutive stacked nt pairs. Must be the same as used for PrecalculateStackingEnergeis.py. Recommended = 5
-e <energy_max> - maximum energy i kcal/mol. Recommended = -15
-u <need_suboptimal> - if True, will try to find suboptimal structures
-d <gt_amount_in_kmer_max> - maximal number of GT pairs in kmer. Recommended for k = 5 is 2

Outputs a list of structures. Each element consists of:
energy in kcal/mol, start_alignment_1, end_alignment1, start_alignment2, end_alignment2, alignment1_sequence1, alignment1_sequence2, structure_in_dor_bracket_notation

##############################################################

MakePretty.py -p <path_to_ph> -g <need_genes>
Takes as input output file from FindPanhandles.py. Makes it human-readable

-p <path_to_ph> - path to FindPanhandles.py output file
-g <need_genes> - set True, if want add gene id and gene name from annotation file

################################################################

RecalculateAfterMutations.py -p <path_to_ph> -t <threads> -r <need_random> -i <title> -m <path_to_mut> -w <path_to_ph_with_mut>

Takes file with panhandles and coords of mutations as input. Intersects panhandle handles and mutations. Recalculates energies of ph after mutations introduction. Recalculates energies of the ph with mutations in the same positions but in random nt (from 2 not taken). 
Alternatevely can take file with intersected panhandles and mutations as input and just perform recalculation of energies and shuffling.

-p <path_to_ph> - path to input table with panahndles
-m <path_to_mut> - path to mutations file
OR
-w <path_to_ph_with_mut> - path to file with panhandles and mutations

-t <threads> - number of threads to run in parallel
-r <need_random> - if need recalculating energies for random mutations
-i <title> - additional title for the files 
  
################################################################
  
./Make_big_bed.R ../python_scripts/folding_pretty_copy/out/panhandles_preprocessed.tsv 'all_panhandles'

Makes bed and bigbed from standard panhandle file

################################################################

./Make_handles.R

Makes bed4 from all (left AND right) handles

#################################################################

./TestMutDensity.R path.to.gr1 path.togr2 path.to.mut

Performs Puasson test for mutations density enrichment in target intervas versus population of intervals

#################################################################

