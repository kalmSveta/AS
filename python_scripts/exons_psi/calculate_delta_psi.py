import pandas as pd
import numpy as np
import glob
import re
from functools import partial

def Parse_one_file(type_, path_to_file_with_name, donorId, cassete_exons_df):
    df = pd.read_csv(path_to_file_with_name[0], sep = '\t', header = None)
    # select only cassette exons
    #df = pd.merge(df, cassete_exons_df, left_on = [0,3,4], right_on = ['chr','exon_start','exon_stop'], how = 'inner')
    # parse exc and inc
    df['exc'] = df.loc[:,8].str.split('; ',5,expand=True).loc[:,1].str.split(' ',expand = True)[1].str.split('"',expand=True)[1]
    df['inc'] = df.loc[:,8].str.split('; ',5,expand=True).loc[:,2].str.split(' ',expand = True)[1].str.split('"',expand=True)[1]
    df = df.loc[df[2] == 'exon']
    df['exon_chr_start_end_strand'] = df[0].map(str) + ':' + df[3].map(str) + '_' + df[4].map(str) + df[6].map(str)
    df = df.loc[:,['inc','exc','exon_chr_start_end_strand']]
    df.inc = df['inc'].apply(int)
    df.exc = df['exc'].apply(int)
    # remove proteins which have too small denominator 
    df = df.loc[df['inc'] + 2 * df['exc'] >= 10]
    # rename columns to be able to merge 
    df.rename(columns={"inc": "inc_" + str(donorId) + '_' + type_, "exc": "exc_" + str(donorId) + '_' + type_}, inplace = True)
    return(df)

def process_one_tissue(tissue, metadata, path_to_files, cassete_exons_df):    
    metadata_tmp = metadata.loc[metadata['tissue'] == tissue]
    print(metadata_tmp.shape)

    cassete_exons_df = cassete_exons_df.loc[:,['chr','exon_start','exon_stop']]
    cassete_exons_df = cassete_exons_df.drop_duplicates()

    files = glob.glob(path_to_files + '*.gff') # list all files in dir
    func = partial(re.sub, path_to_files, '') 
    files = map(func, files)# remove path to files from list
    normal_donorIds = []
    for file_ in files: # list of donorIds, hich have a normal file
        if len(metadata_tmp[metadata_tmp['normal_fileName_gff'].str.contains(file_)].donorId.values) > 0:
            normal_donorIds.append(metadata_tmp[metadata_tmp['normal_fileName_gff'].str.contains(file_)].donorId.values[0])
    tumor_donorIds = []
    for file_ in files:# list of donorIds, hich have a tumor file
        if len(metadata_tmp[metadata_tmp['tumor_fileName_gff'].str.contains(file_)].donorId.values) > 0:
            tumor_donorIds.append(metadata_tmp[metadata_tmp['tumor_fileName_gff'].str.contains(file_)].donorId.values[0])
    donorIds = set(normal_donorIds) & set(tumor_donorIds)# select only donors which have both norm and tum file
    print('Have donors: ' + str(len(donorIds)))
    n_donors = len(donorIds)
    # merge tumor_normal for all donors
    # if inc+2exc <= 10 -> remove exon from tumor and normal
    
    donorIds = list(donorIds)
    donorId = donorIds[0]
    normal_file_name = path_to_files + metadata_tmp[metadata_tmp.donorId == donorId].normal_fileName_gff.values
    tumor_file_name = path_to_files + metadata_tmp[metadata_tmp.donorId == donorId].tumor_fileName_gff.values
    df_normal = Parse_one_file('normal', normal_file_name, donorId, cassete_exons_df)
    df_tumor = Parse_one_file('tumor', tumor_file_name, donorId, cassete_exons_df)
    df_global = pd.merge(df_normal, df_tumor, on = 'exon_chr_start_end_strand', how = 'inner')
    
    i = 0

    for donorId in donorIds[1:]:
        i = i + 1
        print(i, donorId)
        normal_file_name = path_to_files + metadata_tmp[metadata_tmp.donorId == donorId].normal_fileName_gff.values
        tumor_file_name = path_to_files + metadata_tmp[metadata_tmp.donorId == donorId].tumor_fileName_gff.values
        df_normal = Parse_one_file('normal', normal_file_name, donorId, cassete_exons_df)
        df_tumor = Parse_one_file('tumor', tumor_file_name, donorId, cassete_exons_df)
        df_merged = pd.merge(df_normal, df_tumor, on = 'exon_chr_start_end_strand', how = 'inner')
        df_normal = None
        df_tumor = None
        df_global = pd.merge(df_global, df_merged, on = 'exon_chr_start_end_strand', how = 'outer')
        df_merged = None
        df_global.to_csv(tissue +'_not_filtered_not_completed.csv')
	print(df_global.shape)
    print('have gotten the big df')
    df_global.to_csv(tissue +'_not_filtered_completed.csv') # inc, exc
    print('have saved it')
    # select only exons with info in at last half of the matched pairs
    df_global['#null'] = df_global.isnull().sum(axis=1)
    df_global = df_global[df_global['#null']/4 <= n_donors/2]
    # count sum and median delta psi = psi_normal - psi_tumor
    df_global = df_global.fillna(0)
    df_global['sum_exc_tumor'] = df_global.filter(like ='tumor', axis=1).filter(like = 'exc', axis = 1).sum(axis=1)
    df_global['sum_inc_tumor'] = df_global.filter(like ='tumor', axis=1).filter(like = 'inc', axis = 1).sum(axis=1)
    df_global['sum_exc_normal'] = df_global.filter(like ='normal', axis=1).filter(like = 'exc', axis = 1).sum(axis=1)
    df_global['sum_inc_normal'] = df_global.filter(like ='normal', axis=1).filter(like = 'inc', axis = 1).sum(axis=1)
    df_global['sum_psi_tumor'] = df_global['sum_inc_tumor']/(df_global['sum_inc_tumor'] + 2*df_global['sum_exc_tumor'])
    df_global['sum_psi_normal'] = df_global['sum_inc_normal']/(df_global['sum_inc_normal'] + 2*df_global['sum_exc_normal'])
    df_global['sum_delta_psi'] = df_global['sum_psi_normal'] - df_global['sum_psi_tumor']
    df_global['median_exc_tumor'] = df_global.filter(like ='tumor', axis=1).filter(like = 'exc', axis = 1).median(axis=1)
    df_global['median_inc_tumor'] = df_global.filter(like ='tumor', axis=1).filter(like = 'inc', axis = 1).median(axis=1)
    df_global['median_exc_normal'] = df_global.filter(like ='normal', axis=1).filter(like = 'exc', axis = 1).median(axis=1)
    df_global['median_inc_normal'] = df_global.filter(like ='normal', axis=1).filter(like = 'inc', axis = 1).median(axis=1)
    df_global['median_psi_tumor'] = df_global['median_inc_tumor']/(df_global['median_inc_tumor'] + 2*df_global['median_exc_tumor'])
    df_global['median_psi_normal'] = df_global['median_inc_normal']/(df_global['median_inc_normal'] + 2*df_global['median_exc_normal'])
    df_global['median_delta_psi'] = df_global['median_psi_normal'] - df_global['median_psi_tumor']
    df_global.to_csv(tissue + '2.csv')# all delta psi, all exc, inc
    df_global =  df_global.loc[:,["exon_chr_start_end_strand","median_delta_psi",'sum_delta_psi']]
    df_global = df_global[(abs(df_global["median_delta_psi"]) >= 0.05) & (abs(df_global["sum_delta_psi"]) >= 0.05)]
    df_global.to_csv(tissue + '_delta_psi2.csv') #filtered by delta psi, only delta psi
    return(0)
    
    
