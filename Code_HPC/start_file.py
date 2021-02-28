## Code author SP



import os


#import matplotlib
#matplotlib.use('Agg')


############
#### update path variables for data
############


path_data = '/Users/Downloads/Jeff_MERFISH_data' 
path_t = os.getcwd() 
path_t = '/Users/Downloads/Jeff_MERFISH_data/current'
path_python = path_t

path_scseq = os.path.join(path_data + '/GSE113576_matrix.mtx')
path_cb = os.path.join(path_data + '/M22E1_codebook.csv')
path_temp_2 = os.path.join(path_data + '/GSE113576_genes.tsv')

path_Allen = os.path.join(path_data+'/SongLinROIS_deduplicated.json')
brian_csv_path = os.path.join(path_data+'/smFISH_MCT_CZI_Panel_0_spot_table.csv')

#read in the cluster ids
path_sc_cellid = os.path.join(path_data + '/aau5324_Moffitt_Table-S1.xlsx')

#read in Allen data .csvs
path_Allen_data = '/home/Merfish_model/Merfish/Jeff_MERFISH_data/Allen_data'


exec(open(os.path.join(path_python+"/init_file.py")).read())
