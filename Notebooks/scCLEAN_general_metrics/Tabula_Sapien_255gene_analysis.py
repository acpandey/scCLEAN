# Jon Bezney
# processing for scCLEAN - Tabula Sapien 255 gene analysis 
# Reduce the count matrices down to only the 255 genes removed with scCLEAN
# this was ran for each of the 4 types of cell types: Immune, stromal, epithelial, endothelial 

import scanpy as sc
import pandas as pd
import numpy as np 

#read in the 255 protein coding genes targeted with scCLEAN
NVG = pd.read_csv('NVG_scCLEAN.csv', header=None)
NVG = list(NVG[0])

MT = pd.read_csv('MT_scCLEAN.csv', header=None)
MT = list(MT[0])

Ribo = pd.read_csv('Ribo_scCLEAN.csv', header=None)
Ribo = list(Ribo[0])

scclean_genes = NVG+MT+Ribo

#read in the tabula sapien single cell data
#The counts in the attached X matrix are normalized
adata = sc.read_h5ad('Immune_tabula_sapien.h5ad')

#change the index to the gene names and not ensemble IDs
gene_df = pd.DataFrame(adata.var)
gene_df = gene_df.set_index('feature_name')
adata.var = gene_df

#delete the stuff you don't need, keep obs and var dataframes
del adata.uns
del adata.obsm
del adata.layers
del adata.obsp

#round the count matrices to the first decimal place
adata.layers['X_round'] = round(adata.X.copy(),1)
adata.X = adata.layers['X_round']
del adata.layers['X_round']

#keep the OBS of interest
obs_df = pd.DataFrame(adata.obs)
obs_df = obs_df[['cell_type']]
del adata.obs
adata.obs['cell_type'] = obs_df['cell_type']

#keep the var of interest
var_df = pd.DataFrame(adata.var)
var_df = var_df[['dispersions_norm', 'means', 'dispersions']]
del adata.var
adata.var['dispersions_norm'] = var_df['dispersions_norm']
adata.var['means'] = var_df['means']
adata.var['dispersions'] = var_df['dispersions']

#save the information needed
var_df.to_csv('Immune_VAR_tabula_sapien_full.csv')
obs_df.to_csv('Immune_OBS_tabula_sapien_full.csv')

#subset the data matrix to only the genes of interest
#not valid ['SLC25A6', 'H3F3B', 'H3F3A', 'H2AFZ', 'ATP5MPL']

not_found = ['ATP5MPL', 'H2AFZ', 'H3F3A', 'H3F3B', 'SLC25A6']
scclean_rev = [name for name in scclean_genes if name not in not_found]

adata_sub = adata[:, scclean_rev]

#save as h5ad - this is the full not subsetted data 
t=adata_sub.X.toarray()
X_csv = pd.DataFrame(data=t, index=adata_sub.obs_names, columns=adata_sub.var_names)
X_csv.to_csv('Immune_tabula_sapien_255.csv')


