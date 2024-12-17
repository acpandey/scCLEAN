
###################################################################
#Script Name    :  scVI - scArches reference PBMC model                                                                                     
#Description    :  across all cell types, Control comparison                                                              
#Args               :                                                                                           
#Author         :  Jon Bezney                                                
#Email          :  jbezney@stanford.edu                                       
###################################################################

import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
import torch
import scarches as sca

#import the query control information 
#train the query model relative to the reference model 
#the scCLEAN model needs to be run first because the reference model has already been generated

#the count matrix is alread normalized 
rare = sc.read_h5ad('/oak/stanford/groups/smontgom/jbezney/scCLEAN/scArches/rare_pdc_pbmc.h5ad')
#subset to only the PBMC and not the bone marrow samples
rare = rare[rare.obs['tissue']=='PBMCs']


#here is the list of genes not found in both datasets
#need to be removed from both datasets
genes_na = ['IGLL5', 'FAM129C', 'FAM198B', 'PLA2G16', 'KIAA1468', 'SEPT4', 'TMEM56', 'FAM129A', 
 'COL4A3BP', 'CTGF', 'FAM84B', 'HRASLS2', 'FAM45A', 'PQLC3', 'DIRC2', 'KLRC4-KLRK1', 
 'RARRES3', 'FAM173A', 'DOPEY2', 'MUM1', 'HEXDC', 'HIF1A-AS2', 'HKR1', 'ATP5S', 'SSFA2', 
 'ALS2CR12', 'TROVE2', 'SEPT1', 'SEPT10']

#make sure scVI uses the raw counts as input 
rare.layers['normalized'] = rare.X.copy()
rare.X = rare.raw.X

#subset to only the genes found in both
rare_genes_scVI = [gene for gene in list(rare.var_names) if gene not in genes_na]

#read in the query data 
adata = sc.read_10x_h5('/oak/stanford/groups/smontgom/jbezney/scCLEAN/MIRA_topic_models/Control_not_masked_filtered_feature_bc_matrix.h5')
#remove the cells that are not incorporated into the full dataset
annot_ctrl = sc.read_h5ad('/oak/stanford/groups/smontgom/jbezney/scCLEAN/MIRA_topic_models/COntrol_Seurat_annotated.h5ad')
adata = adata[adata.obs_names.isin(list(annot_ctrl.obs.index))]
adata.var_names_make_unique()
adata.obs['batch'] = 'Standard'

#subset to only the genes found in both 
adata = adata[:,rare_genes_scVI]
adata = adata.copy()

#initiate the query model 
model = sca.models.SCVI.load_query_data(
    adata,
    '/oak/stanford/groups/smontgom/jbezney/scCLEAN/scArches/reference_PBMC_scVI_model',
    freeze_dropout = True,
)

model.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0), use_gpu=False)

#save the query model 
model.save('/oak/stanford/groups/smontgom/jbezney/scCLEAN/scArches/query_Standard_PBMC_scVI_model', overwrite=True)
















