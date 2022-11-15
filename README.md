# scCLEAN
Single Cell Notebooks for scCLEAN Method Comparison

# Notebook Descriptions: Single Cell Analysis
| File  | Description |
| :-------------: | :-------------: |
| GIT_QC_Batch_Correction_Clustering_Workflow.ipynb  | Quality Control, Sample Aggregation, Batch Correction, Clustering, Trajectory FLE Graph: This processing workflow was performed on all VSMC and VEC cohorts|
| GIT_VECs_CellRank_Trajectory_Inference_Control.ipynb  | Cell Rank Trajectory Analysis performed on the Control VEC cohort: Used to determine the total number of lineages  |
| GIT_VECs_CellRank_Trajectory_Inference_scCLEAN.ipynb  | Cell Rank Trajectory Analysis performed on the scCLEAN VEC cohort: Used to determine the total number of lineages  |
| GIT_VSMC_CellRank_Trajectory_Inference_Control.ipynb  | Cell Rank Trajectory Analysis performed on the Control VSMC cohort: Used to determine the total number of lineages and their association to their original artery identity|
| GIT_VSMC_CellRank_Trajectory_Inference_scCLEAN.ipynb  | Cell Rank Trajectory Analysis performed on the scCLEAN VSMC cohort: Used to determine the total number of lineages and their association to their original artery identity. In addition, this includes the conversion from Cell Rank to scFates in order to identify higher resolution gene to lineage associations|
| GIT_PBMC_DESC_Deep_Learning_Clustering_Variable_Resolution.ipynb  | PBMC single cell clustering improvement validated using an orthogonal deep learning clustering technique (DESC) that tries to reduce user bioinformatician bias|

# Primary Software
- Scanpy
  - https://github.com/scverse/scanpy
- Pegasus 
  - https://github.com/lilab-bcb/pegasus
- CellRank 
  - https://github.com/theislab/cellrank
- scFates 
  - https://github.com/LouisFaure/scFates
- DESC 
  - https://github.com/eleozzr/desc
  
  



