# scFASTCORMICS

Getting Started:
scFASTCORMICS builds multi-cell population models using single-cell data as input and a generic input reconstruction (a COBRA model) such as Recon3D, or human1 for humans and optionally bulk RNAseq data from the same context to perform parameter optimization. 
Call_scFASTCORMICS gives an example on how to call the tool. Modifications such as path changing, changes in the folder name, imports of the inputs should only be done in the call_scFASTCORMICS scripts. 

Requirements:
scFASTCORMICS requires the Cplex solver, the COBRA toolbox and rFASTCORMICS that can be downloaded https://www.ibm.com/analytics/cplex-optimizer, from https://opencobra.github.io/cobratoolbox/stable/ and  https://github.com/sysbiolux/rFASTCORMICS, respectively. The solver and both tools must be in the path. 

Inputs:
The single-cell data should be stored in a Folder called dataset with subfolder with the called Data1, Data2, etc. if the building is performed for more dataset. 
The Data folders should contain text files, one per cluster, in table form with the row names being the genes, the columns corresponding to the cells in the cluster and the entries to the normalized intensities obtained after running Seurat or any other clustering tool.

The bulk data should be first discretized the map_2_expression function of rFASTCORMICS. To obtain a table with the data identifier as row names and the sample as columns. The entries are discretized values: 1 (expressed), (0) unknown, -1 (not expressed) 

A dictionary should also be provided that match the data identifier to the model identifiers. 

Running Times:
The model building can be run with the default parameters that were optimized for Recon3D or with tailored parameter for the dataset and model. The running time without parameter optimization is around 10 minutes per model while with the parameter optimization step, the running time can be from a few hours for smaller dataset with few clusters (4-5) to several days for the dataset with 15 and more clusters.

Other Remark(s):
scFASTCORMICS will save heatmaps and files on your computer so make sure that you have writing rights.


