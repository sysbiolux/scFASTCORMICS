# Simple running example scFASTCORMICS
## Metabolic network modelling of keratinocyte clonal types

Paper: [Enzo et al, 2021](https://pubmed.ncbi.nlm.nih.gov/33947848/)

Data: [GSE155816](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155816)

In this example we reconstruct a Metabolic Multi-Cell Population Model of timepoint1 using (only) 20 cells per cell-type. This will allow to demonstrate the model building and analysis pipeline with short run time (approx. 20-30 minutes on local desktop). For accurate model reconstruction we suggest to use min 1000 cells per cell type. Thich will increase the run time to approx 1h per cell type, i.e. 7h for this model.

call_scFASTCORMICS_mediumTS: Generates medium constraint Metabolic Multi-Cell Population Model.

driverAnalysisMP_TS: Generates some model stats and performs basic analysis of exchanged metabolites. For this a combined weighted biomass is optimized:

Sum over all clusters i: Cells(i)/TotalCells * biomass(i) -> max
