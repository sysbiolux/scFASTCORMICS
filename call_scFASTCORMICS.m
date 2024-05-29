%% ScFASTCORMICS Driver script
% (c) Maria Pires Pacheco 2022 et al -University of Luxembourg

%% REQUIREMENTS

% FASTCORE, Cobra, and rFASTCORMICS should be added to the path
% FASTCORE and rFASTCORMICS are downloadable on the Sysbiolux/Github
% FASTCORE_1.0 and rFASTCORMICS_2020 scripts
% Last version downloaded: 29/03/2022
% Driver script for the parameter optimization (20 data sets)
% This script is an example of how to run scFASTCORMICS
% The script creates folders to store and organize the data. So, please
% Make sure you have writing rights in the working directory

%% INPUTS: A table contained the genes names and the bulk RNAa-seq values
% A consistent model
% A dictionary between the data IDs and the model IDs
% The single-cell data should be stored in a Folder called Dataset

%run optimization
run_optimization=1;
printLevel=1;

%Path
changeCobraSolver('ibm_cplex')
user_path='C:\Users\maria.pacheco\Documents\GitHub\scFASTCORMICS';

% Data_1
biodbnet = readtable([user_path,'\dico_biodbnet.txt']); % The dictionary should map the model genes id to the data ids
load('simpleRev_recon3.mat','model');
generic_input_reconstruction=model;
clear model
biomass_reaction='biomass_reaction';
load medium_example

% uncomment if you get an error concerning eval
% feature astheightlimit 2000

Best_keep=zeros(20,2);
TIME=zeros(20,1);
function_keep.biomass='biomass_reaction';
medium=Medium.Medium1(~cellfun('isempty',Medium.Medium1));
function_keep.medium=medium;
function_keep.obj='biomass_reaction';
function_keep.not_medium='NONE';

tic

Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_tumor.txt']);

scdataset = dir([user_path,'\Datasets\Data_', num2str(1),'\*.txt']);
set_name = strcat('Data_',num2str(1),'_model_orig');
if run_optimization==1
    
    coverage = [0.0025, 0.005]%, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.05, 0.1, 0.5 0.6];
    REI = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]; % Ex: Percent of 5 = 5 , of 20 = 20
else
    coverage=0.005;
    REI=1;
end
[best, multi_cell_population_model, ExpandedInputModel, A]=scFASTCORMICS(Discretization_Table, set_name, scdataset, biodbnet,user_path, generic_input_reconstruction, coverage,REI, run_optimization, function_keep);


