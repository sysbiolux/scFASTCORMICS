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
run_optimization=0;
printLevel=1;

%Path
changeCobraSolver('ibm_cplex')
user_path='';

% Data_1
biodbnet = readtable([user_path,'\dico_biodbnet.txt']); % The dictionary should map the model genes id to the data ids
load('model.mat')

% uncomment if you get an error concerning eval
% feature astheightlimit 2000

Best_keep=zeros(20,2);
TIME=zeros(20,1);


tic
load('medium_example.mat')
model=removeRxns(model,'biomass_maintenance' );
model=removeRxns(model,'biomass_maintenance_noTrTr' );

biomass_rxn='biomass_reaction'
function_keep='biomass_reaction'
not_medium=''
[model] = constrain_model_rFASTCORMICS(model, medium_example, not_medium, biomass_rxn, function_keep)
A=fastcc_4_rfastcormics(model, 1e-4,1);
model=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A)));

   
[model]=simplyRev_fastbox(model);
generic_input_reconstruction=model;

Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_tumor.txt']);

scdataset = dir([user_path,'\Datasets\Data_', num2str(1),'\*.txt']);
set_name = strcat('Data_',num2str(1),'_model');
if run_optimization==1
    
    coverage = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.05, 0.1, 0.5 0.6];
    REI = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]; % Ex: Percent of 5 = 5 , of 20 = 20
else
    coverage=0.005;
    REI=1;
end

[best,  Final_model, ExpandedInputModel]=scFASTCORMICS(Discretization_Table, set_name, scdataset, biodbnet,user_path, generic_input_reconstruction, coverage,REI, run_optimization, printLevel);

