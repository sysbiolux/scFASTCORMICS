%% ScFASTCORMICS Driver script
% (c) Maria Pires Pacheco & Thomas Sauter 2022-2024 et al -University of Luxembourg
% Example call for Elena's sc model 2024
%% REQUIREMENTS

% FASTCORE, Cobra, and rFASTCORMICS should be added to the path
% FASTCORE and rFASTCORMICS are downloadable on the Sysbiolux/Github
% FASTCORE_1.0 and rFASTCORMICS_2020 scripts
% Last version downloaded: 29/03/2022
% Driver script for the parameter optimization (20 data sets)
% This script is an example of how to run scFASTCORMICS
% The script creates folders to store and organize the data. So, please
% Make sure you have writing rights in the working directory

clearvars -except solverOK, clc, close all

%% INPUTS: A table contained the genes names and the bulk RNAa-seq values
% A consistent model
% A dictionary between the data IDs and the model IDs
% The single-cell data should be stored in a Folder called Dataset

%run optimization
run_optimization=0;
printLevel=1;

%Path
% % % changeCobraSolver('ibm_cplex')
% % % user_path='C:\Users\maria.pacheco\OneDrive - University of Luxembourg\Documents\GitHub\scFASTCORMICS_o';
user_path='C:\Users\thomas.sauter\OneDrive - University of Luxembourg\work_other\Projects\Elena_2024\scMetMod\TEST_NEW_SCFASTCORMICS_VERSION'

% Data_1
% % % biodbnet = readtable([user_path,'\dico_biodbnet.txt']); % The dictionary should map the model genes id to the data ids
biodbnet = readtable(['dico_biodbnet.txt']); % The dictionary should map the model genes id to the data ids
load('simpleRev_recon3.mat','model');
model.lb(find(ismember(model.rxns,'EX_h2o2[e]')))=0;
model.lb(find(ismember(model.rxns,'EX_o2s[e]')))=0;
model.lb(find(ismember(model.rxns,'EX_oh1[e]')))=0;
model.ub(find(ismember(model.rxns,'EX_oh1[e]')))=0;
model.lb(find(ismember(model.rxns,'EX_ppi[e]')))=0;
model.lb(find(ismember(model.rxns,'sink_fe3[c]')))=0;
model.lb(find(ismember(model.rxns,'sink_band[c]')))=0;
model.ub(find(ismember(model.rxns,'sink_band[c]')))=0;

% biomass_reaction / biomass_maintenance / biomass_maintenance_noTrTr
biomass_reaction='biomass_reaction'; %Recon3D biomass
[NUM,TXT,RAW]=xlsread('medium_withconc.xlsx');
medium=[RAW(:,2); 'hdl_hs[e]'; 'nh4[e]'; 'so3[e]'; 'tsul[e]'; 'tag_hs(e)']
T=table(RAW(:,1),RAW(:,2),NUM,'VariableNames',{'medium_rxns_keep','medium_composition','met_Conc_mM'})
not_medium='NONE' % 'NONE'
unpenalizedSystems = {'Transport, endoplasmic reticular';
    'Transport, extracellular';
    'Transport, golgi apparatus';
    'Transport, mitochondrial';
    'Transport, peroxisomal';
    'Transport, lysosomal';
    'Transport, nuclear'};

% uncomment if you get an error concerning eval
feature astheightlimit 2000
function_keep=struct();
function_keep.biomass=biomass_reaction;
function_keep.medium=medium; %'' %medium
function_keep.obj=biomass_reaction;
function_keep.obj=[function_keep.obj; cellfun(@(c)['EX_' c],medium,'uni',false)]; %keep all medium
function_keep.not_medium=not_medium;
unpenalized = model.rxns(ismember(vertcat(model.subSystems{:}),unpenalizedSystems));
% optional_settings.unpenalized = unpenalized; %don't put transporters to core
function_keep.unpenalized = unpenalized; %don't put transporters to core

mm = constrain_model_rFASTCORMICS(model, function_keep.medium, function_keep.not_medium, function_keep.biomass, function_keep.obj);
A=fastcc_4_rfastcormics(mm, 1e-4,1);
mm=removeRxns(mm, mm.rxns(setdiff(1:numel(mm.rxns),A)));
mm.rxns(find(contains(mm.rxns,'biomass')))

temp=setdiff(function_keep.obj,mm.rxns) %remove unconsistent medium components from obj
function_keep.obj=setdiff(function_keep.obj,temp);

%%
generic_input_reconstruction=mm;

function_keep.medium=''; %'' %medium

Best_keep=zeros(20,2);
TIME=zeros(20,1);
tic

% % % Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_tumor.txt']);
Discretization_Table = [];

% scdataset = dir([user_path,'\Datasets\Data_', num2str(1),'\*.txt']);
scdataset = dir(['ClusterData\*.txt']);
set_name = strcat('Data_',num2str(1),'_model_orig');
if run_optimization==1
    coverage = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.05, 0.1, 0.5 0.6];
    REI = [1, 5];%, 10, 20, 30, 40, 50, 60, 70, 80, 90]; % Ex: Percent of 5 = 5 , of 20 = 20
else
    coverage=0.005;
    REI=1;
end
%index_cells=1:10;% do not forget the gene annotations (3 columns) if you want to run all remouve argument
[best, multi_cell_population_model, ExpandedInputModel, A]=scFASTCORMICS(Discretization_Table, set_name, scdataset,biodbnet, user_path, generic_input_reconstruction, coverage,REI, run_optimization, printLevel, function_keep);
multi_cell_population_model.rxns(find(contains(multi_cell_population_model.rxns,'biomass')))

save WORKSPACE

Results.best=best;
Results.multi_cell_population_model=multi_cell_population_model;
Results.ExpandedInputModel=ExpandedInputModel;
Results.A=A;
nrClusters=length(scdataset)
save(['Results' date],'Results','nrClusters')
