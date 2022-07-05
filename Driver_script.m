%% Driver script
% First add the at least Fastcore and Cobra (contains the rFASTCORMICS) folder to the path
%
% FASTCORE and rFASTCORMICS are downloadable on the Sysbiolux/Github
% FASTCORE_1.0 and rFASTCORMICS_2020 scripts
% Last version downloaded: 29/03/2022
%
% To do list: merge codes together as good as possible
% keeping the load times to a minimum
% Scenario: The user just obtained a clustered data set
% the data set needs to be numbered
% mapped
% discretized
% the program works for 1 data set at a time


% clear; clc;
% tic
% feature astheightlimit 2000

%%%%%%%%%%%%%

% Initialise the cluster data set, names and other variables
% User input here
% user_path = 'D:\ScFASTCORMICS';
% scdataset = dir([user_path,'\Datasets\Data_1\*.txt']); % data set name must be same as in Dataset folder
% set_name = 'Data_1_model'; % name does not matter, do not use confusing names
run_toy_model = 0;  % 0 if working with dataset, 1 if working with toy data
% Discretization_Table = readtable('Discretization_table.txt'); % choose the table corresponding your data set
biodbnet = readtable([user_path,'\dico_biodbnet.txt']); % choose a dictionnary that has right genes (genes of recon3d)
Metabolic_Genes_Recon_2 = readtable([user_path,'\Metabolic_Genes_ID.txt']); % list of metabolic genes of recon3d


% choose value (1 or 0) to run plotting scripts (1) or not (0)
printflag = 1;
plot_formula = 1;

%%%%%%%%%%%%%

files_length = length(scdataset);
path = [user_path,'\1Results_',set_name];
dis = [path,'\Discretization_Step'];
extract = [path,'\Discretization_Step\Extracting_all_possible_information'];
mkdir(path);

% Names of the cluster for Jaccard (number of names = number of clusters)
Headers = {'Reactions'};
for h=1:files_length
    Headers{h+1} = ['Cluster ',num2str(h-1)];
end


% Building the composite model
number_of_cluster = files_length;
%load('Recon204.mat','CmodelR204'); 
load('simpleRev_recon3.mat');
%model = fixIrr(Recon3DModel);
[model_composite] = call_build_sc_model_builder(model,number_of_cluster,run_toy_model);
model_composite.description = {'Recon3D based'}; % this line should be added into the model builder script
%model_composite = fixIrr(model_composite);
save([path,'\composite_model.mat'],'model_composite');
%load([path,'\composite_model.mat']);


% Numbering the genes according to their cluster
Numbering;


% Mapping the cluster genes to the composite model (Recon204)
Mapping;


% Discretization: Percent_Cover
Percentile_Cover;


% Extracting additional information from context specific models
Context_Specific_Model_Analysis;


% Determine reactions that should or should not be in the context and composite models according
% to bulk data
Determine_Reactions_should_or_should_not;


% Determine genes that shoul or should not be in the context and composite models according to
% bulk data
Determine_Genes_should_or_should_not;


% Determine genes/reactions that should or should not be there in every
% context specific model
Intersection_Reactions_Genes;


% Jaccard similarity plot (see if cluster can be distinguished by their
% reaction matrix)
% Jaccard is not interesting for the toy model, should be deactivated
if printflag == 1
Jaccard; % is only interesting for real data sets
end


% Best threshold determination with formula
Threshold_Formula; 

% toc
