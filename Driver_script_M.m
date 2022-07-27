function[best]=Driver_script_M(Discretization_Table, set_name, scdataset,biodbnet, user_path, model, Cover_range,REI_range)

%% ScFASTCORMICS Driver script

%% REQUIREMENTS
tic
% FASTCORE, Cobra,rFASTCORMICS should be added to the path
% FASTCORE and rFASTCORMICS are downloadable on the Sysbiolux/Github
% FASTCORE_1.0 and rFASTCORMICS_2020 scripts
% Last version downloaded: 29/03/2022

% The scripts creates folders to store and organize the data. So please
% make sure you have writing rights in the working directory


% uncomment if you get an error concerning eval 
% feature astheightlimit 2000

%% User input 
%add the path to scFASTCORMICS
run_toy_model = 0;  % 0 if working with dataset, 1 if working with toy data

%load a dictionary 
printLevel=1; %if 1 the script will print all the outputs  if zero the intermediary files not be printed

% choose value (1 or 0) to run plotting scripts (1) or not (0)
plot_formula = 1;


%Create a folder 
path = [user_path,'\1Results_',set_name];
mkdir(path);
files_length = length(scdataset);
if files_length==0
'The single-cell data was not found' 
else



% Names of the cluster for Jaccard (number of names = number of clusters)
Headers = {'Reactions'};
for h=1:files_length
    Headers{h+1} = ['Cluster ',num2str(h-1)];
end


% Building the composite model
number_of_cluster = files_length;


% create a composite input generic model that contains a set of internal
% reactions, metabolites and associated genes peowar cluster in the dataset
[model_composite] = call_build_sc_model_builder(model,number_of_cluster,run_toy_model);
if ~isfield(model_composite,'description')
model_composite.description = {'Recon3D based'}; % this line should be added into the model builder script
end

% Numbering the genes according to their cluster
%Assign a _1,_2,_3 etc to each internal reaction, metabolite,genes to assign it to a cluster 
[input_data]=Numbering(scdataset,path,printLevel);
% 
% 
% % Mapping the cluster genes to the composite model
[input_data]=Mapping(input_data, model_composite, path, printLevel,scdataset);
% save mapping
% 
% 
% % Perform the parameter tuning
%Cover_range = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]/100; % Ex: Cover of 5 = 0.05 , of 20 = 0.2
%REI_range = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]; % Ex: Percent of 5 = 5 , of 20 = 20
Cover_range = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.05, 0.1, 0.5 0.6];
REI_range = [1, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90]; % Ex: Percent of 5 = 5 , of 20 = 20
[input_data,Optimization_global]=Percentile_Cover(model_composite,input_data,Cover_range, REI_range,path,scdataset,printLevel);
% 
toc
T=struct();
% % Extracting additional information from context specific models
[Analysis, T]=Context_Specific_Model_Analysis(model_composite,Optimization_global,files_length, path, Cover_range, REI_range, printLevel, T);
% % Identify reactions that are supported or not supported by the bulk data
[T]=Determine_Reactions_should_or_should_not(model_composite,Discretization_Table, biodbnet, Cover_range, REI_range, Optimization_global, printLevel, T)
% % Identify genes that are supported or not supported by the bulk data
%[T_bulk]=Determine_Genes_should_or_should_not(Discretization_Table,  path, Cover_range, REI_range, printLevel,Metabolic_input_model, Analysis);
% 
% 
% % Determine genes/reactions that should or should not be there in every
% % context specific model
[T]=Intersection_Reactions_Genes(model_composite,T, Cover_range, REI_range, Optimization_global, printLevel);
% 
% % Jaccard similarity plot (see if cluster can be distinguished by their
% % reaction matrix)
% % Jaccard is not interesting for the toy model, should be deactivated
% if printflag == 1
%     mkdir ([extract,'\Jaccard_Similarity']);
% 
% %[c]=Jaccard(Cover_range, REI_range, Analysis, Headers,extract); % is only interesting for real data sets
% end
% 
% 
% % Best threshold determination with formula
[best]=Threshold_Formula(model_composite,path,T, plot_formula,Cover_range, REI_range,Optimization_global);
toc
TIME=toc
save(set_name)
end
end
% 
