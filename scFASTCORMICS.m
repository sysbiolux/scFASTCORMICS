function[best,  multicell_model,Expanded_input_model ]=scFASTCORMICS(Discretization_Table, set_name, scdataset,biodbnet, user_path, generic_input_reconstruction, Cover_range,REI_range, run_opti, printLevel)
% (c) Maria Pires Pacheco 2022 et al -University of Luxembourg

best=[];
plot_formula = 1;
run_toy_model=0;
%Create a folder for the results
path = [user_path,'\1Results_',set_name];
mkdir(path);
files_length = length(scdataset);
if files_length==0
    % checks the inputs'
    'There is no single-cell data';
    return
else    
    % Building the Expanded Input model
    number_of_cluster = files_length;
    
    % create an expanded input model that contains a set of internal
    % reactions, metabolites and associated genes per cluster in the dataset
    [Expanded_input_model] = call_build_expanded_input_model(generic_input_reconstruction,number_of_cluster,run_toy_model);
 
    
    % Numbering the genes according to their cluster
    % Assign a _1,_2,_3 etc to the identifier of each internal reaction, metabolite, and gene to assign it to a cluster
    [input_data]=Numbering(scdataset,path,printLevel);
    
    % Mapping the cluster genes to the composite model
    [input_data]=Mapping(input_data, Expanded_input_model, path, printLevel,scdataset);
    
    % Perform the parameter tuning
    if run_opti==1
        [~,Optimization_global]=Parameter_optimisation(Expanded_input_model,input_data,Cover_range, REI_range,path,scdataset,printLevel);
        T=struct();
        % Extracting additional information from the expanded input model
        [~, T]=Multi_cell_population_model_Analysis(Expanded_input_model,Optimization_global,files_length, path, Cover_range, REI_range, printLevel, T);
        %  Identify reactions that are supported or not supported by the
        %  bulk data (BE and BU)
        [T]=Determine_BE_BU(Expanded_input_model,Discretization_Table, biodbnet, Cover_range, REI_range, Optimization_global, T);
           
        [T]=Intersection_Reactions_Genes(Expanded_input_model,T, Cover_range, REI_range, Optimization_global);
        
        %Determine optimal REI and cover
        [best]=Determine_optimal_REI_and_coverage(Expanded_input_model,path,T, plot_formula,Cover_range, REI_range,Optimization_global);
        toc
        TIME=toc;
    end
    if run_opti==0
        best.Best_Cover_Threshold_Reaction_Formula=Cover_range;
        best.Best_REI_Threshold_Reaction_Formula=REI_range;
    end
    
   [~,Final_model]=Build_Multi_cell_population_model(Expanded_input_model,input_data, best.Best_Cover_Threshold_Reaction_Formula,best.Best_REI_Threshold_Reaction_Formula,path,scdataset,printLevel);
    multicell_model = removeRxns(Expanded_input_model,Expanded_input_model.rxns(setdiff(1:numel(Expanded_input_model.rxns),find(Final_model.A))));
end
save(set_name)

end
%
