function [Expanded_Input_model]=call_build_expanded_input_model(generic_input_reconstruction,number_of_cluster,run_toy_model)
%  (c) Maria Pires Pacheco 2022 et al -University of Luxembourg

fields_to_keep={'S','rxns', 'mets', 'genes', 'rules', 'grRules', 'b', 'rxnGeneMat','c', 'lb', 'ub', 'rev', 'subSystems'};
[Expanded_Input_model]= build_expanded_input_model(generic_input_reconstruction, number_of_cluster,run_toy_model,fields_to_keep);
%% sanity check
[Expanded_Input_model]=fixIrr_rFASTCORMICS(Expanded_Input_model);
if ~isfield(Expanded_Input_model,'description')
    Expanded_Input_model.description = {'Expanded model'}; % this line should be added into the model builder script
end
end
