function [model_composite]=call_build_sc_model_builder(input_model,number_of_cluster,run_toy_model)

% add Fastcore and Cobra to the path
% script that calls the build_scRNAseq_model

fields_to_keep={'S','rxns', 'mets', 'genes', 'rules', 'grRules', 'b', 'rxnGeneMat','c', 'lb', 'ub', 'rev', 'subSystems'};
[model_composite]= build_scRNAseq_model(input_model, number_of_cluster,run_toy_model,fields_to_keep);
%% sanity check
[model_composite]=fixIrr(model_composite);
% A=fastcc(model_composite, 1e-4);
% if numel(A)< numel(model_composite.rxns)
%     'composite_model is inconsistent';
%     cvvvvvvvvvvvvvvvvvvx
end
%T=table(model_composite.rxns,printRxnFormula(model_composite, model_composite.rxns));