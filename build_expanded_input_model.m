function[expanded_input_model]= build_expanded_input_model(generic_input_reconstruction, number_of_cluster,run_toy_model, fields_to_keep)
% build_expanded_input_model 
%   DONT USE THIS FUNCTION DIRECTLY IN YOUR SCRIPT, BUT USE THE FUNCTION
%   call_build_expanded_input_model!!!!
%
%
%   This function is called by the call_build_expanded_input_model function. 
%   and returns generic consistent model object - BUT irreversable rxns and uptake
%   rxns still need to be checked before using the expanded model!! 
%       -> see call_build_expanded_input_model
% 
%   INPUT:
%   - generic_input_reconstruction: generic, consistent model - Recon3d model with fastcc run on it 
%   - number_of_clusters:           number of different clusters for which you want
%                                   to expand the model 
%   - run_toy_model:                option to run a small toy model - to try stuff
%   - fields_to_keep:               fields to keep in the model 
%
%   OUTPUT: 
%   - generic consistent model object
%
% (c) Leonie Thomas 2025 - University of Luxembourg
% (c) Maria Pires Pacheco 2022 et al -University of Luxembourg
generic_input_reconstruction.rxns = strrep(generic_input_reconstruction.rxns,'_','');
if run_toy_model==1
    %%
    %     % The toymodel serves as a sanity check
    %     reference_model=createModel;
    %
    %     reference_model=addReaction(reference_model,'R1', '-> A[e]');
    %     reference_model=addReaction(reference_model,'R2', 'A[e] -> B[c]');
    %     reference_model=addReaction(reference_model,'R3', 'A[e] -> C[m]');
    %     reference_model=addReaction(reference_model,'R4','B[c] -> D[e]');
    %     reference_model=addReaction(reference_model,'R5', 'C[m] -> E[c]');
    %     reference_model=addReaction(reference_model,'R6', 'D[e] ->');
    %     reference_model=addReaction(reference_model,'R7', 'E[e] ->');
    %     reference_model=addReaction(reference_model,'R8', 'E[c] <=> E[e]');
    %
    %     reference_model.rev=zeros(numel(reference_model.rxns),1);
    %     reference_model.rev(reference_model.lb<0)=1;
    %     reference_model.rules={'x(1)';'x(2)';'x(3)';'x(4)';'x(5)';'x(6)';'x(7)';'x(8)'};
    %     reference_model.genes={'8001';'8002';'8003';'8004';'8005';'8006';'8007';'8008'};
    %     reference_model = buildRxnGeneMat(reference_model);
    %model=reference_model;
    %     Original toy model
    reference_model=createModel;
    
    reference_model=addReaction(reference_model,'R1', '-> A[e]');
    reference_model=addReaction(reference_model,'R2', 'A[c] + B[c] -> C[c] + 1D[m]');
    reference_model=addReaction(reference_model,'R3', ' C[g] + 1E[g] ->  G[g]');
    reference_model=addReaction(reference_model,'R4', ' 1D[m] + H[m]  -> I[m] + J[m] ');
    reference_model=addReaction(reference_model,'R5', ' <=> 1E[e] ');
    reference_model=addReaction(reference_model,'R6', '<=> H[e] ');
    reference_model=addReaction(reference_model,'R7', '<=> G[e] ');
    reference_model=addReaction(reference_model,'R8', '<=> J[e]');
    reference_model=addReaction(reference_model,'R9', '-> B[n] ');
    reference_model=addReaction(reference_model,'R10', '<=> I[e]');
    reference_model=addReaction(reference_model,'R11', 'I[m] + K[e] <=> I[e] + K[m]');
    reference_model=addReaction(reference_model,'R12', 'B[c] + K[m] <=> B[n] + K[e] ');
    reference_model=addReaction(reference_model,'R13', 'J[m] <=> J[e]');
    reference_model=addReaction(reference_model,'R14', 'G[g] <=> G[e] ');
    reference_model=addReaction(reference_model,'R15', 'H[m] <=> H[e]');
    reference_model=addReaction(reference_model,'R16', 'K[m] <=> ');
    reference_model=addReaction(reference_model,'R17', '1E[g] <=> 1E[e] ');
    reference_model=addReaction(reference_model,'R18', 'A[c] <=> A[e] ');
    reference_model=addReaction(reference_model,'R19', 'C[c] <=> C[g] ');
    reference_model=addReaction(reference_model,'R20', 'K[e] <=> ');
    
    reference_model.rev=zeros(numel(reference_model.rxns),1);
    reference_model.rev(reference_model.lb<0)=1;
    reference_model.rules={'x(1)';'x(2)';'x(3)';'x(4)';'x(5)';'x(6)';'x(7)';'x(8)';'x(9)';'x(10)';'x(11)';'x(12)';'x(13)';'x(14)';'x(15)';'x(16)';'x(17)';'x(18)'; 'x(19)';'x(20)'};
    reference_model.genes={'GA';'GB';'GC';'GD';'GE';'GF';'GG';'GH';'GI';'GJ';'GK';'GL';'GM';'GN';'GO';'GP';'GQ';'GR';'GS'; 'GT'};
    reference_model = buildRxnGeneMat(reference_model);
    model=reference_model;
    
else
    model=generic_input_reconstruction;
    %% check consistency of the generic model
    model=fixIrr_rFASTCORMICS(model);
    A=fastcc_fastcore(model, 1e-4);
    model=removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A)));
    
end
model.genes=regexprep(model.genes,'\.[0-9]+$','');
% find exchange reactions(exRxns),exchange metabolites(ex_mets),
exRxns=find(sum(abs(model.S),1)==1);
[exmets,~]=find(model.S(:, exRxns));
exmets=exmets(contains(model.mets(exmets),'[e]'));
exmets=unique(model.mets(exmets));
[~,exmetsID,~]=intersect(model.mets, exmets);
[~,r]=find(model.S(exmetsID,:));
exRxns=intersect(exRxns,r);
exRxns=model.rxns(exRxns);

% find metabolites that will be in the [u]
u_mets=unique(regexprep(exmets,'\[e]','[u]'));
% find internal reactions(internal_rxns), internal metabolites (internal_mets),
internal_mets=setdiff(model.mets,exmets);
[~,internal_metsID,~]=intersect(model.mets, internal_mets);
[~,internal_rxnsID]=find(model.S(internal_metsID, :));
internal_rxns=model.rxns(internal_rxnsID);
internal_rxns=unique(internal_rxns);
[~,intercluster_transport_ID]=find(model.S(exmetsID, :));
intercluster_transport=model.rxns(intercluster_transport_ID);
internal_rxns=setdiff(internal_rxns,intercluster_transport);
% find intercluster_transport(internal_rxns), and metabolites that are transported between the clusters and the intersticial  (internal_2_u_mets),
intercluster_transport=setdiff(intercluster_transport,exRxns);
[~,intercluster_transport_ID,~]=intersect(model.rxns, intercluster_transport);
[internal_2_u_mets_ID,~]=find(model.S(:, intercluster_transport_ID));
internal_2_u_mets=model.mets(internal_2_u_mets_ID);
internal_2_u_mets=intersect(internal_2_u_mets,internal_mets);

%% create a empty expanded_input_model
[expanded_input_model, fields]=create_empty_expanded_input_model(model, fields_to_keep);
%% multiply the internal metabolites and intercluster_transport

fields=setdiff(fields,'genes');
fields=setdiff(fields,'rxnGeneMat');
fields=setdiff(fields,'grRules');
fields=setdiff(fields,'mets');
fields=setdiff(fields,'rxns');
fields=setdiff(fields,'S');

[expanded_input_model]=compute_expanded_input_model(model,expanded_input_model,internal_mets,internal_rxns,intercluster_transport,internal_2_u_mets,u_mets,exmets,exRxns,number_of_cluster,fields);

end
