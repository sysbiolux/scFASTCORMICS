% this script is to build a minimal toy model which helps me to understand
% the build expanded model function from maria 
% this model is adapted from the build_expanded_input_model.m script
% within the scFastcormics package, a few reactions are taken out to make
% the toy model smaller and therefore easier to understand

% his should help to get started with the spatial model, a small grid model
if small_model
    reference_model=createModel;
    reference_model=addReaction(reference_model,'R4', ' 1D[m] + H[m]  -> I[m] + J[m] ');
    reference_model=addReaction(reference_model,'R6', '<=> H[e] ');
    reference_model=addReaction(reference_model,'R8', '<=> J[e]');
    reference_model=addReaction(reference_model,'R9', '-> B[n] ');
    reference_model=addReaction(reference_model,'R10', '<=> I[e]');
    reference_model=addReaction(reference_model,'R11', 'I[m] + K[e] <=> I[e] + K[m]');
    reference_model=addReaction(reference_model,'R12', 'B[c] + K[m] <=> B[n] + K[e] ');
    reference_model=addReaction(reference_model,'R13', 'J[m] <=> J[e]');
    reference_model=addReaction(reference_model,'R15', 'H[m] <=> H[e]');
    reference_model=addReaction(reference_model,'R16', 'K[m] <=> ');
    reference_model.rev=zeros(numel(reference_model.rxns),1);
    reference_model.rev(reference_model.lb<0)=1;
    reference_model.rules={'x(1)';'x(2)';'x(3)';'x(4)';'x(5)';'x(6)';'x(7)';'x(8)';'x(9)';'x(10) | x(11)'};%;'x(12)';'x(13)';'x(14)';'x(15)'};% ;'x(16)';'x(17)';'x(18)'; 'x(19)';'x(20)'};
    reference_model.genes={'GA';'GB';'GC';'GD';'GE';'GF';'GG';'GH';'GI';'GJ';'GK'};%;'GL';'GM';'GN';'GO'};%;'GP';'GQ';'GR';'GS'; 'GT'};
    reference_model = buildRxnGeneMat(reference_model);
    toy_model=reference_model;
    model = toy_model;

    graphObj=createMetIntrcNetwork(model,model.mets);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does

    clear toy_model reference_model

else 
    load('simpleRev_recon3.mat','model');
end