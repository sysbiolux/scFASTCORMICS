function [Expanded_Input_model]=call_build_expanded_input_model(generic_input_reconstruction,number_of_cluster,run_toy_model)
% call_build_expanded_input_model: 
%   This function performs the expansion of a consistent model to a multi
%   cell population model. The model is copied
%   as many times as there are different cell populations in the dataset
%   and then these are placed into a shared environment (u compartment) to
%   enable their exchange with each other. 
%
%   INPUT:
%   - generic_input_reconstruction: generic, consistent model - e.g. Recon3D 
%   - number_of_clusters: number of different clusters for which model shall be expanded
%   - run_toy_model: option to run a small toy model - to try stuff
%
%   OUTPUT: 
%   - generic consistent expanded model 
%
%  (c) Maria Pires Pacheco 2022 et al -University of Luxembourg

fields_to_keep={'S','rxns', 'mets', 'metFormulas','genes', 'rules', 'grRules', 'b', 'rxnGeneMat','c', 'lb', 'ub', 'rev', 'subSystems','metFormulas'};
[Expanded_Input_model]= build_expanded_input_model(generic_input_reconstruction, number_of_cluster,run_toy_model,fields_to_keep);
%% sanity check
[Expanded_Input_model]=fixIrr_rFASTCORMICS(Expanded_Input_model);
if ~isfield(Expanded_Input_model,'description')
    Expanded_Input_model.description = {'Expanded model'}; % this line should be added into the model builder script
end

%% fix new Exchangers to be secretion only, TS 14.6.24

mmm=generic_input_reconstruction;
temp=findEX_Rxns(mmm,'biomass');
temp=find(ismember(mmm.rxns,temp));
selExc=temp;
tempM=nonzeros(mmm.S(:,temp));
selUpt=selExc((tempM.*mmm.ub(temp)>0) | (tempM.*mmm.lb(temp)>0));
temp2=selUpt;
% T=table(mmm.rxns(temp2),mmm.lb(temp2),mmm.ub(temp2),printRxnFormula(mmm,mmm.rxns(temp2)))
disp(' ')
disp(['Number of uptakes in generic_input_reconstruction: ' num2str(numel(temp2))])
prevUpt=mmm.rxns(temp2);
prevUpt=strrep(prevUpt,'EX_','Ex_');

mmm=Expanded_Input_model;
temp=findEX_Rxns(mmm,'biomass');
temp=find(ismember(mmm.rxns,temp));
selExc=temp;
tempM=nonzeros(mmm.S(:,temp));
selUpt=selExc((tempM.*mmm.ub(temp)>0) | (tempM.*mmm.lb(temp)>0));
temp2=selUpt;
% T=table(mmm.rxns(temp2),mmm.lb(temp2),mmm.ub(temp2),printRxnFormula(mmm,mmm.rxns(temp2)))
disp(['Number of uptakes in Expanded_Input_model: ' num2str(numel(temp2))])
expUpt=mmm.rxns(temp2);
newUpt=setdiff(expUpt,prevUpt);
disp(['Number of additional uptakes in Expanded_Input_model: ' num2str(numel(newUpt))])
lostUpt=setdiff(prevUpt,expUpt);
disp(['Number of lost uptakes in Expanded_Input_model: ' num2str(numel(lostUpt))])

disp('... fixing new exchangers ...')
for counter=1:numel(newUpt)
    temp=find(ismember(mmm.rxns,newUpt(counter)));
    if mmm.lb(temp)==0 & mmm.ub(temp)>0 & nonzeros(mmm.S(:,temp))>0
        mmm.lb(temp)=-mmm.ub(temp);
        mmm.ub(temp)=0;
    end
%     if mmm.lb(temp)<0 & mmm.ub(temp)==0 & nonzeros(mmm.S(:,temp))<0
%         mmm.ub(temp)=-mmm.lb(temp);
%         mmm.lb(temp)=0;
%     end
end
% T=table(mmm.rxns(temp2),mmm.lb(temp2),mmm.ub(temp2),printRxnFormula(mmm,mmm.rxns(temp2)))

% adding lost uptakes back
disp('... adding lost uptakes back ...')
for counter=1:numel(lostUpt)
    temp=find(ismember(mmm.rxns,lostUpt(counter)));
    mmm.lb(temp)=-1000;
end

temp=findEX_Rxns(mmm,'biomass');
temp=find(ismember(mmm.rxns,temp));
selExc=temp;
tempM=nonzeros(mmm.S(:,temp));
selUpt=selExc((tempM.*mmm.ub(temp)>0) | (tempM.*mmm.lb(temp)>0));
temp2=selUpt;
% T=table(mmm.rxns(temp2),mmm.lb(temp2),mmm.ub(temp2),printRxnFormula(mmm,mmm.rxns(temp2)))
disp(['Number of uptakes in Expanded_Input_model: ' num2str(numel(temp2))])
expUpt=mmm.rxns(temp2);
newUpt=setdiff(expUpt,prevUpt);
disp(['Number of additional uptakes in Expanded_Input_model: ' num2str(numel(newUpt))])
lostUpt=setdiff(prevUpt,expUpt);
disp(['Number of lost uptakes in Expanded_Input_model: ' num2str(numel(lostUpt))])

Expanded_Input_model=mmm;

end
