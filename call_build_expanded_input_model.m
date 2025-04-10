function [Expanded_Input_model]=call_build_expanded_input_model(generic_input_reconstruction,number_of_cluster,run_toy_model)
% call_build_expanded_input_model: 
%   This function performs the expansion of a consistent model into a multi
%   cell population model. In practice that means that the model is copied
%   as many times as there are different cell populations in the dataset
%   and then they are placed into a shared environment (u compartment) to
%   enable their exchange with each other. 
%
%   INPUT:
%   - generic_input_reconstruction: generic, consistent model - Recon3d model with fastcc run on it 
%   - number_of_clusters: number of different clusters for which you want
%   to expand the model 
%   - run_toy_model: option to run a small toy model - to try stuff
%
%   OUTPUT: 
%   - generic consistent expanded model 
%
%  (c) Leonie Thomas 2025 - University of Luxembourg
%  (c) Maria Pires Pacheco 2022 et al -University of Luxembourg

fields_to_keep={'S','rxns', 'mets', 'metFormulas','genes', 'rules', 'grRules', 'b', 'rxnGeneMat','c', 'lb', 'ub', 'rev', 'subSystems','metFormulas'};
[Expanded_Input_model]= build_expanded_input_model(generic_input_reconstruction, number_of_cluster,run_toy_model,fields_to_keep);
%% sanity check

% creating the model.v slot + turning arround reversable rxns which are
% restricted to carry negative flux by the bounds and have a -1 in the S
% matrix to be reactions which carry positive flux & 1 in the S matrix!
[Expanded_Input_model]=fixIrr_rFASTCORMICS(Expanded_Input_model);
if ~isfield(Expanded_Input_model,'description')
    Expanded_Input_model.description = {'Expanded model'}; % this line should be added into the model builder script
end

%% fix new Exchangers to be secretion only, TS 14.6.24
% put this in a separate function ? 
%%% get uptakes from generic model
mmm=generic_input_reconstruction;
% get exchange rxns abbrevation
temp=findEX_Rxns(mmm,'biomass'); % there is no rxns that is called biomass!! so non of the rxns needs to be excluded + biomass rxns are no exchange rxns...
% get exchange rxns indices
temp=find(ismember(mmm.rxns,temp));
selExc=temp;
% get all the nonzero entries from the exchange rxns colums 
tempM=nonzeros(mmm.S(:,temp));
% selecting for the exchange rxns which are uptakes, two options 
% A -> | -1 in the S matrix | positive bound 
% -> A |  1 in the S matrix | negative bound
% catching both cases in line below 
selUpt=selExc((tempM.*mmm.ub(temp)>0) | (tempM.*mmm.lb(temp)>0));
temp2=selUpt;
% T=table(mmm.rxns(temp2),mmm.lb(temp2),mmm.ub(temp2),printRxnFormula(mmm,mmm.rxns(temp2)))
disp(' ')
disp(['Number of uptakes in generic_input_reconstruction: ' num2str(numel(temp2))])
prevUpt=mmm.rxns(temp2);
prevUpt=strrep(prevUpt,'EX_','Ex_');

%%% get uptakes from expanded model
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
% why do we lose uptake ? have so many new ones ? 
% Number of uptakes in Expanded_Input_model: 1035
% Number of additional uptakes in Expanded_Input_model: 986
% Number of lost uptakes in Expanded_Input_model: 4


disp('... fixing new exchangers ...') % what does fixing mean here? 
% we are closing them ? arent we ? 

% we can get rid of this for loop, vectorize it!
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
% his can also be vectorized 
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
