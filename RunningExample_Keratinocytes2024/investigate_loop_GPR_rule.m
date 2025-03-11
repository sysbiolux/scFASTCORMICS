load workspace_loop_GPRrulesMapper_rFASTCORMICS


mappinga = mapping; 
mappingb = mapping;

numel(model_composite.rxns(match))

tic
for ki=1:numel(model_composite.rxns(match))
        for j=1:20
            mappinga(match(ki),j)= GPRrulesMapper_rFASTCORMICS(cell2mat(model_composite.rules(match(ki))),...
                                                              data(:,j));
        end
end
toc

tic
rules = regexprep(model_composite.rules,'x\(([0-9]*)\)','x($1,:)');
for ki=1:numel(model_composite.rxns(match))
            mappingb(match(ki),:)= GPRrulesMapper_rFASTCORMICS(cell2mat(rules(match(ki))),...
                                                          data(:,1:20));
end
toc

% are the two mappings the same ? 
sum(mappinga ~= mappingb)


%% toy data - to check if the matrix wise computation of the GPR are correct
% - there are online posts saying that a eval can not be used in a parfor
% loop - https://nl.mathworks.com/matlabcentral/answers/579876-how-can-we-handle-the-conditional-use-of-a-parfor-loop-with-the-eval-function
% - the weird thing is that it works for my toy example but not for the big
% one 
% error message: Error using GPRrulesMapper_rFASTCORMICS (line 4) Error: The input was too complicated or too big for MATLAB to parse.

addpath(genpath("path to fastcormics to be able to runGPRulesMapper"))
data = [ 1 0 5 8;
         0 3 0 8;
         0 4 6 8;
         1 5 2 7; 
         1 6 3 1];

correct_results = [ 0 3 5 8;
                    1 4 2 1;
                    1 4 6 8];

rules = ["( x(1) | x(2) ) & x(3)";...
         "( x(1) | x(2) | x(3) ) & x(4) & x(5)";...
         "( x(4) & x(1) ) | ( x(2) | x(3) )" ]

disp("changed rules:")
rules = regexprep(rules,'x\(([0-9]*)\)','x($1,:)')'
mapping = zeros(numel(rules), size(data,2));

tic
for k=1:numel(rules)
            mapping(k,:)= GPRrulesMapper_rFASTCORMICS(cell2mat(rules(k)),...
                                                          data);
end
toc


mapping = zeros(numel(rules), size(data,2));
tic
parfor k=1:numel(rules)
            mapping(k,:)= GPRrulesMapper_rFASTCORMICS(cell2mat(rules(k)),...
                                                          data);
end
toc

sum(mapping ~= correct_results)


%% this is a script to investigate how to perform parallel computing in
% matlab 

ver % displays which toolboxes are installed 


parfor i = 1:100
    c(i) = max(eig(rand(1000)));
end









