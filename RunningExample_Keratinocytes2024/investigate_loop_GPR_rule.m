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



