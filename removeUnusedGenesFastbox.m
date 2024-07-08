function[model]=removeUnusedGenesFastbox(model, tag)

%(c) Maria Pires Pacheco and Thomas Sauter, 2023 -University of Luxembourg


%The functions removes genes that are not used in the model


% INPUTS
%  model                        COBRA model Structure
%  tag == 1 if the model is a Recon model and the .1 has to be removed

% OUTPUT 
%  model                        COBRA model Structure
genid=1:numel(model.genes);
if tag==1
    model.genes=regexprep(model.genes,'\.[0-9]+$','');
end
if isfield(model, 'rxnGeneMat')
    if numel(model.genes)==size(model.rxnGeneMat,2)
        Unused_genes=( sum(model.rxnGeneMat,1)==0);
        [~, IA, ~] = intersect(genid,genid(~Unused_genes));
        if   isfield(model, 'rxnGeneMat_back')
            Unused_genes_back=( sum(model.rxnGeneMat_back,1)==0);
            Unused_genes=Unused_genes& Unused_genes_back;
            [~, IA, ~] = intersect(genid,genid(~Unused_genes));
        end
        
        if isempty(IA)
            for i=1:numel(IA)
                index=strfind(model.rules,strcat('x(',num2str(IA(i)),')')) ;
                I = model.rules(not(cellfun('isempty',index)));
                model.rules(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(IA(i)),')'), strcat('x(',num2str(i),'*)'));
                if isfield(model.rules,'rules_back')
                    index=strfind(model.rules_back,strcat('x(',num2str(IA(i)),')')) ;
                    I = model.rules_back(not(cellfun('isempty',index)));
                    model.rules_back(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(IA(i)),')'), strcat('x(',num2str(i),'*)'));                    
                end                
            end
            %%            
            model.rules=strrep(model.rules, '*','');
            model.rules=strrep(model.rules, ' ','');
            if isfield(model,'rules_back')
                model.rules_back=strrep(model.rules_back, '*','');
                model.rules_back=strrep(model.rules_back, ' ','');
            end
            model.genes=model.genes(~Unused_genes);
        end
    else
        disp('the number of genes model does not match the rxnGeneMat')
        return
    end
end
model=buildRxnGeneMat(model);
end