function[model]=removeUnusedGenesFastbox(model, tag)
% removeUnusedGenesFastbox Removes genes not utilized in the model.
%
%   model = removeUnusedGenesFastbox(model, tag) filters out genes from the
%   COBRA model that are not used in the model's reactions.
%   If 'tag' is set to 1, the function also removes postfix ('\.[0-9]+$')
%   from gene names (in case of Recon it removes ".1").
%
% INPUTS
%  model         COBRA model structure
%  tag           1 to remove postfix from the genenames, 0 if nothing
%                should be removed
%                ('\.[0-9]+$', for recon ".1")
%
% OUTPUT 
%  model         Updated COBRA model structure with unused genes removed
%                from the following fields:
%                - 'genes'
%                - 'rxnGeneMat'
%                - 'rules'
%                - 'rxnGeneMat_back' (if present)
%                - 'rules_back' (if present)
%
% NOTE: There is also a removeUnusedGenes function in the cobra toolbox. 
% The cobra toolbox function is slower because it executes the 
% removeFieldEntriesForType function. In contrast to the function implemented
% in the cobratoolbox, removeUnusedGenesFastbox defined in scFASTCORMICS 
% only goes through genes,rxnGeneMat,rxnGeneMat_back, rules, and rules_back slot
% of the model, not through all slots of the model! 
%
%(c) Leonie Thomas, 2025 - University of Luxembourg
%(c) Maria Pires Pacheco and Thomas Sauter, 2023 -University of Luxembourg

genid=1:numel(model.genes);
if tag==1
    model.genes=regexprep(model.genes,'\.[0-9]+$','');
    disp("Remove postfix after the period for every gene!!")
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
        
        if ~isempty(IA)
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