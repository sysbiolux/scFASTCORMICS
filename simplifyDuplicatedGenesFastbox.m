function[model, changed_rules_ON]=simplifyDuplicatedGenesFastbox(model, taggene,display_changed_mixed_rules)
% simplifyDuplicatedGenesFastbox: Removes genes not utilized in the model &
% gets rid of duplicated genes in the model. 
%
% In the first step, it uses the rxnGeneMat slot of the model to remove
% genes with no associated rxns in the model -> see removeUnusedGenesFastbox function
% 
% In a second step, it checks for duplicated genes in the model. 
% Since the indices of genes in the .gene slot change by deleting duplicates,
% the rules & rxnGeneMat slot need to be updated accordingly, so their entries
% index the correct genes.
%
% INPUTS
%  model                                COBRA model Structure
%  taggene                              1 to remove postfix from the genenames
%                                       ('\.[0-9]+$', for recon ".1")
%  display_changed_mixed_rules          1 to see which MIXED rules were
%                                       changed (default = 0)
% 
% OUTPUT 
%  model            Updated COBRA model structure with unused/duplicated genes 
%                   removed from the following fields:
%                   - 'genes'jd
%                   - 'rxnGeneMat'
%                   - 'rxnGeneMat_back' (if present)
%                   ... and ...
%                   - 'rules'
%                   - 'rules_back' (if present)
%                   updated & simplified accordingly.
%
% changed_rules_ON  Vector specifying which of the MIXED rules where changed by the
%                   function (1 changed, 0 unchanged).
%                   BUT: This is only the mixed rules! not the ones which
%                   exclusively entail & or | !!!  This vector can be used to
%                   understand how the more complex rules are simplified
%                   using this function.
%
%(c) Leonie Thomas, 2025 - University of Luxembourg
%(c) Maria Pires Pacheco and Thomas Sauter, 2023 -University of Luxembourg

arguments
    model (1,1) struct 
    taggene (1,1) double =0
    display_changed_mixed_rules (1,1) double =0
end

%%%%%%%%%%%%%%%%%%% Removing unused genes from the Model

% checking that there is a rxnGeneMat slot -> used to determine unused Genes
if isfield(model, 'rxnGeneMat')
    if numel(model.genes)~=size(model.rxnGeneMat,2)
        disp('simplifyDuplicatedGenes_fastbox was not run') 
        return
    end
end

[model]=removeUnusedGenesFastbox(model,taggene);

%%%%%%%%%%%%%%%%%%% Removing duplicated genes from the Model

[genes,~,~] = unique(model.genes);

%----- replace all duplicated genes in the rules --------

% - unique genes are selected
% - since we loose some gene ids, the ids the rules slot refers to needs to
%   be updated -> enter a for loop over all genes 
% - for non unique genes multiple all indices of all copies needs to be
%   replaced 
% - in the rules and in the rules_back object!

if numel(genes) < numel(model.genes)
    %find non unique genes
    [genes,ia,ic] = unique(model.genes);
    for i=1:numel(ia)
        % looping over all unique genes 
        dupidx = find(ic == i);
        % find indices for the first gene in model.genes
        if numel(dupidx)>1  % is it duplicated ?       
            first=min(dupidx); % take the first gene
            index=strfind(model.rules,strcat('x(',num2str(first),')'));
            I = model.rules(not(cellfun('isempty',index)));
            if~isempty(I)
                model.rules(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(first),')'), strcat('x(',num2str(i),'*)'));
            else
                
            end
            if isfield(model,'rules_back')
                I = model.rules_back(not(cellfun('isempty',index)));% % %
                if~isempty(I)
                    model.rules_back(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(first),')'), strcat('x(',num2str(i),'*)'));
                else
                    
                end
            end
            dupidx=setdiff(dupidx,first);
            for ii=1:numel(dupidx)
                index=strfind(model.rules,strcat('x(',num2str(dupidx(ii)),')'));
                I = model.rules(not(cellfun('isempty',index)));
                
                if ~isempty(I)
                    model.rules(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(dupidx(ii)),')'), strcat('x(',num2str(i),'*)'));
                else
                    
                end
                if isfield(model,'rules_back')
                    I = model.rules_back(not(cellfun('isempty',index)));%
                    
                    if ~isempty(I)
                        model.rules_back(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(dupidx(ii)),')'), strcat('x(',num2str(i),'*)'));
                    else
                        
                    end
                end
            end
        else
            
            
            first=min(dupidx);
            index=strfind(model.rules,strcat('x(',num2str(first),')'));
            I = model.rules(not(cellfun('isempty',index)));
            
            if  ~isempty(I)                
                model.rules(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(first),')'), strcat('x(',num2str(i),'*)'));
            else
                                
            end
            
            if isfield(model,'rules_back')
                I = model.rules_back(not(cellfun('isempty',index)));% % %
                
                if  ~isempty(I)                    
                    model.rules_back(not(cellfun('isempty',index)))=strrep(I, strcat('x(',num2str(first),')'), strcat('x(',num2str(i),'*)'));
                else                    
                    
                end
            end
            
        end
    end
    
    model.rules=strrep(model.rules, '*','');
    model.rules=strrep(model.rules, ' ','');
    model.genes=genes;
    if isfield(model,'rules_back')
        model.rules_back=strrep(model.rules_back, '*','');
        model.rules_back=strrep(model.rules_back, ' ','');
    end
    
    
    % ----- replace all duplicated genes in the rules --------
    
    % some of the rules will now entail multiple entries for the same gene id,
    % since we replaced all the duplicated genes with the same number, 
    % next step is to get rid of these duplicates
    
    
    %%first we handle all the cases where there is only & in the fomula
    Index=strfind(model.rules, '&') ;
    match=(cellfun('isempty',Index) & ~strcmp(model.rules, '')); % find statement without or
    In=strfind( model.rules, '|');
    match2=(not(cellfun('isempty',In)));
    match=match & match2;
    A=model.rules(match);
    newA = cellfun(@(s) strsplit(s, '|'), A, 'UniformOutput', false); %or use regexp
    new2=cell(size(newA,1),1);
    for i=1:numel(newA)
        %Cu = cellfun(@(y) unique(y{:} ), newA, 'UniformOutput',false)
        C=unique(newA{i});
        AllC = {cat(2, C{:})};
        AllC= strrep(AllC, ')x', ') | x');
        AllC= strrep(AllC, ')(x', ') | (x');       
        new2(i)=AllC;
    end
    
    model.rules(match)=new2;
        
    %%now we handle the formulas with only | in them   
    Index=strfind(model.rules, '&') ;
    
    match=not(cellfun('isempty',Index));
    In=strfind( model.rules, '|');
    match2=((cellfun('isempty',In)));
    match=match & match2;
    
    
    A=model.rules(match);
    newA = cellfun(@(s) strsplit(s, '|'), A, 'UniformOutput', false); %or use regexp
    new2=cell(size(newA,1),1);
    for i=1:numel(newA)
        C=unique(newA{i});
        AllC = {cat(2, C{:})};
        AllC= strrep(AllC, ')x', ') & x');
        AllC= strrep(AllC, ')(x', ') & (x');        
        new2(i)=AllC;
    end
    model.rules(match)=new2;
    
    %now do the same for rules_back!
    
    if isfield(model, 'rules_back')
        Index=strfind(model.rules_back, '&') ;
        match=(cellfun('isempty',Index) & ~strcmp(model.rules, '')); % find statement only with OR
        In=strfind( model.rules_back, '|');
        match2=(not(cellfun('isempty',In)));
        match=match & match2;
        A=model.rules_back(match);
        newA = cellfun(@(s) strsplit(s, '|'), A, 'UniformOutput', false); %or use regexp
        new2=cell(size(newA,1),1);
        for i=1:numel(newA)
            %Cu = cellfun(@(y) unique(y{:} ), newA, 'UniformOutput',false)
            C=unique(newA{i});
            AllC = {cat(2, C{:})};
            AllC= strrep(AllC, ')x', ') | x');
            AllC= strrep(AllC, ')(x', ') | (x');         
            new2(i)=AllC;
        end
        
        model.rules_back(match)=new2;
        Index=strfind(model.rules, '&') ;        %
        match=not(cellfun('isempty',Index));
        In=strfind( model.rules_back, '|');
        match2=((cellfun('isempty',In)));
        match=match & match2;
        
        
        A=model.rules_back(match);
        newA = cellfun(@(s) strsplit(s, '|'), A, 'UniformOutput', false); %or use regexp
        new2=cell(size(newA,1),1);
        for i=1:numel(newA)
            C=unique(newA{i});
            AllC = {cat(2, C{:})};
            AllC= strrep(AllC, ')x', ') & x');
            AllC= strrep(AllC, ')(x', ') & (x');
            
            new2(i)=AllC;
        end
    end%
end


model=buildRxnGeneMat(model);
mixed_rules=find(contains(model.rules, '&') & contains(model.rules, '|'));
OR_rules=find(~contains(model.rules, '&') & contains(model.rules, '|'));
AND_rules=find(contains(model.rules, '&') & ~contains(model.rules, '|'));

% after correcting the rules with only & and the rules with only | in them,
% there are still the mixed rules which have both and and or
% looping over all the mixed rules!


changed_rules_ON=zeros(numel(model.genes),1);
model_old_rules = model.rules;

for i=1:numel(mixed_rules)
    
    [~,g]=find(model.rxnGeneMat(mixed_rules(i),:));
    OrGenes=zeros(numel(g),1);
    nb_AND = count(model.rules(mixed_rules(i)),'&');
    nb_OR= count(model.rules(mixed_rules(i)),'|');
    
    if  numel(g)< (nb_OR + nb_AND +1)
        for ii=1:numel(g)
            x=zeros(numel(model.genes),1);
            xm=x;
            x(g(ii))=1;
            xm(g(ii))=-1;
            rr=mixed_rules(i);
            mapping=zeros(numel(model.rxns),1);
            
            % this for loop finds out, which of the genes are part of a & or an |
            % operation by setting the expression of this gene to 1/-1 and the rest to
            % 0, by running the GPRrulemapper, you get the information wether the rxns
            % is active or not -> if the gene is an or gene then the -1/1 will not
            % appear in the mapping of the rxns, since the other gene in the rules
            % which is connected by a or is set to 0
            
            for k=1:numel(model.rxns(rr))
                mapping(rr(k),1)= GPRrulesMapper_rFASTCORMICS(cell2mat(model.rules(rr(k))),x);
            end
            mapping_m=zeros(numel(model.rxns),1);
            
            for k=1:numel(model.rxns(rr))
                mapping_m(rr(k),1)= GPRrulesMapper_rFASTCORMICS(cell2mat(model.rules(rr(k))),xm);
            end
            if sum(mapping(rr))==1 && sum(mapping_m(rr)==-1)==0
                
                OrGenes(ii)=1;
            end
            
            if sum(mapping(rr))==1 && sum(mapping_m(rr)==-1)==1
                OrGenes(ii)=2;
            end
        end
        
        AndGenes=~OrGenes;
        AndGenes2=zeros(size(AndGenes));
        x=zeros(numel(model.genes),1);
        x(g(AndGenes))=1;
        xm=x;
        
        % now after defining all and and or genes 
        % the rules are assembled back together to form shorter rules
        % without duplicates
        AndGenes=find(AndGenes);
        for iii=1:numel(AndGenes)
            xm(g(AndGenes(iii)))=-1;
            mapping=zeros(numel(model.rxns),1);
            for k=1:numel(model.rxns(rr))
                mapping(rr(k),1)= GPRrulesMapper_rFASTCORMICS(cell2mat(model.rules(rr(k))),x);
            end
            for k=1:numel(model.rxns(rr))
                mapping_m(rr(k),1)= GPRrulesMapper_rFASTCORMICS(cell2mat(model.rules(rr(k))),xm);
            end
            if sum(mapping(rr))==1 && sum(mapping_m(rr) ~=0)==0
                AndGenes2(iii)=1;
            end
        end
        if sum(AndGenes2~=OrGenes)== numel(g)
            
            % to make it more clear, examples from the Keratinocyte Example are put into the comments
            gAnd=g(AndGenes);
            if ~isempty(gAnd)
                % turn '(x(410)&x(411))|x(414)|(x(410)&x(411))|x(412)' -> '(x(410)&x(411))| x(412)| x(414)' 
    
                rules=strcat('(x(', num2str(gAnd(1)),')&x(',num2str(gAnd(2)),'))');
                if numel(AndGenes)>2
                   
                    
                else
                    OrGenes=find(OrGenes);
                    for o=1:numel(OrGenes)
                        rules=strcat(rules,'| x(', num2str(g(OrGenes(o))),')');
                    end
                    model.rules(rr)=cellstr(rules);
                    changed_rules_ON(mixed_rules(i))=1;   
                end
            else
                % turns 'x(1798)|x(839)|x(1219)|(x(838)&x(837))|x(837)|x(311)|x(838)|x(839)'
                % to -> 'x(311)|x(837)| x(838)| x(839)| x(1219)| x(1798)'
                OrGenes=find(OrGenes);
                OrGenes=find(OrGenes);
                rules=strcat('x(', num2str(g(OrGenes(1))),')|x(',num2str(g(OrGenes(2))),')');
                for o=3:numel(OrGenes)
                    rules=strcat(rules,'| x(', num2str(g(OrGenes(o))),')');
                end
                model.rules(rr)=cellstr(rules);
                changed_rules_ON(mixed_rules(i))=1;               
            end
            
        else
            
            SuperOR=find(OrGenes==2);
            % rules like: '(x(431))|(x(431)&x(1570))&(x(1332))&(x(1196))'
            % -> x(431) -> cause x(431) is essential, if the other genes
            % are expressed does not really play a role
            if~isempty(SuperOR)
                rules=strcat('x(', num2str(g(SuperOR(1))),')');
                model.rules(rr)=cellstr(rules);
                changed_rules_ON(mixed_rules(i))=1;
                
            end
        end
        
    end
end

if display_changed_mixed_rules 
    [model.rules(find(changed_rules_ON)), model_old_rules(find(changed_rules_ON))]
end

mixed_rules_after=find(contains(model.rules, '&') & contains(model.rules, '|'));
to_check=setdiff(1:numel(model.rxns), mixed_rules_after);
for iv=1:numel(to_check)
    [~,g]=find(model.rxnGeneMat(to_check(iv), :));
    if numel(unique(g))< numel(g)
       'warning';
        
    end
end
% after handling the mixed rules -> there might be new rules which were mixed rules before
% but now after simplifying them they entail only &/| statements 
% -> these we handle now
for i=1:numel(AND_rules)
    g=unique(find(model.rxnGeneMat(AND_rules(i),:)));
    r=strcat('x(',num2str(g(1)),')');
    if numel(g)>1
       for ii=2:numel(g) 
           r =strcat(r,'&x(',num2str(g(ii)),')');
       end
    end
        model.rules(AND_rules(i))=cellstr(r);

end

for i=1:numel(OR_rules)
    g=unique(find(model.rxnGeneMat(OR_rules(i),:)));
    r=strcat('x(',num2str(g(1)),')');
    if numel(g)>1
       for ii=2:numel(g) 
           r =strcat(r,'|x(',num2str(g(ii)),')') ;
       end
    end

    model.rules(OR_rules(i))=cellstr(r);
end
if isfield(model, 'genes')
model=removeUnusedGenes(model);
end

% -- now that we fixed all the rules, and the rxnGeneMat slot
model = creategrRulesField(model);
model = buildRxnGeneMat(model);
end

