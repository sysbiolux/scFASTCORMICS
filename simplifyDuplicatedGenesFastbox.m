function[model, changed_rules_ON]=simplifyDuplicatedGenesFastbox(model, taggene)
%% %(c) Maria Pires Pacheco 2023 et al -University of Luxembourg
% INPUTS
% model_target cobra model
%   S                           m x n stoichiometric matrix
%   lb                          n x 1 flux lower bound
%   ub                          n x 1 flux upper bound
%   rxns                        n x 1 cell array of reaction abbreviations
%   c                           n x 1

%  taggene == 1 if the model is a Recon model and the .1 has to be removed
if isfield(model, 'rxnGeneMat')
    if numel(model.genes)~=size(model.rxnGeneMat,2)
        disp('simplifyDuplicatedGenes_fastbox was not run') 
        return
    end
end

[model]=removeUnusedGenesFastbox(model,taggene);
[genes,~,~] = unique(model.genes);
if numel(genes) < numel(model.genes)
    %find non unique genes
    [genes,ia,ic] = unique(model.genes);
    for i=1:numel(ia)        
        dupidx = find(ic == i);
        if numel(dupidx)>1            
            first=min(dupidx);
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
    Index=strfind(model.rules, '&') ;
    match=(cellfun('isempty',Index) & ~strcmp(model.rules, '')); % find statement only with OR
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
    %
    model.rules(match)=new2;
    Index=strfind(model.rules, '&') ;
    %
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

changed_rules_ON=zeros(numel(mixed_rules),1);
for i=1:numel(mixed_rules)
    %%
    
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
            
            gAnd=g(AndGenes);
            if ~isempty(gAnd)
                rules=strcat('(x(', num2str(gAnd(1)),')&x(',num2str(gAnd(2)),'))');
                if numel(AndGenes)>2
                 
                    
                else
                    OrGenes=find(OrGenes);
                    for o=1:numel(OrGenes)
                        rules=strcat(rules,'| x(', num2str(g(OrGenes(o))),')');
                    end
                    model.rules(rr)=cellstr(rules);
                    
                end
            else
                OrGenes=find(OrGenes);
                OrGenes=find(OrGenes);
                
                rules=strcat('x(', num2str(g(OrGenes(1))),')|x(',num2str(g(OrGenes(2))),')');
                for o=3:numel(OrGenes)
                    rules=strcat(rules,'| x(', num2str(g(OrGenes(o))),')');
                end
                model.rules(rr)=cellstr(rules);
                model.rules(rr)=cellstr(rules);
                changed_rules_ON(mixed_rules(i))=1;               
            end
            
        else
            SuperOR=find(OrGenes==2);
            if~isempty(SuperOR)
                rules=strcat('x(', num2str(g(SuperOR(1))),')');
                model.rules(rr)=cellstr(rules);
                changed_rules_ON(mixed_rules(i))=1;
                
            end
        end
        %GLCNACT2g
        
    end
end


mixed_rules_after=find(contains(model.rules, '&') & contains(model.rules, '|'));
to_check=setdiff(1:numel(model.rxns), mixed_rules_after);
for iv=1:numel(to_check)
    [~,g]=find(model.rxnGeneMat(to_check(iv), :));
    if numel(unique(g))< numel(g)
       'warning';
        
    end
end
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
model = creategrRulesField(model);
model = buildRxnGeneMat(model);
end

