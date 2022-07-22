function[T]=Intersection_Reactions_Genes(model_composite,T, Cover_range, REI_range, Optimization_global, printLevel)
% From the composite model, we have the names of present and not present
% reactions/genes. Now, we want to intersect these reactions from the composite
% model with those from every context specific model.

%% Intersection of context model reactions and composite model reactions

comp_rxns_present = T.Reactions_Present;
comp_rxns_not_present =T.Reactions_Not_Present;


Table_spec_rxns_present = cell(100,3);
Table_spec_rxns_not_present = cell(100,3);
Table_spec_genes_present = cell(100,3);
Table_spec_genes_not_present = cell(100,3);

kk=  0;

for c=1:numel(Cover_range)
    for p=1:numel(REI_range)
        
        kk = kk + 1;
       % A=Optimization(kk).A
        Context_Specific_model_rxns=model_composite.rxns(find(Optimization_global(1).A(:,kk)));

       % load([extract,'\A_value_of_the_context_Specific_Model_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat']);
        %load([extract,'\Context_Specific_Model_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat']);
        
        %%%% Determination of bulk rxns in context specific model
        spec_model_rxns = cellfun(@(S) S(1:end-2), Context_Specific_model_rxns, 'Uniform', 0);
        spec_rxns_present = intersect(spec_model_rxns, comp_rxns_present);
        spec_rxns_not_present = intersect(spec_model_rxns, comp_rxns_not_present);
        
        Table_spec_rxns_present{kk,1} = Cover_range(c);
        Table_spec_rxns_present{kk,2} = REI_range(p);
        Table_spec_rxns_present{kk,3} = numel(spec_rxns_present);
        
        Table_spec_rxns_not_present{kk,1} = Cover_range(c);
        Table_spec_rxns_not_present{kk,2} = REI_range(p);
        Table_spec_rxns_not_present{kk,3} = numel(spec_rxns_not_present);
        %%%%%%%%%
        
        %%%% Determination of bulk genes in context specific model
%        spec_model_genes = Context_Specific_model.genes((sum(Context_Specific_model.rxnGeneMat,1)~=0));
 %       spec_model_genes = cellfun(@(S) S(1:end-2), spec_model_genes, 'Uniform', 0);
        
  %      spec_genes_present = intersect(string(spec_model_genes), string(comp_genes_present));
   %     spec_genes_not_present = intersect(string(spec_model_genes), string(comp_genes_not_present));
        
    %    Table_spec_genes_present{kk,1} = Cover_range(c); 
     %   Table_spec_genes_present{kk,2} = REI_range(p);
      %  Table_spec_genes_present{kk,3} = numel(spec_genes_present);
        
       % Table_spec_genes_not_present{kk,1} = Cover_range(c);
        %Table_spec_genes_not_present{kk,2} = REI_range(p);
        %Table_spec_genes_not_present{kk,3} = numel(spec_genes_not_present);
        %%%%%%%%%
        
    end
end
T.Table_spec_rxns_present=Table_spec_rxns_present;
T.Table_spec_rxns_not_present=Table_spec_rxns_not_present;
T.Table_spec_genes_present=Table_spec_genes_present;
T.Table_spec_genes_not_present=Table_spec_genes_not_present;
if printLevel==1
    
% writecell(Table_spec_rxns_present,[dis,'\Table_intersect_rxns_present.txt'])
% writecell(Table_spec_rxns_not_present,[dis,'\Table_intersect_rxns_not_present.txt'])
% 
% writecell(Table_spec_genes_present,[dis,'\Table_intersect_genes_present.txt'])
% writecell(Table_spec_genes_not_present,[dis,'\Table_intersect_genes_not_present.txt'])
 end

