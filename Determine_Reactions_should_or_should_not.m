function[T]=Determine_Reactions_should_or_should_not(model_composite,Discretization_Table, biodbnet, Cover_range, REI_range, Optimization_global, printLevel, T)
% Determine the reactions in context specific and composite model that should
% or should not be there according to bulk data
Discretization_Table = Discretization_Table((~isnan(table2array(Discretization_Table(:,2)))),:);

dico_biodbnet = [string(biodbnet.GeneID), biodbnet.GeneSymbol];
[r,~] = find(cellfun('isempty',dico_biodbnet));
dico_biodbnet(r,:) = [];
dico_biodbnet = cell2table(cellstr(dico_biodbnet));

Values_Discretization_Table = table2array(Discretization_Table(:,4:end));
Rownames_Bulk_data = cellstr(string(Discretization_Table.GeneID));


Genes_Composite_Model = model_composite.genes;
Genes_Composite_Model = cellfun(@(S) S(1:end-2), Genes_Composite_Model, 'Uniform', 0);
model_composite_tmp=model_composite;
model_composite_tmp.genes=Genes_Composite_Model;
%writetable(rxns_should_bulk_data,[dis,'\Table_Reaction_Should_Bulk.txt']);
%writetable(rxns_should_not_bulk_data,[dis,'\Table_Reaction_Should_not_Bulk.txt']);

Table_Reaction_Expressed =zeros(numel(REI_range)*numel(Cover_range),3);
Table_Reaction_Never_Expressed = zeros(numel(REI_range)*numel(Cover_range),3);

kk = 0;
[mapping] = map_expression_2_data_rFASCTCORMICS(model_composite_tmp, Values_Discretization_Table, dico_biodbnet, Rownames_Bulk_data);
Rxns_Not_Present = model_composite.rxns(sum(mapping,2) == (size(mapping,2)*(-1)));
Rxns_Present = model_composite.rxns(sum(mapping,2) == size(mapping,2));

Rxns_Not_Present = cellfun(@(S) S(1:end-2), Rxns_Not_Present, 'Uniform', 0);
Rxns_Present = cellfun(@(S) S(1:end-2), Rxns_Present , 'Uniform', 0);
Rxns_Not_Present = unique(Rxns_Not_Present);
Rxns_Present = unique(Rxns_Present);
T(1).Rxns_Present=Rxns_Present;
T(1).Rxns_Not_Present=Rxns_Not_Present;
model_composite_tmp.rxns=cellfun(@(S) S(1:end-2), model_composite_tmp.rxns, 'Uniform', 0);
for c=1:numel(Cover_range)
    for p=1:numel(REI_range)
        
        kk = kk + 1;
        Context_Specific_model_rxns=model_composite_tmp.rxns(find(Optimization_global(1).A(:,kk)));
        
        %load([extract,'\Context_Specific_Model_Cover_' ,num2str(Cover_range(c)),'_Percentile_', num2str(REI_range(p)),'.mat'])
        %load([extract,'\A_value_of_the_context_Specific_Model_Cover_' ,num2str(Cover_range(c)),'_Percentile_', num2str(REI_range(p)),'.mat'])
        
        
        Rxns_Present_c=intersect(Context_Specific_model_rxns,Rxns_Present);
        Rxns_Not_Present_c=intersect(Context_Specific_model_rxns,Rxns_Not_Present);

        
        Table_Reaction_Expressed(kk,1) = Cover_range(c);
        Table_Reaction_Expressed(kk,2) = REI_range(p);
        Table_Reaction_Expressed(kk,3) = numel(Rxns_Present_c);
        
        Table_Reaction_Never_Expressed(kk,1) = Cover_range(c);
        Table_Reaction_Never_Expressed(kk,2) = REI_range(p);
        Table_Reaction_Never_Expressed(kk,3) = numel(Rxns_Not_Present_c);
        
    end
end

Table_Reaction_Expressed = array2table(Table_Reaction_Expressed );
Table_Reaction_Never_Expressed = array2table (Table_Reaction_Never_Expressed);

Header_Reaction_Expressed = {'Cover', 'Percentile', 'Number_Reactions_that_should_be_there'};
Header_Reaction_Never_Expressed = {'Cover', 'Percentile', 'Number_Reactions_that_should_not_be_there'};

Table_Reaction_Expressed.Properties.VariableNames = Header_Reaction_Expressed ;
Table_Reaction_Never_Expressed.Properties.VariableNames = Header_Reaction_Never_Expressed;
T.Table_Reaction_Expressed=Table_Reaction_Expressed;
T.Table_Reaction_Never_Expressed=Table_Reaction_Never_Expressed;
%writetable(Table_Reaction_Expressed, [dis,'\Table_Number_of_reactions_should_be_there.txt']);
%writetable(Table_Reaction_Never_Expressed, [dis,'\Table_Number_of_reactions_should_not_be_there.txt']);


%% Determine the number of reactions which should or should not be there according to bulk data

model = model_composite;

Rownames_Bulk_data = cellstr(string(Discretization_Table.GeneID));
[mapping] = map_expression_2_data_rFASCTCORMICS(model, Values_Discretization_Table, dico_biodbnet, Rownames_Bulk_data);

Rxns_Not_Present = model_composite.rxns((sum(mapping,2) == size(mapping,2)*(-1)));
Rxns_Present = model_composite.rxns((sum(mapping,2) == size(mapping,2)));

Rxns_Not_Present = cellfun(@(S) S(1:end-2), Rxns_Not_Present, 'Uniform', 0);
Rxns_Present = cellfun(@(S) S(1:end-2), Rxns_Present, 'Uniform', 0);

Rxns_Not_Present = unique(Rxns_Not_Present);
Rxns_Present = unique(Rxns_Present);
T.Reactions_Present=Rxns_Present;
T.Reactions_Not_Present=Rxns_Not_Present;

Table_Reactions_Presence = {'Should', 'Should not be here'; numel(Rxns_Present), numel(Rxns_Not_Present)};
T.Reactions_Presence=Table_Reactions_Presence;

end