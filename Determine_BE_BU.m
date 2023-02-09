function[T]=Determine_BE_BU(Expanded_Input_model,Discretization_Table, biodbnet, coverage_range, REI_range, Optimization_global, T)
% (c) Maria Pires Pacheco 2022 et al -University of Luxembourg

% Determine the reactions in the expanded input model and multi-cell population models that should
% or should not be there according to bulk data
Discretization_Table = Discretization_Table((~isnan(table2array(Discretization_Table(:,2)))),:);

dico_biodbnet = [string(biodbnet.GeneID), biodbnet.GeneSymbol];
[r,~] = find(cellfun('isempty',dico_biodbnet));
dico_biodbnet(r,:) = [];
dico_biodbnet = cell2table(cellstr(dico_biodbnet));

Values_Discretization_Table = table2array(Discretization_Table(:,4:end));
Rownames_Bulk_data = cellstr(string(Discretization_Table.GeneID));


Genes_Expanded_Input_model = Expanded_Input_model.genes;
Genes_Expanded_Input_model = cellfun(@(S) S(1:end-2), Genes_Expanded_Input_model, 'Uniform', 0);
Expanded_Input_model_tmp=Expanded_Input_model;
Expanded_Input_model_tmp.genes=Genes_Expanded_Input_model;
%writetable(rxns_should_bulk_data,[dis,'\Table_Reaction_Should_Bulk.txt']);
%writetable(rxns_should_not_bulk_data,[dis,'\Table_Reaction_Should_not_Bulk.txt']);

Table_Reaction_Expressed =zeros(numel(REI_range)*numel(coverage_range),3);
Table_Reaction_Never_Expressed = zeros(numel(REI_range)*numel(coverage_range),3);

kk = 0;
[mapping] = map_expression_2_data_rFASTCORMICS(Expanded_Input_model_tmp, Values_Discretization_Table, dico_biodbnet, Rownames_Bulk_data);
Rxns_Not_Present = Expanded_Input_model.rxns(sum(mapping,2) <= (size(mapping,2)*(-0.9)));
Rxns_Present = Expanded_Input_model.rxns(sum(mapping,2) > (size(mapping,2) * 0.9));

Rxns_Not_Present = cellfun(@(S) S(1:end-2), Rxns_Not_Present, 'Uniform', 0);
Rxns_Present = cellfun(@(S) S(1:end-2), Rxns_Present , 'Uniform', 0);
Rxns_Not_Present = unique(Rxns_Not_Present);
Rxns_Present = unique(Rxns_Present);
T(1).Rxns_Present=Rxns_Present;
T(1).Rxns_Not_Present=Rxns_Not_Present;
Expanded_Input_model_tmp.rxns=cellfun(@(S) S(1:end-2), Expanded_Input_model_tmp.rxns, 'Uniform', 0);
for c=1:numel(coverage_range)
    for p=1:numel(REI_range)
        
        kk = kk + 1;
        Multi_cell_population_model=Expanded_Input_model_tmp.rxns(find(Optimization_global(1).A(:,kk)));
        
     
        Rxns_Present_c=intersect(Multi_cell_population_model,Rxns_Present);
        Rxns_Not_Present_c=intersect(Multi_cell_population_model,Rxns_Not_Present);

        
        Table_Reaction_Expressed(kk,1) = coverage_range(c);
        Table_Reaction_Expressed(kk,2) = REI_range(p);
        Table_Reaction_Expressed(kk,3) = numel(Rxns_Present_c);
        
        Table_Reaction_Never_Expressed(kk,1) = coverage_range(c);
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

model = Expanded_Input_model;

Rownames_Bulk_data = cellstr(string(Discretization_Table.GeneID));
[mapping] = map_expression_2_data_rFASTCORMICS(model, Values_Discretization_Table, dico_biodbnet, Rownames_Bulk_data);

Rxns_Not_Present = Expanded_Input_model.rxns((sum(mapping,2) == size(mapping,2)*(-1)));
Rxns_Present = Expanded_Input_model.rxns((sum(mapping,2) == size(mapping,2)));

Rxns_Not_Present = cellfun(@(S) S(1:end-2), Rxns_Not_Present, 'Uniform', 0);
Rxns_Present = cellfun(@(S) S(1:end-2), Rxns_Present, 'Uniform', 0);

Rxns_Not_Present = unique(Rxns_Not_Present);
Rxns_Present = unique(Rxns_Present);
T.Reactions_Present=Rxns_Present;
T.Reactions_Not_Present=Rxns_Not_Present;

Table_Reactions_Presence = {'Should', 'Should not be here'; numel(Rxns_Present), numel(Rxns_Not_Present)};
T.Reactions_Presence=Table_Reactions_Presence;

end
