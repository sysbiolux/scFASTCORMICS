% Determine the reactions in context specific and composite model that should
% or should not be there according to bulk data

Discretization_Table = Discretization_Table(find(~isnan(table2array(Discretization_Table(:,2)))),:);

dico_biodbnet = [string(biodbnet.GeneID), biodbnet.GeneSymbol];
[r,c] = find(cellfun('isempty',dico_biodbnet));
dico_biodbnet(r,:) = [];
dico_biodbnet = cell2table(cellstr(dico_biodbnet));

Values_Discretization_Table = table2array(Discretization_Table(:,4:end));
Rownames_Bulk_data = cellstr(string(Discretization_Table.GeneID));

% extract = [path,'\Discretization_Step\Extracting_all_possible_information'];

Genes_Composite_Model = model_composite.genes;
Genes_Composite_Model = cellfun(@(S) S(1:end-2), Genes_Composite_Model, 'Uniform', 0);

rxns_should_bulk_data = Discretization_Table(find(sum(Values_Discretization_Table,2) == size(Values_Discretization_Table,2)),:);
rxns_should_not_bulk_data = Discretization_Table(find(sum(Values_Discretization_Table,2) == size(Values_Discretization_Table,2)*(-1)),:);

writetable(rxns_should_bulk_data,[dis,'\Table_Reaction_Should_Bulk.txt']);
writetable(rxns_should_not_bulk_data,[dis,'\Table_Reaction_Should_not_Bulk.txt']);

Table_Reaction_Expressed = [];
Table_Reaction_Never_Expressed = [];

kk = 0;

for c=1:numel(Cover)
    for p=1:numel(percent)
        
        kk = kk + 1;
        load([extract,'\Context_Specific_Model_Cover_' ,num2str(Cover(c)),'_Percentile_', num2str(percent(p)),'.mat'])
        load([extract,'\A_value_of_the_context_Specific_Model_Cover_' ,num2str(Cover(c)),'_Percentile_', num2str(percent(p)),'.mat'])
        
        Context_Specific_model.genes = Genes_Composite_Model;
        
        [mapping] = map_expression_2_data_rFASCTCORMICS(Context_Specific_model, Values_Discretization_Table, dico_biodbnet, Rownames_Bulk_data);
        
        mapping_reactions = [path,'\Discretization_Step\mapping_reactions_',num2str(i)];
        save(mapping_reactions,'mapping')
        
        Rxns_Not_Present = Context_Specific_model.rxns(find(sum(mapping,2) == size(mapping,2)*(-1)));
        Rxns_Present = Context_Specific_model.rxns(find(sum(mapping,2) == size(mapping,2)));
        
        Rxns_Not_Present = cellfun(@(S) S(1:end-2), Rxns_Not_Present, 'Uniform', 0);
        Rxns_Present = cellfun(@(S) S(1:end-2), Rxns_Present , 'Uniform', 0);
        
        Rxns_Not_Present = unique(Rxns_Not_Present);
        Rxns_Present = unique(Rxns_Present);
        
        Table_Reaction_Expressed(kk,1) = Cover(c);
        Table_Reaction_Expressed(kk,2) = percent(p);
        Table_Reaction_Expressed(kk,3) = numel(Rxns_Present);
        
        Table_Reaction_Never_Expressed(kk,1) = Cover(c);
        Table_Reaction_Never_Expressed(kk,2) = percent(p);
        Table_Reaction_Never_Expressed(kk,3) = numel(Rxns_Not_Present);
        
    end
end

Table_Reaction_Expressed = array2table(Table_Reaction_Expressed );
Table_Reaction_Never_Expressed = array2table (Table_Reaction_Never_Expressed);

Header_Reaction_Expressed = {'Cover', 'Percentile', 'Number_Reactions_that_should_be_there'};
Header_Reaction_Never_Expressed = {'Cover', 'Percentile', 'Number_Reactions_that_should_not_be_there'};

Table_Reaction_Expressed.Properties.VariableNames = Header_Reaction_Expressed ;
Table_Reaction_Never_Expressed.Properties.VariableNames = Header_Reaction_Never_Expressed;

writetable(Table_Reaction_Expressed, [dis,'\Table_Number_of_reactions_should_be_there.txt']);
writetable(Table_Reaction_Never_Expressed, [dis,'\Table_Number_of_reactions_should_not_be_there.txt']);


%% Determine the number of reactions which should or should not be there according to bulk data

model = model_composite;
model.genes = Genes_Composite_Model;

Rownames_Bulk_data = cellstr(string(Discretization_Table.GeneID));
[mapping] = map_expression_2_data_rFASCTCORMICS(model, Values_Discretization_Table, dico_biodbnet, Rownames_Bulk_data);

Rxns_Not_Present = model_composite.rxns(find(sum(mapping,2) == size(mapping,2)*(-1)));
Rxns_Present = model_composite.rxns(find(sum(mapping,2) == size(mapping,2)));

Rxns_Not_Present = cellfun(@(S) S(1:end-2), Rxns_Not_Present, 'Uniform', 0);
Rxns_Present = cellfun(@(S) S(1:end-2), Rxns_Present, 'Uniform', 0);

Rxns_Not_Present = unique(Rxns_Not_Present);
Rxns_Present = unique(Rxns_Present);

Table_Reactions_Presence = {'Should', 'Should not be here'; numel(Rxns_Present), numel(Rxns_Not_Present)};

writecell(Table_Reactions_Presence, [dis,'\Table_Reactions_Presence_Count.txt']);
writecell(Rxns_Not_Present, [dis,'\Table_Reactions_Not_Present.txt']);
writecell(Rxns_Present, [dis,'\Table_Reactions_Present.txt']);
