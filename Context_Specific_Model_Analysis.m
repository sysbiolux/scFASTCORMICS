%% Create a list where you see directly in how many clusters each reaction was present

%extract = [path,'\Discretization_Step\Extracting_all_possible_information'];

Cluster_Numbers = 1:1:files_length; % Indicates the number of clusters
mkdir ([extract,'\Reactions_Present_in_each_Cluster']);
mkdir ([extract,'\Reaction_Matrix']);
mkdir ([extract,'\Gene_Matrix']);

kk = 0;

Table_Unique_Genes_Asso = [];
Table_Unique_Reactions_Context_Specific_Model = [];

for c=1:numel(Cover)
    for p=1:numel(percent)
        load([extract,'\Context_Specific_Model_Cover_' ,num2str(Cover(c)),'_Percentile_', num2str(percent(p)),'.mat'])
        load([extract,'\A_value_of_the_context_Specific_Model_Cover_' ,num2str(Cover(c)),'_Percentile_', num2str(percent(p)),'.mat'])
        
        Reactions_Context_Specific_Model = model_composite.rxns(A);
        Reactions_Context_Specific_Model = cellfun(@(S) S(1:end-2), Reactions_Context_Specific_Model , 'Uniform', 0);
        Reactions_Context_Specific_Model = array2table(Reactions_Context_Specific_Model);
        
        if isempty(Reactions_Context_Specific_Model) == 0
            Trans_Reactions = contains(Reactions_Context_Specific_Model.Reactions_Context_Specific_Model, 'Trans');
            Trans_Reactions = array2table(Trans_Reactions);
            Reaction_Table = [Reactions_Context_Specific_Model, Trans_Reactions];
            
            Reaction_Table(ismember(Reaction_Table.Trans_Reactions, 1),:)=[];
            Reactions_Context_Specific_Model = Reaction_Table(:,1);
            
            [reactions_count, reactions_names] = histcounts(categorical(Reactions_Context_Specific_Model.Reactions_Context_Specific_Model));
            Table_Reactions_in_each_Cluster = table(reactions_names', reactions_count');
            Header_Table_Reactions_in_each_Cluster = {'Reactions' , 'Present_in_so_many_clusters'};
            Table_Reactions_in_each_Cluster.Properties.VariableNames =  Header_Table_Reactions_in_each_Cluster ;
            
            name_1 = ['\Reactions_in_each_Cluster_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
            Table_Reactions_in_each_Cluster_List = [extract,'\Reactions_Present_in_each_Cluster',name_1,'.txt'];
            writetable(Table_Reactions_in_each_Cluster,Table_Reactions_in_each_Cluster_List);
        end
    

        %% Creating the matrix which can be used for the Jaccard similarity plot
        
        Reactions_Context_Specific_Model_2 = model_composite.rxns(A);
        
        %if isempty(Reactions_Context_Specific_Model) == 0
        % Starting to create the matrix which can be used for the Jaccard similarity plot
        % model_composite.rxns (A) still has _1 behind the reaction
        % unique_reactions_Context_SPecific_Model only contains the unique reactions without _1
        % First searching in model_composite (A) the index
        % Secondly with the help of the index search for the reaction name
        % Finally create a matrix that shows which reaction are present in which cluster
        
        Reactions_Context_Specific_Model_2 = cellfun(@(S) S(1:end-2), Reactions_Context_Specific_Model_2, 'Uniform', 0);
        unique_Reactions_Context_Specific_Model = unique(Reactions_Context_Specific_Model_2);
        unique_Reactions_Context_Specific_Model_matrix = zeros(numel(unique_Reactions_Context_Specific_Model),length(Cluster_Numbers));
        
        for clu=1:numel(Cluster_Numbers)
            Index_of_Reactions_in_Context_Specific_Model = find(contains(model_composite.rxns(A),strcat('_',num2str(Cluster_Numbers(clu)))));
            Index_of_the_reactions_Cluster_in_the_Context = Reactions_Context_Specific_Model_2(Index_of_Reactions_in_Context_Specific_Model);
            unique_Reactions_Context_Specific_Model_matrix(ismember(unique_Reactions_Context_Specific_Model, Index_of_the_reactions_Cluster_in_the_Context),clu)=1;
        end
        
        unique_Reactions_Context_Specific_Model_matrix = array2table(unique_Reactions_Context_Specific_Model_matrix);
        Table_test = [unique_Reactions_Context_Specific_Model, unique_Reactions_Context_Specific_Model_matrix];
        
        % Eliminate the Trans and Ex reactions (not needed)
        if isempty(Reactions_Context_Specific_Model_2) == 0
            Trans_Reactions = contains(Table_test.Var1, 'Trans_');
            Trans_Reactions = array2table(Trans_Reactions);
            Reaction_Table_1 = [Table_test, Trans_Reactions];
            
            Ex_Reactions = contains(Table_test.Var1, 'Ex_');
            Ex_Reactions = array2table(Ex_Reactions);
            Reaction_Table_2 = [Reaction_Table_1, Ex_Reactions];
            
            Reaction_Table_2(ismember(Reaction_Table_2.Trans_Reactions, 1),:)=[];
            Reaction_Table_2(ismember(Reaction_Table_2.Ex_Reactions, 1),:)=[];
            
            Reaction_Matrix = Reaction_Table_2(:,1:end-2);
            
            name_2 = ['\Reaction_Matrix_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
            Reaction_Mat = [extract,'\Reaction_Matrix',name_2,'.txt'];
            writetable(Reaction_Matrix,Reaction_Mat);
        else
            name_2 = ['\Reaction_Matrix_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
            Reaction_Mat = [extract,'\Reaction_Matrix',name_2,'.txt'];
            writetable(Table_test, Reaction_Mat);
        end
  


        %% Genes per cluster

        [~,g]=(find(model_composite.rxnGeneMat(A,:)));
        Genes_in_each_Cluster = unique(model_composite.genes(g));
        
        % All_Genes_Composite_Model = Gene names with _1
        % All_Genes_Composite_Model_Pure = Gene names without _1
        
        All_Genes_Composite_Model = model_composite.genes;
        All_Genes_Composite_Model_Pure= cellfun(@(S) S(1:end-2), All_Genes_Composite_Model  , 'Uniform', 0);
        Unique_All_Genes_Composite_Model_Pure = unique(All_Genes_Composite_Model_Pure);
        Gene_matrix = zeros(numel(Unique_All_Genes_Composite_Model_Pure),length(Cluster_Numbers));
        
        for clu=1:numel(Cluster_Numbers)
            Index_of_Genes_in_Context_Specific_Model = find(contains(Genes_in_each_Cluster,strcat('_',num2str(Cluster_Numbers(clu)))));
            Name_of_the_gene = All_Genes_Composite_Model_Pure(Index_of_Genes_in_Context_Specific_Model);
            Gene_matrix(ismember(Unique_All_Genes_Composite_Model_Pure, Name_of_the_gene),clu)=1;
        end
        
        Finished_Gene_matrix = array2table(Gene_matrix);
        Gene_Table = [Unique_All_Genes_Composite_Model_Pure,Finished_Gene_matrix];
        
        name_3 = ['\Gene_Matrix_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
        Gene_Mat = [extract,'\Gene_Matrix',name_3,'.txt'];
        writetable(Gene_Table,Gene_Mat);
        
  

        %% Gene associated to reactions

        kk = kk + 1;
        
        % First extract the RXN_Gene_Matrix from the context-specific model
        % Extract the unique Gene_Associated to reactions
        
        Rxn_Gene_Matrix_Context_Specific_Model = Context_Specific_model.rxnGeneMat;
        Genes_Asso = Context_Specific_model.genes(find(sum(Rxn_Gene_Matrix_Context_Specific_Model,1) ~=0));
        Genes_Asso = cellfun(@(S) S(1:end-2), Genes_Asso, 'Uniform', 0);
        Genes_Asso_Unique = unique(Genes_Asso);
        
        Table_Unique_Genes_Asso(kk,1) = Cover(c);
        Table_Unique_Genes_Asso(kk,2) = percent(p);
        Table_Unique_Genes_Asso(kk,3) = numel(Genes_Asso_Unique);
        

        %% Reactions in the context-specific model
      
        Rxns_Context_Specific_Model = Context_Specific_model.rxns;
        Rxns_Context_Specific_Model = cellfun(@(S) S(1:end-2), Rxns_Context_Specific_Model  , 'Uniform', 0);
        Unique_Rxns_Context_Specific_Model = unique(Rxns_Context_Specific_Model);
        
        Table_Unique_Reactions_Context_Specific_Model(kk,1) = Cover(c);
        Table_Unique_Reactions_Context_Specific_Model(kk,2) = percent(p);
        Table_Unique_Reactions_Context_Specific_Model(kk,3) = numel(Unique_Rxns_Context_Specific_Model);
        
    end
end

Table_Unique_Genes_Asso = array2table(Table_Unique_Genes_Asso);
Header_Gene_Asso = {'Cover' , 'Percentile' , 'Unique_Gene_Associated'};
Table_Unique_Genes_Asso.Properties.VariableNames = Header_Gene_Asso;
writetable(Table_Unique_Genes_Asso, [path,'\Discretization_Step\Table_Unique_Gene_Asso.txt']);

Table_Unique_Reactions_Context_Specific_Model = array2table(Table_Unique_Reactions_Context_Specific_Model);
Header_Unique_Reactions_Context_Specific_Model  = {'Cover', 'Percentile', 'Unique_Reactions_in_Context_Specific_Model'};
Table_Unique_Reactions_Context_Specific_Model.Properties.VariableNames = Header_Unique_Reactions_Context_Specific_Model;
writetable(Table_Unique_Reactions_Context_Specific_Model, [path,'\Discretization_Step\Table_Unique_Reactions_in_Context_Specific_Model.txt']);

