function[Analysis,T]=Context_Specific_Model_Analysis(model_composite,Optimization_global,files_length, path, Cover_range,REI_range, printLevel,T)
%% Create a list where you see directly in how many clusters each reaction was present

extract = [path,'\Discretization_Step\Extracting_all_possible_information'];

Cluster_Numbers = 1:1:files_length; % Indicates the number of clusters
if printLevel==1
mkdir ([extract,'\Reactions_Present_in_each_Cluster']);
mkdir ([extract,'\Reaction_Matrix']);
mkdir ([extract,'\Gene_Matrix']);
end
kk = 0;

k=0;
Analysis=struct();

for c=1:numel(Cover_range)
    for p=1:numel(REI_range)
        k=k+1;
        A=find(Optimization_global(1).A(:,k));
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
            if printLevel==1
            name_1 = ['\Reactions_in_each_Cluster_Cover_',num2str(Cover_range(c)),'_Percentile_',num2str(REI_range(p))];
            Table_Reactions_in_each_Cluster_List = [extract,'\Reactions_Present_in_each_Cluster',name_1,'.txt'];
            writetable(Table_Reactions_in_each_Cluster,Table_Reactions_in_each_Cluster_List);
            end
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
            Index_of_Reactions_in_Context_Specific_Model = (contains(model_composite.rxns(A),strcat('_',num2str(Cluster_Numbers(clu)))));
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
            if printLevel==1
            name_2 = ['\Reaction_Matrix_Cover_',num2str(Cover_range(c)),'_Percentile_',num2str(REI_range(p))];
            Reaction_Mat = [extract,'\Reaction_Matrix',name_2,'.txt'];
            writetable(Reaction_Matrix,Reaction_Mat);
            end
             Analysis(k).Reaction_mat=Reaction_Matrix;

        else
            if printLevel==1
            name_2 = ['\Reaction_Matrix_Cover_',num2str(Cover_range(c)),'_Percentile_',num2str(REI_range(p))];
            Reaction_Mat = [extract,'\Reaction_Matrix',name_2,'.txt'];
            
            writetable(Table_test, Reaction_Mat);
            end
            Analysis(k).Reaction_mat=Reaction_Matrix;
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
            Index_of_Genes_in_Context_Specific_Model = (contains(Genes_in_each_Cluster,strcat('_',num2str(Cluster_Numbers(clu)))));
            Name_of_the_gene = All_Genes_Composite_Model_Pure(Index_of_Genes_in_Context_Specific_Model);
            Gene_matrix(ismember(Unique_All_Genes_Composite_Model_Pure, Name_of_the_gene),clu)=1;
        end
        
        Finished_Gene_matrix = array2table(Gene_matrix);
        Gene_Table = [Unique_All_Genes_Composite_Model_Pure,Finished_Gene_matrix];
        if printLevel==1
        name_3 = ['\Gene_Matrix_Cover_',num2str(Cover_range(c)),'_Percentile_',num2str(REI_range(p))];
        Gene_Mat = [extract,'\Gene_Matrix',name_3,'.txt'];
        writetable(Gene_Table,Gene_Mat);
        end
        Analysis(k).gene_mat=Gene_Table;
  

        %% Gene associated to reactions

        kk = kk + 1;
        
        % First extract the RXN_Gene_Matrix from the context-specific model
        % Extract the unique Gene_Associated to reactions
        
       
        
        
        

        %% Reactions in the context-specific model
      
        
    end
end


%writetable(Table_Unique_Genes_Asso, [path,'\Discretization_Step\Table_Unique_Gene_Asso.txt']);


end