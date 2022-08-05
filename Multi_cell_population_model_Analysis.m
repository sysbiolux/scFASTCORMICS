function[Analysis,T]=Multi_cell_population_model_Analysis(Expanded_input_model,Optimization_global,files_length, path, coverage_range,REI_range, printLevel,T)
% (c) Maria Pires Pacheco 2022 et al -University of Luxembourg
%% Create a list where you see directly in how many clusters each reaction was present
%Compute a reaction and a gene matrix

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

for c=1:numel(coverage_range)
    for p=1:numel(REI_range)
        k=k+1;
        A=find(Optimization_global(1).A(:,k));
        Reactions_Multi_cell_population_model = Expanded_input_model.rxns(A);
        Reactions_Multi_cell_population_model = cellfun(@(S) S(1:end-2), Reactions_Multi_cell_population_model , 'Uniform', 0);
        Reactions_Multi_cell_population_model = array2table(Reactions_Multi_cell_population_model);
        
        if isempty(Reactions_Multi_cell_population_model) == 0
            Trans_Reactions = contains(Reactions_Multi_cell_population_model.Reactions_Multi_cell_population_model, 'Trans');
            Trans_Reactions = array2table(Trans_Reactions);
            Reaction_Table = [Reactions_Multi_cell_population_model, Trans_Reactions];
            
            Reaction_Table(ismember(Reaction_Table.Trans_Reactions, 1),:)=[];
            Reactions_Multi_cell_population_model = Reaction_Table(:,1);
            
            [reactions_count, reactions_names] = histcounts(categorical(Reactions_Multi_cell_population_model.Reactions_Multi_cell_population_model));
            Table_Reactions_in_each_Cluster = table(reactions_names', reactions_count');
            Header_Table_Reactions_in_each_Cluster = {'Reactions' , 'Present_in_so_many_clusters'};
            Table_Reactions_in_each_Cluster.Properties.VariableNames =  Header_Table_Reactions_in_each_Cluster ;
            if printLevel==1
                name_1 = ['\Reactions_in_each_Cluster_Cover_',num2str(coverage_range(c)),'_Percentile_',num2str(REI_range(p))];
                Table_Reactions_in_each_Cluster_List = [extract,'\Reactions_Present_in_each_Cluster',name_1,'.txt'];
                writetable(Table_Reactions_in_each_Cluster,Table_Reactions_in_each_Cluster_List);
            end
        end
        
        
        %% Creating the matrix which can be used for the Jaccard similarity plot
        
        Reactions_Multi_cell_population_model_2 = Expanded_input_model.rxns(A);         
        Reactions_Multi_cell_population_model_2 = cellfun(@(S) S(1:end-2), Reactions_Multi_cell_population_model_2, 'Uniform', 0);
        unique_Reactions_Multi_cell_population_model = unique(Reactions_Multi_cell_population_model_2);
        unique_Reactions_Multi_cell_population_model_matrix = zeros(numel(unique_Reactions_Multi_cell_population_model),length(Cluster_Numbers));
        
        for clu=1:numel(Cluster_Numbers)
            Index_of_Reactions_in_Multi_cell_population_model = (contains(Expanded_input_model.rxns(A),strcat('_',num2str(Cluster_Numbers(clu)))));
            Index_of_the_reactions_Cluster_in_the_Context = Reactions_Multi_cell_population_model_2(Index_of_Reactions_in_Multi_cell_population_model);
            unique_Reactions_Multi_cell_population_model_matrix(ismember(unique_Reactions_Multi_cell_population_model, Index_of_the_reactions_Cluster_in_the_Context),clu)=1;
        end
        
        unique_Reactions_Multi_cell_population_model_matrix = array2table(unique_Reactions_Multi_cell_population_model_matrix);
        Table_test = [unique_Reactions_Multi_cell_population_model, unique_Reactions_Multi_cell_population_model_matrix];
        
        % Eliminate the Trans and Ex reactions (not needed)
        if isempty(Reactions_Multi_cell_population_model_2) == 0
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
                name_2 = ['\Reaction_Matrix_Cover_',num2str(coverage_range(c)),'_Percentile_',num2str(REI_range(p))];
                Reaction_Mat = [extract,'\Reaction_Matrix',name_2,'.txt'];
                writetable(Reaction_Matrix,Reaction_Mat);
            end
            Analysis(k).Reaction_mat=Reaction_Matrix;
            
        else
            if printLevel==1
                name_2 = ['\Reaction_Matrix_Cover_',num2str(coverage_range(c)),'_Percentile_',num2str(REI_range(p))];
                Reaction_Mat = [extract,'\Reaction_Matrix',name_2,'.txt'];
                
                writetable(Table_test, Reaction_Mat);
            end
            Analysis(k).Reaction_mat=Reaction_Matrix;
        end
        
        
        
        %% Genes per cluster
        
        [~,g]=(find(Expanded_input_model.rxnGeneMat(A,:)));
        Genes_in_each_Cluster = unique(Expanded_input_model.genes(g));
        
        % All_Genes_Composite_Model = Gene names with _1
        % All_Genes_Composite_Model_Pure = Gene names without _1
        
        All_Genes_Composite_Model = Expanded_input_model.genes;
        All_Genes_Composite_Model_Pure= cellfun(@(S) S(1:end-2), All_Genes_Composite_Model  , 'Uniform', 0);
        Unique_All_Genes_Composite_Model_Pure = unique(All_Genes_Composite_Model_Pure);
        Gene_matrix = zeros(numel(Unique_All_Genes_Composite_Model_Pure),length(Cluster_Numbers));
        
        for clu=1:numel(Cluster_Numbers)
            Index_of_Genes_in_Multi_cell_population_model = (contains(Genes_in_each_Cluster,strcat('_',num2str(Cluster_Numbers(clu)))));
            Name_of_the_gene = All_Genes_Composite_Model_Pure(Index_of_Genes_in_Multi_cell_population_model);
            Gene_matrix(ismember(Unique_All_Genes_Composite_Model_Pure, Name_of_the_gene),clu)=1;
        end
        
        Finished_Gene_matrix = array2table(Gene_matrix);
        Gene_Table = [Unique_All_Genes_Composite_Model_Pure,Finished_Gene_matrix];
        if printLevel==1
            name_3 = ['\Gene_Matrix_Cover_',num2str(coverage_range(c)),'_Percentile_',num2str(REI_range(p))];
            Gene_Mat = [extract,'\Gene_Matrix',name_3,'.txt'];
            writetable(Gene_Table,Gene_Mat);
        end
        Analysis(k).gene_mat=Gene_Table;       
              
        kk = kk + 1;       
    end
end
end