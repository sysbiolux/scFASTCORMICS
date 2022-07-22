function[T_bulk]=Determine_Genes_should_or_should_not(Discretization_Table,  path, Cover_range, REI_range,  printLevel,Metabolic_input_model, Analysis)
% Create from the bulk data : a list with expressed genes (in at least 90% of the samples) + a list of never expressed genes
dis = [path,'\Discretization_Step'];
Discretization_Table_Values = table2array(Discretization_Table(:,4:end));
Sum_Rows_Discretization_Table = sum(Discretization_Table_Values, 2);

expressed = ceil(0.9*size(Discretization_Table_Values,2));
n_expressed = -size(Discretization_Table_Values,2);

Expressed_Genes = Discretization_Table.GeneID((Sum_Rows_Discretization_Table >= expressed));
Expressed_Genes = cellstr(string(Expressed_Genes));


Metabolic_Expressed_Genes = Expressed_Genes((ismember((Expressed_Genes), (Metabolic_input_model))));
T_bulk.Gene_ID_Expressed_Bulk_Data=Metabolic_Expressed_Genes;
T_bulk.Gene_ID_Expressed_Bulk_Data=Metabolic_Expressed_Genes;

if printLevel==1
    writecell(Metabolic_Expressed_Genes, [dis,'\Gene_ID_Expressed_Bulk_Data.txt']);
end
Never_Expressed_Genes = Discretization_Table.GeneID((Sum_Rows_Discretization_Table == n_expressed));
Never_Expressed_Genes = cellstr(string(Never_Expressed_Genes));


Metabolic_Never_Expressed_Genes = Never_Expressed_Genes((ismember(Never_Expressed_Genes, Metabolic_input_model)));
T_bulk.Gene_ID_Never_Expressed_Bulk_Data=Metabolic_Never_Expressed_Genes;
if printLevel==1
    writecell(Metabolic_Never_Expressed_Genes, [dis,'\Gene_ID_Never_Expressed_Bulk_Data.txt']);
end

%% Compare your expressed and never expressed gene list with your sc-data
kk = 0;
Table_Expressed_Genes = zeros(numel(Cover_range)*numel(REI_range),3);
Table_Gene_Never_Expressed = zeros(numel(Cover_range)*numel(REI_range),3);
Header_Count = {'Count_of_Genes'};

for c=1:numel(Cover_range)
    for p=1:numel(REI_range)
        kk = kk + 1;
        
        
        % Gene_Matrix = readtable([extract,'\Gene_Matrix\Gene_Matrix_Cover_' ,num2str(Cover_range(c)),'_REI_', num2str(REI(p)),'.txt']);
        Gene_Matrix=Analysis(kk).gene_mat;
        Gene_Matrix_Values = table2array(Gene_Matrix(:,2:end));
        Sum_Rows = array2table(sum(Gene_Matrix_Values,2)); % Sum of the clusters to see how often the gene is present
        Sum_Rows.Properties.VariableNames = Header_Count;
        Gene_List = [Gene_Matrix, Sum_Rows];
        
        % Expressed genes between bulk and sc-RNA
        % First make a comparison between the sc-RNA data and the bulk data
        % Change line 67 depending if you work with gene ID or gene symbol
        Number_Expressed_Genes = array2table(ismember(Gene_List.Var1,Metabolic_Expressed_Genes));
        Header_1 = {'Bulk_Sc_Comparison'};
        Number_Expressed_Genes.Properties.VariableNames = Header_1;
        Gene_List_Expressed = [Gene_List, Number_Expressed_Genes];
        
        % Eliminate all the genes which are only present in the sc-data but not in the bulk
        
        Gene_List_Expressed = Gene_List_Expressed((Gene_List_Expressed.Bulk_Sc_Comparison),:);
        
        % Check now if the gene is present in at least one cluster
        Gene_List_Expressed = Gene_List_Expressed.Var1((Gene_List_Expressed.Count_of_Genes >= 1));
        
        
        Table_Expressed_Genes(kk,1) = Cover_range(c);
        Table_Expressed_Genes(kk,2) = REI_range(p);
        Table_Expressed_Genes(kk,3) = numel(Gene_List_Expressed);
        
        
        % Never Expressed genes between bulk and sc-RNA
        % First make a comparaison between the sc-RNA data and the bulk data
        
        Number_Never_Expressed_Genes = array2table(ismember(Gene_List.Var1,Metabolic_Never_Expressed_Genes));
        Number_Never_Expressed_Genes.Properties.VariableNames = Header_1;
        Gene_List_Never_Expressed = [Gene_List, Number_Never_Expressed_Genes];
        
        % Eliminate all the genes which are only present in the sc-data but not in the bulk
        
        Gene_List_Never_Expressed = Gene_List_Never_Expressed((Gene_List_Never_Expressed.Bulk_Sc_Comparison),:);
        
        % Check now if the gene is present in at least one cluster
        
        Gene_List_Never_Expressed = Gene_List_Never_Expressed.Var1((Gene_List_Never_Expressed.Count_of_Genes >= 1));
        
        Table_Gene_Never_Expressed(kk,1) = Cover_range(c);
        Table_Gene_Never_Expressed(kk,2) = REI_range(p);
        Table_Gene_Never_Expressed(kk,3) = numel(Gene_List_Never_Expressed);
    end
end


Table_Expressed_Genes = array2table(Table_Expressed_Genes);
Table_Gene_Never_Expressed = array2table(Table_Gene_Never_Expressed);

Header_Table_Expressed_Genes = {'Cover' , 'Percentile' , 'Number_of_genes_should_be_there'};
Header_Table_Gene_Never_Expressed  = {'Cover' , 'Percentile' , 'Number_of_genes_should_not_be_there'};

Table_Expressed_Genes.Properties.VariableNames = Header_Table_Expressed_Genes;
Table_Gene_Never_Expressed.Properties.VariableNames = Header_Table_Gene_Never_Expressed;
if printLevel==1
writetable(Table_Expressed_Genes, [dis,'\Table_Number_Genes_Should_be_there.txt']);
writetable(Table_Gene_Never_Expressed, [dis,'\Table_Number_Genes_Should_not_be_there.txt']);
end
end
