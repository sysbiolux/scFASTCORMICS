%% Discretization of the data : Cover and Percentile

mkdir ([path,'\Discretization_Step'])

Cover = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]/100; % Ex: Cover of 5 = 0.05 , of 20 = 0.2
%Cover = [0.0025 0.0075, 0.0125, 0.015, 0.025, 0.035, 0.045];
%Cover = [0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.02, 0.05, 0.1, 0.5];
percent = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90]; % Ex: Percent of 5 = 5 , of 20 = 20

kkk = 0;
kk = 0;

percentile_keep = zeros(3,1000);

% performing mapping takes ~1h30

% loading in the mapped_cluster_data

all_table = cell(files_length,1);

for i=1:files_length
    
    % First step : Load each data + add the correct header
    
    table_mapped = readtable([path,'\Mapped\Mapped_',scdataset(i).name]);
    
    all_table{i,1} = table_mapped;
    
    gene_id = table_mapped.Gene_Numeration;
    
    % Second step : Discretization starting with percentile followed by cover
    
    [I,ia,ib] = intersect(model_composite.genes,gene_id); % I = Intersect , ia = Index of model_genes and ib = index of gene annotation
    
    data = zeros(numel(model_composite.genes),(size(table_mapped,2)-1));
    gen = cell(numel(model_composite.genes),1);
    
    data(ia,:) = table2array(table_mapped(ib,2:end));
    gen(ia) = cellstr(gene_id(ib));
    
    data = array2table(data);
    data = table2array(data);
    
    mapping = zeros(numel(model_composite.rxns), size(data,2));
    
    for j=1:size(data,2)
        for ki=1:numel(model_composite.rxns)
            mapping(ki,j)= GPRrulesMapper_rFASTCORMICS(cell2mat(model_composite.rules(ki)),data(:,j));
        end
    end
    name = [path,'\Discretization_Step\mapping',num2str(i)];
    save(name,'mapping');
    
end


% takes about 5min  --  before taking the mapping loop out, the code would
% have ran for more than a day

mkdir ([dis,'\Unique_Core_Genes']);
mkdir ([dis,'\Unique_Core_Reactions']);
mkdir ([dis,'\Core_Genes']);
mkdir ([dis,'\Core_Reactions']);
mkdir ([dis,'\Core_Genes_Names_Model']);
mkdir ([dis,'\Core_Reactions_Each_Cluster']);
mkdir ([dis,'\Core_Reactions_Names_Model']);

Table_Core_Genes = [];
Table_Core_Reactions = [];
Table_Unique_Core_Genes = [];
Table_Unique_Core_Reactions = [];

for c=1:numel(Cover)
    for p=1:numel(percent)
        C_keep = []; % Core reactions
        Core_Genes_keep = []; % Core genes
        Unique_Core_Genes = []; % Unique core genes
        Unique_Core_Reactions = []; % Unique core reactions
        
        kkk = kkk + 1;
        kk = kk + 1;
        
        for i=1:files_length
            
            % First step : Load each data + add the correct header
            
            table_mapped = all_table{i,1};
            
            gene_id = table_mapped.Gene_Numeration;
            
            table_1 = table2array(table_mapped(:,2:end));
            table_1 = table_1(:);
            table_num = table_1;
            
            table_3_num = table_num(table_num > 0);
            
            Percentile = prctile(table_3_num(:),percent(p));
            percentile_keep(1,kkk) = Percentile;
            percentile_keep(2,kkk) = Cover(c);
            percentile_keep(3,kkk) = percent(p);
            
            name = ['\mapping',num2str(i),'.mat'];
            load([dis,name]);
            mapping(mapping < Percentile)= 0;
            
            [I,ia,ib] = intersect(model_composite.genes,gene_id); % I = Intersect , ia = Index of model_genes and ib = index of gene annotation
            
            data = zeros(numel(model_composite.genes),(size(table_mapped,2)-1));
            gen = cell(numel(model_composite.genes),1);
            
            data(ia,:) = table2array(table_mapped(ib,2:end));
            gen(ia) = cellstr(gene_id(ib));
            
            data = array2table(data);
            data = table2array(data);
            
            data(data < Percentile) = 0;
            
            C = find(sum(mapping>0,2) >= (size(data,2))*Cover(c));
            C_keep = [C_keep;C]; % C = index of the core reactions
            
            % Extracting the number of core genes of each threshold
            
            Core_Genes = find(sum(data>0,2) >= (size(data,2))*Cover(c));
            Core_Genes_keep = [Core_Genes_keep ; Core_Genes];
            Table_Core_Genes(kk,1) = Cover(c);
            Table_Core_Genes(kk,2) = percent(p);
            Table_Core_Genes(kk,3) = numel(Core_Genes_keep);
            
            % Creating a list with the names of the genes
            
            Core_Genes_Names = model_composite.genes(Core_Genes_keep);
            Core_Genes_Names = array2table(Core_Genes_Names);
            
            % Number of Unique core genes
            Core_Genes_Names_Unique = table2array(Core_Genes_Names);
            Core_Genes_Names_Unique = cellfun(@(S) S(1:end-2), Core_Genes_Names_Unique  , 'Uniform', 0);
            Core_Genes_Names_Unique = array2table(unique(Core_Genes_Names_Unique ));
            
            Unique_Core_Genes = unique([Unique_Core_Genes; Core_Genes_Names_Unique]);
            Table_Unique_Core_Genes(kk,1) = Cover(c);
            Table_Unique_Core_Genes(kk,2) = percent(p);
            Table_Unique_Core_Genes (kk,3) = numel(Unique_Core_Genes);
            
            % Number of Unique core reactions
            name_6 = [dis,'\Core_Reactions_Each_Cluster\Cluster_',num2str(i),'_Core_reactions_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat'];
            save (name_6, 'C')
            
            C_keep = unique(C_keep);
            
            Core_Reactions_Names_Unique = model_composite.rxns(C_keep);
            Core_Reactions_Names_Unique = cellfun(@(S) S(1:end-2), Core_Reactions_Names_Unique  , 'Uniform', 0);
            Core_Reactions_Names_Unique = array2table(unique(Core_Reactions_Names_Unique));
            
            Unique_Core_Reactions = unique([Unique_Core_Reactions ; Core_Reactions_Names_Unique]);
            Table_Unique_Core_Reactions(kk,1) = Cover(c);
            Table_Unique_Core_Reactions(kk,2) = percent(p);
            Table_Unique_Core_Reactions(kk,3) = numel(Unique_Core_Reactions);
            
            % Extracting the number of core reactions
            
            Table_Core_Reactions(kk,1) = Cover(c);
            Table_Core_Reactions(kk,2) = percent(p);
            Table_Core_Reactions(kk,3) = numel(C_keep);
            
            % Creating a list with the names of the reactions
            
            Core_Reactions_Names = model_composite.rxns(C_keep);
            Core_Reactions_Names = array2table(Core_Reactions_Names);
        end
        
        
        name_1 = ['\Unique_Core_Genes_Name_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
        Unique_Core_Genes_Name_List = [dis,'\Unique_Core_Genes',name_1,'.txt'];
        writetable(Core_Genes_Names_Unique,Unique_Core_Genes_Name_List);
        
        name_2 = ['\Core_Genes_Name_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
        Core_Genes_Name_List = [dis,'\Core_Genes_Names_Model',name_2,'.txt'];
        writetable(Core_Genes_Names,Core_Genes_Name_List);
        
        name_3 = [dis,'\Core_Genes\Core_genes_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat'];
        save(name_3,'Core_Genes')
        
        name_4 = ['\Unique_Core_Reaction_Name_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
        Unique_Core_Reaction_Name_List = [dis,'\Unique_Core_Reactions',name_4,'.txt'];
        writetable(Core_Reactions_Names_Unique,Unique_Core_Reaction_Name_List);
        
        name_5 = ['\Core_Reactions_Name_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p))];
        Core_Reaction_Name_List = [dis,'\Core_Reactions_Names_Model',name_5,'.txt'];
        writetable(Core_Reactions_Names,Core_Reaction_Name_List);
        
        name_7 = [dis,'\Core_Reactions\Core_Reactions_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat'];
        save(name_7,'C_keep');
        
    end
end


% Creating the tables containing the number of core reactions and core genes

Table_Core_Genes = array2table(Table_Core_Genes);
Table_Core_Reactions = array2table (Table_Core_Reactions);
Table_Unique_Core_Genes = array2table (Table_Unique_Core_Genes);
Table_Unique_Core_Reactions = array2table (Table_Unique_Core_Reactions);

Header_Genes = {'Cover' , 'Percentile' , 'Core_Genes'};
Header_Reactions = {'Cover' , 'Percentile' , 'Core_Reactions'};
Header_Unique_Genes = {'Cover' , 'Percentile' , 'Unique_Core_Genes'};
Header_Unique_Reactions = {'Cover' , 'Percentile' , 'Unique_Core_Reactions'};

Table_Core_Genes.Properties.VariableNames = Header_Genes;
Table_Core_Reactions.Properties.VariableNames = Header_Reactions;
Table_Unique_Core_Genes.Properties.VariableNames = Header_Unique_Genes;
Table_Unique_Core_Reactions.Properties.VariableNames = Header_Unique_Reactions;

writetable(Table_Core_Genes, [path,'\Discretization_Step\Table_Core_Genes.txt']);
writetable(Table_Core_Reactions, [path,'\Discretization_Step\Table_Core_Reactions.txt']);
writetable(Table_Unique_Core_Genes, [path,'\Discretization_Step\Table_Unique_Core_Genes.txt']);
writetable(Table_Unique_Core_Reactions, [path,'\Discretization_Step\Table_Unique_Core_Reactions.txt']);


%% Generating a context-specific model with the help of FASTCORE
% Takes a long time to run >7h
clear C_keep

mkdir ([dis,'\Extracting_all_possible_information']);
C_add = find(contains(model_composite.rxns,'biomassreaction'));
trex_ind = find(ismember(string(model_composite.subSystems),'Transport, extracellular'));
drug_ind = find(ismember(string(model_composite.subSystems),'Drug metabolism'));

for c=1:numel(Cover)
    for p=1:numel(percent)
        
        load([dis,'\Core_Reactions\Core_Reactions_Cover_' ,num2str(Cover(c)),'_Percentile_', num2str(percent(p)),'.mat'])
        C_keep = union(C_add, C_keep);
        C_keep = setdiff(C_keep,trex_ind);
        C_keep = setdiff(C_keep,drug_ind);
        A = fastcore_original(C_keep, model_composite, 1e-4);
        A_keep(A,1) = 1 ;
        
        model_composite.genes = cellstr(model_composite.genes);
        Context_Specific_model = removeRxns(model_composite,model_composite.rxns(setdiff(1:numel(model_composite.rxns),A)));
        Context_Specific_model.genes = cellstr(Context_Specific_model.genes);
        
        name = ['\Context_Specific_Model_','Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat'];
        name1 = ['\A_value_of_the_context_Specific_Model_','Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.mat'];
        save([dis,'\Extracting_all_possible_information',name],'Context_Specific_model');
        save([dis,'\Extracting_all_possible_information',name1],'A');
        
    end
end
