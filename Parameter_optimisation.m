function[input_data, Optimization_global]=Parameter_optimisation(Expanded_input_model,input_data,coverage_range, REI_range,path, ~, printLevel)
% Discretization of the data : Cover and REI

mkdir ([path,'\Discretization_Step'])
files_length= numel(input_data);

Nb_core=zeros(numel(coverage_range)*numel(REI_range));
for i=1:files_length
    % First step : Load each data + add the correct header
    table_mapped = input_data(i).table_final;
    gene_id = table_mapped.Gene_Numeration;
    % Second step : Discretization starting with REI followed by coverage
    [~,ia,ib] = intersect(Expanded_input_model.genes,gene_id); % I = Intersect , ia = Index of model_genes and ib = index of gene annotation
    data = zeros(numel(Expanded_input_model.genes),(size(table_mapped,2)-1));
    data(ia,:) = table2array(table_mapped(ib,2:end));
    mapping = zeros(numel(Expanded_input_model.rxns), size(data,2));
    [r,~]=find(sum(Expanded_input_model.rxnGeneMat,2)==1);
    %map genes expression to the reactions
    for ii=1:numel(r)
        [~,g]=find(Expanded_input_model.rxnGeneMat(r(ii),:));
        mapping(r(ii),:)=data(g,:);
    end
    match=find(~strcmp(Expanded_input_model.rules,'')& sum(Expanded_input_model.rxnGeneMat,2)>1);
    for j=1:size(data,2)
        for ki=1:numel(Expanded_input_model.rxns(match))
            mapping(match(ki),j)= GPRrulesMapper_rFASTCORMICS(cell2mat(Expanded_input_model.rules(match(ki))),data(:,j));
        end
    end
    input_data(i).mapping=mapping;
    if printLevel==1
        name = [path,'\Discretization_Step\mapping',num2str(i)];
        save(name,'mapping');
    end
end

dis = [path,'\Discretization_Step'];
mkdir ([dis,'\Unique_Core_Genes']);
mkdir ([dis,'\Unique_Core_Reactions']);
mkdir ([dis,'\Core_Genes']);
mkdir ([dis,'\Core_Reactions']);
mkdir ([dis,'\Core_Genes_Names_Model']);
mkdir ([dis,'\Core_Reactions_Each_Cluster']);
mkdir ([dis,'\Core_Reactions_Names_Model']);
Optimization_global=struct();
k=0;
C_keep=zeros(numel(Expanded_input_model.rxns),(numel(coverage_range)*numel(REI_range)));
A_keep=zeros(numel(Expanded_input_model.rxns),(numel(coverage_range)*numel(REI_range)));
Tresh=zeros(numel(coverage_range)*numel(REI_range),2);
%% Compute the core reaction set of each REI and cover setting
for c=1:numel(coverage_range)
    for p=1:numel(REI_range)
        k=k+1;
        
        for i=1:files_length
            % Load each data + add the correct header
            table_mapped = input_data(i).table_final;
            table_1 = table2array(table_mapped(:,2:end));
            table_1 = table_1(:);
            table_num = table_1;
            table_1_num = table_num(table_num > 0);
            
            REI = prctile(table_1_num(:),REI_range(p));
            
            mapping=input_data(i).mapping;
            mapping(mapping < REI)= 0;
            
            % Extracting the number of core genes of each threshold
            C = sum(mapping>0,2) >= (size(mapping,2)*coverage_range(c));
            C_keep(C,k)=1; % C = index of the core reactions
        end
    end
end
%% Generating a multi-cell population model with the help of FASTCORE

mkdir ([dis,'\Extracting_all_possible_information']);
C_add = find(contains(Expanded_input_model.rxns,'biomassreaction'));
trex_ind = find(ismember(string(Expanded_input_model.subSystems),'Transport, extracellular'));
drug_ind = find(ismember(string(Expanded_input_model.subSystems),'Drug metabolism'));
k=0;
Expanded_input_model.genes = cellstr(Expanded_input_model.genes);
%    Run Fastcore for each cover and REI setting
for c=1:numel(coverage_range)
    for p=1:numel(REI_range)
        k=k+1;
        
        C = union(C_add, find(C_keep(:,k)));
        C = setdiff(C,trex_ind);
        C = setdiff(C,drug_ind);
        A = fastcore_original(C, Expanded_input_model, 1e-4);
        A_keep(A,k) = 1 ;
        Tresh(k,1)=coverage_range(c);
        Tresh(k,2)=REI_range(p);
        C_all= unique(find(C_keep(:,k)));
        
        Core_Reactions_Names_Unique = Expanded_input_model.rxns(find(C_all));
        Core_Reactions_Names_Unique = cellfun(@(S) S(1:end-2), Core_Reactions_Names_Unique  , 'Uniform', 0);
        Nb_core(k,1)=numel(unique(Core_Reactions_Names_Unique));
    end
end

Optimization_global(1).A=A_keep;
Optimization_global(1).thresh=Tresh;
Optimization_global(1).C_keep=C_keep;
Optimization_global(1).Nb_core=Nb_core;

end