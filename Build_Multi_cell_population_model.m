function[input_data, Optimization_global]=Build_Multi_cell_population_model(model_composite,input_data,Cover_range, REI_range,path, ~, printLevel,function_keep)
% %% Discretization of the data : Cover and REI
%

mkdir ([path,'/Discretization_Step'])
files_length= numel(input_data);
% performing mapping takes ~1h30
% loading in the mapped_cluster_data
Nb_core=zeros(numel(Cover_range)*numel(REI_range));
disp('... mapping ...')
tic
for i=1:files_length
    i
    % First step : Load each data + add the correct header
    table_mapped = input_data(i).table_final;
    gene_id = table_mapped.Gene_Numeration;
    % Second step : Discretization starting with percentile followed by cover
    [~,ia,ib] = intersect(model_composite.genes,gene_id); % I = Intersect , ia = Index of model_genes and ib = index of gene annotation
    data = zeros(numel(model_composite.genes),(size(table_mapped,2)-1));
try %old input data format
    data(ia,:) = table2array(table_mapped(ib,2:end)); %%% isnumeric(table2array(table_mapped(ib,2:end)))??
catch %new input data format
    data(ia,:) = str2double(table_mapped{ib,2:end});
end
    mapping = zeros(numel(model_composite.rxns), size(data,2));
    [r,~]=find(sum(model_composite.rxnGeneMat,2)==1);
    for ii=1:numel(r)
        [~,g]=find(model_composite.rxnGeneMat(r(ii),:));
        mapping(r(ii),:)=data(g,:);
    end
    match=find(~strcmp(model_composite.rules,'')& sum(model_composite.rxnGeneMat,2)>1);

%     for j=1:size(data,2)
%         for ki=1:numel(model_composite.rxns(match))
%             mapping(match(ki),j)= GPRrulesMapper_rFASTCORMICS(cell2mat(model_composite.rules(match(ki))),data(:,j));
%         end
%     end


rules = regexprep(model_composite.rules,'x\(([0-9]*)\)','x($1,:)');
for ki=1:numel(model_composite.rxns(match))
            mapping(match(ki),:)= GPRrulesMapper_rFASTCORMICS(cell2mat(rules(match(ki))),...
                                                              data);
end


% % % load(['C:\Users\thomas.sauter\OneDrive - University of Luxembourg\work_other\Projects\Elena_2024\scMetMod\1ResultsData_1_model_orig\Discretization_Step\mapping' num2str(i) '.mat'])

input_data(i).mapping=mapping;
    if printLevel==1
        name = [path,'/Discretization_Step/mapping',num2str(i)];
        save(name,'mapping');
    end
end
toc

% takes about 5min  --  before taking the mapping loop out, the code would
% have ran for more than a day
dis = [path,'/Discretization_Step'];

mkdir ([dis,'/Unique_Core_Genes']);
mkdir ([dis,'/Unique_Core_Reactions']);
mkdir ([dis,'/Core_Genes']);
mkdir ([dis,'/Core_Reactions']);
mkdir ([dis,'/Core_Genes_Names_Model']);
mkdir ([dis,'/Core_Reactions_Each_Cluster']);
mkdir ([dis,'/Core_Reactions_Names_Model']);
Optimization_global=struct();
k=0;
C_keep=zeros(numel(model_composite.rxns),(numel(Cover_range)*numel(REI_range)));
A_keep=zeros(numel(model_composite.rxns),(numel(Cover_range)*numel(REI_range)));
Tresh=zeros(numel(Cover_range)*numel(REI_range),2);

for c=1:numel(Cover_range)
    for p=1:numel(REI_range)
        k=k+1;
        
        for i=1:files_length
            
            % First step : Load each data + add the correct header
            table_mapped = input_data(i).table_final;

            if isnumeric(table2array(table_mapped(ib,2:end)))
                table_1 = table2array(table_mapped(:,2:end)); %%%
            else
                table_1 = str2double(table_mapped{:,2:end});
            end
            
            table_1 = table_1(:);
            table_num = table_1;
            
            table_1_num = table_num(table_num > 0);
            
            Percentile = prctile(table_1_num(:),REI_range(p));
            
            mapping=input_data(i).mapping;
            mapping(mapping < Percentile)= 0;
            
            
            
            
            C = sum(mapping>0,2) >= (size(mapping,2)*Cover_range(c))*0.9;
            C_keep(C,k)=1; % C = index of the core reactions
            
            % Extracting the number of core genes of each threshold
            
            
            
        end
    end
    
end
%% Generating a context-specific model with the help of FASTCORE
% Takes a long time to run >7h

mkdir ([dis,'/Extracting_all_possible_information']);
C_add = find(contains(model_composite.rxns,strrep(function_keep.obj,'_','')));
trex_ind = find(ismember(string(model_composite.subSystems),'Transport, extracellular'));
drug_ind = find(ismember(string(model_composite.subSystems),'Drug metabolism'));
k=0;
model_composite.genes = cellstr(model_composite.genes);

for c=1:numel(Cover_range)
    for p=1:numel(REI_range)
        k=k+1;
        C = union(C_add, find(C_keep(:,k)));
        C = setdiff(C,trex_ind);
        C = setdiff(C,drug_ind);
        A = fastcore_original(C, model_composite, 1e-4);
        A_keep(A,k) = 1 ;
        Tresh(k,1)=Cover_range(c);
        Tresh(k,2)=REI_range(p);
        C_all= unique(find(C_keep(:,k)));
        
        Core_Reactions_Names_Unique = model_composite.rxns(find(C_all));
        Core_Reactions_Names_Unique = cellfun(@(S) S(1:end-2), Core_Reactions_Names_Unique  , 'Uniform', 0);
        Nb_core(k,1)=numel(unique(Core_Reactions_Names_Unique));
    end
end
Optimization_global(1).A=A_keep;
Optimization_global(1).thresh=Tresh;
Optimization_global(1).C_keep=C_keep;
Optimization_global(1).Nb_core=Nb_core;

end
