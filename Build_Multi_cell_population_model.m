function[input_data, Optimization_global]=Build_Multi_cell_population_model(model_composite,input_data,Cover_range, REI_range,path, ~, printLevel,function_keep)
% Build_Multi_cell_population_model
%   This function takes the generic consistent expanded model and performs
%   3 main tasks: 
%       1) Mapping of the gene expression to protein abundance using the
%       GPR rules 
%       2) Determines the set of core reactions based on a threshold set (coverage and REI) on
%       the mapped protein abundance
%       3) Builds consistent context-specific multi-cell population model for every set of core reaction by running fastcore
% 
%   INPUT: 
%   model_composite:    consistent generic multi population model 
%   input_data:         struct with two fields, in this function the
%                       input_data field is used, which entails one table
%                       for each cluster/cell population in the model. The
%                       first column entails the gene names with a postfix
%                       depending on the cluster number and the rest of the
%                       columns are the single cell expression data
%   Cover_range:        range of different values of coverage for which the 
%                       set of core reaction should be defined  
%   REI_range:          range of different values of REI for which the set 
%                       of core reaction should be defined
%   
%
%   OUTPUT:
%   input_data:         input data, unchanged
%   Optimization_global:final consistent, context-specific multi-population
%                       model
%   
% Discretization of the data : Cover and REI
%

mkdir ([path,'/Discretization_Step'])
files_length= numel(input_data);
% performing mapping takes ~1h30
% loading in the mapped_cluster_data
Nb_core=zeros(numel(Cover_range)*numel(REI_range));
disp('... mapping ...')
tic
% which of the rxns rules are not empyt ? therefore need to be mapped ? 
match=find(~strcmp(model_composite.rules,''));
% change the rules of the rxns, in order to be able to loop through the
% rxns for all cells/columns in one go x(1,:) | x(2,:) instead of x(1)|x(2)
% the rules are adjusted here, this adjustment helps to be able to loop
% over each of the cluster data in one go, instead of per cell, since
% we are computing the rules on the whole matrix in one go instead of
% every single cell separatly!
rules = regexprep(model_composite.rules,'x\(([0-9]*)\)','x($1,:)');


concat_data = [];
files_length = 60;
for cluster=1:files_length
    disp("cluster nr processed:")
    cluster
    % First step : Load each data + add the correct header
    table_mapped = input_data(1).table_final;
    gene_id = table_mapped.Gene_Numeration;
    % Second step : Discretization starting with percentile followed by cover
    [~,ia,ib] = intersect(model_composite.genes,gene_id); % I = Intersect , ia = Index of model_genes and ib = index of gene annotation
    data = zeros(numel(model_composite.genes),(size(table_mapped,2)-1));
try %old input data format
    data(ia,:) = table2array(table_mapped(ib,2:end)); %%% isnumeric(table2array(table_mapped(ib,2:end)))??
catch %new input data format
    data(ia,:) = str2double(table_mapped{ib,2:end});
end
    input_data(cluster).data = [data,data,data,data,data,data,data,data,data,data];
    input_data(cluster).mapping = zeros(numel(model_composite.rxns), size(input_data(cluster).data,2));
    input_data(cluster).store_column_nr_of_clustercells_in_concat_object = size(concat_data,2)+1:size(concat_data,2)+ size(input_data(cluster).data,2);
    concat_data = [concat_data,input_data(cluster).data];
end
size(concat_data)
mapping = zeros(numel(model_composite.rxns), size(concat_data,2));
for ki=1:numel(model_composite.rxns(match))
            mapping(match(ki),:)= GPRrulesMapper_rFASTCORMICS(cell2mat(rules(match(ki))),...
                                                                  concat_data);
end


% 
% for cluster=1:files_length
% %     [r,~]=find(sum(model_composite.rxnGeneMat,2)==1);
% %     tic
% %     for ii=1:numel(r)
% %         [~,g]=find(model_composite.rxnGeneMat(r(ii),:));
% %         mapping(r(ii),:)=data(g,:);
% %     end
% %     toc
% %     match=find(~strcmp(model_composite.rules,'')& sum(model_composite.rxnGeneMat,2)>1);
% %     for j=1:size(data,2)
% %         for ki=1:numel(model_composite.rxns(match))
% %             mapping(match(ki),j)= GPRrulesMapper_rFASTCORMICS(cell2mat(model_composite.rules(match(ki))),data(:,j));
% %         end
% %     end
% %     input_data(i).data = data;
% end
files_length
for cluster=1:files_length
    % mapping back of the global mapping performed on the concatenated
    % clusters - back to the individual inputdata.(clusternr) objects!
    input_data(cluster).mapping=mapping(:,input_data(cluster).store_column_nr_of_clustercells_in_concat_object);
%     if printLevel==1
%         name = [path,'/Discretization_Step/mapping',num2str(cluster)];
%         save(name,'mapping');
%     end
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
