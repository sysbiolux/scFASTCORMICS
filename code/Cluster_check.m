%% Counting cells for cluster cell ratio 
table1 = readtable('D:\scFASTCORMICS\Datasets\Data_1\scdata_cluster_0.txt');
table2 = readtable('D:\scFASTCORMICS\Datasets\Data_1\scdata_cluster_1.txt');
table3 = readtable('D:\scFASTCORMICS\Datasets\Data_1\scdata_cluster_2.txt');
table4 = readtable('D:\scFASTCORMICS\Datasets\Data_1\scdata_cluster_3.txt');
table5 = readtable('D:\scFASTCORMICS\Datasets\Data_1\scdata_cluster_4.txt');

cluster_names = table1.Properties.VariableNames(4:end);
count = sum(contains(cluster_names,'Fibroblast'));

res = 100*[83, 12, 11, 9, 6, 4, 3,]/128;



%% Similarity matrix of model metrics
% Loading the model and construct a table containing metrics for each
% cluster (metrics that do not belong to a cluster are excluded)
load('D:\scFASTCORMICS\Results_Data_1_model\Discretization_Step\Extracting_all_possible_information\Context_Specific_Model_Cover_0.01_Percentile_60.mat');
dico = readtable('D:\scFASTCORMICS\dico_biodbnet.txt');

% Reactions similarity matrix - Cluster exclusive reactions
rxn_list = {};
trans_ind = find(contains(Context_Specific_model.rxns, 'Trans_'));
rx_clust_lim = trans_ind(1)-1;
rx = Context_Specific_model.rxns(1:rx_clust_lim);
num_clust = max(str2double(extractAfter(rx,'_')));
for i=1:num_clust
    rx_num = erase(rx(find(contains(rx,['_',num2str(i)]))),['_',num2str(i)]);
    rxn_list(i) = {rx_num};
end
rxn_sim_matrix = zeros(num_clust);
for j=1:num_clust
    cycle1 = rxn_list{j};
    for k=1:num_clust
        cycle2 = rxn_list{k};
        same_rxns = numel(intersect(cycle1,cycle2));
        rxn_sim_matrix(j,k) = same_rxns;
    end
end

rxns_exclusive = {};
for i=1:num_clust
    base = rxn_list{i};
    for j=1:num_clust
        subj = rxn_list{j};
        diff = setdiff(base,subj);
        if numel(diff) ~= 0
            base = diff;
        end
    end
    rxns_exclusive{i} = base;
end


% Metabolites similarity matrix
mets_list = {};
mets_clust_lim = 15909;
met = Context_Specific_model.mets(1:mets_clust_lim);
for i=1:num_clust
    met_num = erase(met(find(contains(met,['_',num2str(i)]))),['_',num2str(i)]);
    mets_list(i) = {met_num};
end
met_sim_matrix = zeros(num_clust);
for j=1:num_clust
    cycle1 = mets_list{j};
    for k=1:num_clust
        cycle2 = mets_list{k};
        same_mets = numel(intersect(cycle1,cycle2));
        met_sim_matrix(j,k) = same_mets;
    end
end

mets_exclusive = {};
for i=1:num_clust
    base = mets_list{i};
    for j=1:num_clust
        subj = mets_list{j};
        diff = setdiff(base,subj);
        if numel(diff) ~= 0
            base = diff;
        end
    end
    mets_exclusive{i} = base;
end


% Genes similarity matrix
gl = table2array(readtable('D:\scFASTCORMICS\Results_Data_1_model\Discretization_Step\Core_Genes_Names_Model\Core_Genes_Name_Cover_0.01_Percentile_60.txt'));
genes_list = {};
for i=1:num_clust
    gen_num = gl(find(contains(string(gl(:,2)),num2str(i))),1);
    genes_list(i) = {gen_num};
end
gen_sim_matrix = zeros(num_clust);
for j=1:num_clust
    cycle1 = genes_list{j};
    for k=1:num_clust
        cycle2 = genes_list{k};
        same_genes = numel(intersect(cycle1,cycle2));
        gen_sim_matrix(j,k) = same_genes;
    end
end

genes_exclusive = {};
for i=1:num_clust
    base = genes_list{i};
    for j=1:num_clust
        subj = genes_list{j};
        diff = setdiff(base,subj);
        if numel(diff) ~= 0
            base = diff;
        end
    end
    genes_exclusive{i} = base;
end

for i=1:num_clust
    symbol_list{i} = dico(find(ismember(dico{:,1},genes_exclusive{i})),2);
end


%% Gene exclusive in cluster c3 and c4  + the intersection with c0
c0c3_genes = intersect(genes_list{1},genes_list{4});
c0c4_genes = intersect(genes_list{1},genes_list{5});

c0c3_fin = [c0c3_genes; genes_exclusive{4}];
c0c4_fin = [c0c4_genes; genes_exclusive{5}];

c0c3_symbol = dico(find(ismember(dico{:,1}, c0c3_fin)),2);
c0c4_symbol = dico(find(ismember(dico{:,1}, c0c4_fin)),2);

%% Discretization with thresholds without the composite model constriction
Cover = 0.01;
percent = 60;

gene_table = readtable('D:\scFASTCORMICS\Datasets\Data_1\scdata_cluster_4.txt');
gene_mat = table2array(gene_table(:,4:end));
gene_mat = gene_mat(:);
gene_num = gene_mat(gene_mat > 0);
Percentile = prctile(gene_num(:),percent);
gene_mat_per = table2array(gene_table(:,4:end));
gene_mat_per(gene_mat_per < Percentile) = 0;
gene_mat_cover = find(sum(gene_mat_per>0,2) >= (size(gene_mat_per,2))*Cover);
gene_table_clear = gene_table(gene_mat_cover,:);
gene_selec = gene_table_clear(:,3);
gene_selec = extractBefore(gene_selec{:,1},'_ENSG');
gene_selec4 = extractAfter(gene_selec,'_');

dis_genes_list = {gene_selec0, gene_selec1, gene_selec2, gene_selec3, gene_selec4};

genes_exclusive = {};
for i=1:numel(dis_genes_list)
    base = dis_genes_list{i};
    for j=1:numel(dis_genes_list)
        subj = dis_genes_list{j};
        diff = setdiff(base,subj);
        if numel(diff) ~= 0
            base = diff;
        end
    end
    genes_exclusive{i} = base;
end
