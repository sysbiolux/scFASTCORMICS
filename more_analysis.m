%%
load('D:\scFASTCORMICS\Results_Data_1_model_Clow\Discretization_Step\Extracting_all_possible_information\Context_Specific_Model_Cover_0.0075_Percentile_60.mat')
cluster_num = 5;
% [m,~] = find(model_composite.S(:,10));
% meta = model_composite.mets(m);
% unique(strtok(meta,'['))
% extractAfter(meta,']_')

%% General model analysis
subsys = Context_Specific_model.subSystems;
subsys = string(subsys);
subsys_uni = unique(subsys); % "" empty cell, rxns that do not belong to a pathway (yet to determine)
[subsys_count, subsys_names] = histcounts(categorical(subsys));
table_subsys = table(subsys_names', subsys_count');
table_subsys.Properties.VariableNames = {'subsys_names','subsys_count'};
table_subsys = sortrows(table_subsys,{'subsys_count'},{'descend'});
%writetable(table_subsys,'.\Context_Subsys_table.txt')

% Bar graph of subsystem count
[~,I] = sort(table_subsys.subsys_count,'descend'); %sort by descreasing subsystems counts
figure
bar(table_subsys.subsys_count(I));
xlim([1 numel(table_subsys.subsys_names)]);
% xticklabels(T.pathway_names(I)) %enable this to see the names, the figure will be very small
xticks(1:numel(table_subsys.subsys_names));
xticklabels(table_subsys.subsys_names);
xtickangle(60);
title('Number of reactions per subsystem');

%% Reactions per cluster 
total_rxns = numel(Context_Specific_model.rxns);
cycle_rxns = Context_Specific_model.rxns(1:19909);
cluster_rxns = {};
for i=1:cluster_num
    cluster_rxns{i} = cycle_rxns(find(contains(cycle_rxns, ['_',num2str(i)])));
end


%% Reactions metabolites analysis
% load('D:\scFASTCORMICS\Results_Data_1_model\Discretization_Step\Extracting_all_possible_information\Context_Specific_Model_Cover_0.01_Percentile_60.mat');
% load('D:\scFASTCORMICS\Results_Data_1_model\composite_model.mat');

ind = find(contains(Context_Specific_model.mets,'[u]'));
model_mets = Context_Specific_model.mets(ind);
[~,r] = find(Context_Specific_model.S(ind,:));
model_rxns = Context_Specific_model.rxns(r);
ind_trans = find(contains(Context_Specific_model.rxns(r),'Trans_'));
rxns_need = setdiff(r,r(ind_trans));
% printRxnFormula(Context_Specific_model,Context_Specific_model.rxns(rxns_need));

% find only transport reactions that transitions between compartments
% check in S matrix for reactions that only contain [u] metabolites
u_only = Context_Specific_model.S(:,rxns_need);
rxns_need2 = [];
for i=1:size(u_only,2)
   row_mets = Context_Specific_model.mets(find(u_only(:,i)));
   row_mets_count = numel(row_mets);
   row_u_count = sum(contains(row_mets,'[u]'));
   if row_mets_count ~= row_u_count
       rxns_need2(i,1) = rxns_need(i);
   else
       rxns_need2(i,1) = 0;
   end
end
% remove zeros
rxns_need2 = rxns_need2(rxns_need2~=0);
table_rxnformula_u = [Context_Specific_model.subSystems(rxns_need2),Context_Specific_model.rxns(rxns_need2),printRxnFormula(Context_Specific_model,Context_Specific_model.rxns(rxns_need2))];

% pathway count
cmodel_subsys = Context_Specific_model.subSystems(rxns_need2);
cmodel_rxnstrans = Context_Specific_model.rxns(rxns_need2);


%separate into clusters
sum(contains(Context_Specific_model.rxns(rxns_need2),'_1'))
clust1 = find(contains(Context_Specific_model.rxns(rxns_need2),'_1'));
clust1_ind = find(ismember(Context_Specific_model.rxns,cmodel_rxnstrans(clust1)));
clust1_subsys = Context_Specific_model.subSystems(clust1_ind);

sum(contains(Context_Specific_model.rxns(rxns_need2),'_2'))
clust2 = find(contains(Context_Specific_model.rxns(rxns_need2),'_2'));
clust2_ind = find(ismember(Context_Specific_model.rxns,cmodel_rxnstrans(clust2)));
clust2_subsys = Context_Specific_model.subSystems(clust2_ind);

sum(contains(Context_Specific_model.rxns(rxns_need2),'_3'))
clust3 = find(contains(Context_Specific_model.rxns(rxns_need2),'_3'));
clust3_ind = find(ismember(Context_Specific_model.rxns,cmodel_rxnstrans(clust3)));
clust3_subsys = Context_Specific_model.subSystems(clust3_ind);

sum(contains(Context_Specific_model.rxns(rxns_need2),'_4'))
clust4 = find(contains(Context_Specific_model.rxns(rxns_need2),'_4'));
clust4_ind = find(ismember(Context_Specific_model.rxns,cmodel_rxnstrans(clust4)));
clust4_subsys = Context_Specific_model.subSystems(clust4_ind);

sum(contains(Context_Specific_model.rxns(rxns_need2),'_5'))
clust5 = find(contains(Context_Specific_model.rxns(rxns_need2),'_5'));
clust5_ind = find(ismember(Context_Specific_model.rxns,cmodel_rxnstrans(clust5)));
clust5_subsys = Context_Specific_model.subSystems(clust5_ind);

% List of metabolites present in the reactions
[m,~] = find(Context_Specific_model.S(:,rxns_need2));
meta = Context_Specific_model.mets(m);


%% Reaction pathway analysis

%subsys = Context_Specific_model.subSystems;
subsys = Context_Specific_model.subSystems(rxns_need2);
subsys = model.subSystems;
subsys = string(subsys);
subsys_uni = unique(subsys); % "" empty cell, rxns that do not belong to a pathway (yet to determine)
[subsys_count, subsys_names] = histcounts(categorical(subsys));
table_subsys = table(subsys_names', subsys_count');
table_subsys.Properties.VariableNames = {'subsys_names','subsys_count'};
table_subsys = sortrows(table_subsys,{'subsys_count'},{'descend'});
writetable(table_subsys,'.\Context_Subsys_table.txt')

% Bar graph of subsystem count
[~,I] = sort(table_subsys.subsys_count,'descend'); %sort by descreasing subsystems counts
figure
bar(table_subsys.subsys_count(I));
xlim([1 numel(table_subsys.subsys_names)]);
% xticklabels(T.pathway_names(I)) %enable this to see the names, the figure will be very small
xticks(1:numel(table_subsys.subsys_names));
xticklabels(table_subsys.subsys_names);
xtickangle(60);
title('Number of reactions per subsystem');


mets_comp = extractAfter(Context_Specific_model.mets,'[');
mets_comp = extractBefore(mets_comp,']');
% mets_comp = extractAfter(mets_comp,'[');
mets_clust1 = sum(contains(Context_Specific_model.mets,'_1'));
mets_clust2 = sum(contains(Context_Specific_model.mets,'_2'));
mets_clust3 = sum(contains(Context_Specific_model.mets,'_3'));
mets_clust4 = sum(contains(Context_Specific_model.mets,'_4'));
mets_clust5 = sum(contains(Context_Specific_model.mets,'_5'));

mets_c_count = sum(contains(mets_comp,'c')); %
mets_e_count = sum(contains(mets_comp,'e')); %
mets_g_count = sum(contains(mets_comp,'g')); %
mets_i_count = sum(contains(mets_comp,'i')); %
mets_l_count = sum(contains(mets_comp,'l')); %
mets_m_count = sum(contains(mets_comp,'m')); %
mets_n_count = sum(contains(mets_comp,'n')); %
mets_r_count = sum(contains(mets_comp,'r')); %
mets_u_count = sum(contains(mets_comp,'u')); %
mets_x_count = sum(contains(mets_comp,'x')); %


%% Optimization 
% load model
% load('D:\scFASTCORMICS\Results_Data_1_model\Discretization_Step\Extracting_all_possible_information\Context_Specific_Model_Cover_0.01_Percentile_60.mat');

% check if biomass functions exist
checkObjective(Context_Specific_model);
cmodel_rxns = Context_Specific_model.rxns;

% check if biomass reactions exist
bm = [string(cmodel_rxns(find(contains(cmodel_rxns,'biomassreaction')))), find(contains(cmodel_rxns,'biomassreaction'))];

% Partial weight of clusters in the context specific model
total_cells = 128+106+18+17+10;
clust_cells = [128, 106, 18, 17, 10];
clust_ratio = clust_cells/total_cells;
sum(clust_ratio)
Context_Specific_model.c(ismember(Context_Specific_model.rxns, 'biomassreaction_1'))=1;
Context_Specific_model.c(ismember(Context_Specific_model.rxns, 'biomassreaction_2'))=1;
Context_Specific_model.c(ismember(Context_Specific_model.rxns, 'biomassreaction_3'))=1;
Context_Specific_model.c(ismember(Context_Specific_model.rxns, 'biomassreaction_4'))=1;
Context_Specific_model.c(ismember(Context_Specific_model.rxns, 'biomassreaction_5'))=1;
Context_Specific_model.c(find(Context_Specific_model.c))
cind = find(Context_Specific_model.c);
Context_Specific_model.c(cind) = clust_ratio;
checkObjective(Context_Specific_model);


changeCobraSolver('ibm_cplex');

% FVA 
tic
[minFlux,maxFlux] = fluxVariability(Context_Specific_model,100,'max',cmodel_rxnstrans);
toc

table_flux = [cmodel_rxnstrans printRxnFormula(Context_Specific_model, Context_Specific_model.rxns(rxns_need2)) num2cell(full(minFlux)) num2cell(full(maxFlux)) Context_Specific_model.subSystems(rxns_need2)];

biomassreactions = [cmodel_rxnstrans printRxnFormula(Context_Specific_model, Context_Specific_model.rxns(rxns_need2)) num2cell(full(minFlux)) num2cell(full(maxFlux)) Context_Specific_model.subSystems(rxns_need2)];


%remove no flux reactions
nonzeroflux_ind = unique([find(minFlux);find(maxFlux)]);
%nonzeroflux_ind = unique([find(rminFlux);find(rmaxFlux)]);
biomassreaction_nonzero = biomassreactions(nonzeroflux_ind,:);

biomassreaction1 = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),'_1')),:);
biomassreaction2 = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),'_2')),:);
biomassreaction3 = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),'_3')),:);
biomassreaction4 = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),'_4')),:);
biomassreaction5 = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),'_5')),:);

% determine outflowing reactions
smatrix_nonzero = Context_Specific_model.S(:,find(ismember(Context_Specific_model.rxns,biomassreaction_nonzero(:,1))));
[m_nonzero,~] = find(smatrix_nonzero==1);
%all_mets = Context_Specific_model.mets(m_nonzero);


% USE SCRIPT metsbolites_heatmap INSTEAD OF THIS PART

% all_flowout_mets = {};
% all_flowin_mets = {};
% for i=1:cluster_num
%     bm_nonzero_clust = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),['_',num2str(i)])),1);
%     smatrix_nonzero = Context_Specific_model.S(:,find(ismember(Context_Specific_model.rxns,bm_nonzero_clust)));
%     [m_out,~] = find(smatrix_nonzero == 1);
%     [m_in,~] = find(smatrix_nonzero == -1);
%     csm_mets_out = Context_Specific_model.mets(m_out);
%     csm_mets_in = Context_Specific_model.mets(m_in);
%     u_mets_out = csm_mets_out(find(contains(csm_mets_out,'[u]')));
%     u_mets_in = csm_mets_in(find(contains(csm_mets_in,'[u]')));
%     all_flowout_mets{i} = u_mets_out;
%     all_flowin_mets{i} = u_mets_in;
% end
% 
% % the intersect defines the interactions between clusters
% % rows: mets exiting cluster; column: entering cluster
% for i=1:cluster_num
%     for j=1:cluster_num
%         mets_int_list = extractBefore(intersect(all_flowout_mets{i}, all_flowin_mets{j}),'[u]');
%         mets_full_ind = find(ismember(extractBefore(model.mets,'['),mets_int_list));
%         mets_fullname = unique(model.metNames(mets_full_ind));
%         mets_cellmat{i,j} = mets_fullname;
%     end
% end
% 
% 
% % how many metabolites are involved in a cluster for all the reactions
% met_count_all = {};
% for indice=1:cluster_num
%     metrxns_clust_ind = find(contains(Context_Specific_model.rxns,['_',num2str(indice)]));
%     [met_ind,~] = find(Context_Specific_model.S(:,metrxns_clust_ind));
%     met_count_all{indice} = unique(extractBefore(Context_Specific_model.mets(met_ind),'['));
% end
% 
% % how many metabolites are involved in a cluster's active transport reactions
% met_count_transport = {};
% for indice=1:cluster_num
%     metrxns_clust = biomassreaction_nonzero(find(contains(biomassreaction_nonzero(:,1),['_',num2str(indice)])),1);
%     [met_ind,~] = find(Context_Specific_model.S(:,find(ismember(Context_Specific_model.rxns,metrxns_clust))));
%     met_count_transport{indice} = unique(extractBefore(Context_Specific_model.mets(met_ind),'['));
% end
% 
% all_mets = {};
% for i=1:5
%     for j=1:5
%         all_mets = [all_mets;mets_cellmat{i,j}];
%     end
% end
% all_mets = unique(all_mets);
% mets_uni_table = unique(extractBefore(model.mets(find(ismember(model.metNames,all_mets))),'['));
% mets_uni_table = append(mets_uni_table,'[u]');
% 
% mets_in_table=[];
% mets_out_table=[];
% for i=1:numel(mets_uni_table)
%     for j=1:5
%         ind_in = ismember(mets_uni_table(i),all_flowin_mets{j});
%         ind_out = ismember(mets_uni_table(i),all_flowout_mets{j});
%         if ind_in>0
%             mets_in_table(i,j)=1;
%         else mets_in_table(i,j)=0;
%         end
%         if ind_out>0
%             mets_out_table(i,j)=-1;
%         else mets_out_table(i,j)=0;
%         end
%         if mets_out_table(i,j)==0 && mets_in_table(i,j)==0
%             mets_out_table(i,j)=-2;
%             mets_in_table(i,j)=-2;
%         end
%     end
% end
% mets_tot_table = mets_in_table + mets_out_table;
% headers = {'Epithelial 1','Epithelial 2','Epithelial 3','T-cell','B-cell'}';
% removelist = [17;27;30;31;42;44];
% mets_tot_table(removelist,:) = [];
% mets_uni_table(removelist) = [];
% mets_uni_table = extractBefore(mets_uni_table,'[');
% mets_uni_table = strrep(mets_uni_table,'_',' ');
% 
% figure('WindowState','maximized')
% I = clustergram(mets_tot_table);
% %caxis([-1 1])
% title('Metabolite inter-cluster flux')
% ylabel('Metabolites','FontSize',16)
% xlabel('Clusters','FontSize',16)
% xticks([1:5]); 
% xticklabels({'Epithelial 1','Epithelial 2','Epithelial 3','T-cell','B-cell'})
% xtickangle(0)
% yticks([1:50])
% yticklabels(mets_uni_table)
% ytickangle(30)
% colorbar;
% ax = gca;
% c = ax.FontSize;
% ax.FontSize = 11;
% 
% defaultPosition = get(0,'DefaultFigurePosition');
% screensize = get(groot, 'Screensize');
% screensize = [1,1,800,800];
% set(0, 'DefaultFigurePosition', screensize);
%               
% cgo_J = clustergram(mets_tot_table,...
%     'RowLabels', mets_uni_table,...
%     'RowLabelsRotate', 360,...
%     'ColumnLabels', headers,...
%     'ColumnLabelsRotate', 360,...
%     'Cluster', 'Column', ...
%     'symmetric','False',...
%     'Colormap', redbluecmap)
% addTitle(cgo_J,{'Metabolite inter-cluster flux'});
% set(0,'ShowHiddenHandles','on')
% allhnds = get(0,'Children');
% h = findall(allhnds, 'Tag', 'HeatMapAxes');
% set(h, 'FontSize', 9)
% plot(cgo_J)
% 
% %Obtain list of full metabolite names
% m_list = mets_cellmat{1,4};
% m_list = extractBefore(m_list,'[u]');
% r3d_m = extractBefore(model.mets,'[');
% r3d_mfull = model.metNames;
% 
% m_r3d_ind = find(ismember(r3d_m,m_list));
% m_r3dmfull = r3d_mfull(m_r3d_ind);
% mfull_table = cell2table([r3d_m(m_r3d_ind),model.metNames(m_r3d_ind)]);
% mfull_uni_table = unique(mfull_table,'stable');
% 

%% Analysis from AdvSys Course Day 4

% model_keep_generic = zeros(numel(model_composite.rxns),2);
% model_keep_generic(:,1) = 1;
% rxns_present = ismember(model_composite.rxns,Context_Specific_model.rxns);
% model_keep_generic(:,2) = rxns_present;
% 
% figure
% plot(model_keep_generic);
% 
% title('reaction presence')
% ylabel('reactions')
% colorbar
% xticks(1:2)
% xticklabels({'cancer','healthy'})


% Reaction count for each pathway
Pathways = table(unique(string(Context_Specific_model.subSystems)));
[pathways, ~, ub] = unique(string(Context_Specific_model.subSystems));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);

[I, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.ContextSpecific(ia) = T.path_counts(ib) ;

% Bar graph of subsystem count in context specific model
[~,I] = sort(Pathways.ContextSpecific,'descend'); %sort by descreasing subsystems counts
figure
bar(Pathways.ContextSpecific(I));
xlim([1 numel(Pathways.ContextSpecific)]);
% xticklabels(T.pathway_names(I1)) %enable this to see the names, the figure will be very small
xticks(1:numel(Pathways.ContextSpecific));
xticklabels(Pathways.Var1(I));
xtickangle(60);
ylabel('Reactions')
title('Context Specific model: Number of reactions per subsystem');


% Reaction count of each pathway in generic model
Pathways_composite = table(unique(string(model_composite.subSystems)));
[pathways_composite, ~, ub1] = unique(string(model_composite.subSystems));
path_counts_composite = histc(ub1, 1:length(pathways_composite));
T = table(pathways_composite, path_counts_composite);

[I, ia, ib] = intersect(Pathways_composite.Var1, T.pathways_composite);
Pathways_composite.composite(ia) = T.path_counts_composite(ib) ;

% Bar graph of subsystem count in context specific model
[~,I1] = sort(Pathways_composite.composite,'descend'); %sort by descreasing subsystems counts
figure
bar(Pathways_composite.composite(I1));
xlim([1 numel(Pathways_composite.composite)]);
% xticklabels(T.pathway_names(I1)) %enable this to see the names, the figure will be very small
xticks(1:numel(Pathways_composite.composite));
xticklabels(Pathways_composite.Var1(I1));
xtickangle(60);
ylabel('Reactions')
title('Generic model: Number of reactions per subsystem');

% compare models' pathways
figure
hold on
bar(Pathways.ContextSpecific(I1),'FaceAlpha',1)
bar(Pathways_composite.composite(I1) ,'FaceAlpha',0.7)
legend('Cancer','Healthy')
title('Sorted pathway presence')


%% Pathways per cluster analysis
csm_unisub = unique(string(Context_Specific_model.subSystems));
csm_unisub = csm_unisub(2:end);
csm_c1_ind = find(contains(Context_Specific_model.rxns,'_1'));

subsys = Context_Specific_model.subSystems(csm_c1_ind);
subsys = string(subsys);
subsys_uni = unique(subsys); % "" empty cell, rxns that do not belong to a pathway (yet to determine)
[subsys_count, subsys_names] = histcounts(categorical(subsys));
sub_diff = cellstr(setdiff(csm_unisub,subsys_names));
subsys_names = [subsys_names,sub_diff'];
subsys_count = [subsys_count, zeros(1,numel(sub_diff))];
table_subsys = table(subsys_names', subsys_count');
table_subsys.Properties.VariableNames = {'subsys_names','subsys_count'};
table_subsys1 = sortrows(table_subsys,{'subsys_count'},{'descend'});
%writetable(table_subsys,'.\Context_Subsys_table.txt')


% csm_unisub = unique(string(C1.subSystems));
% csm_unisub = csm_unisub(2:end);
% csm_c1_ind = find(contains(C1.rxns,'_5'));
% 
% subsys = C1.subSystems(csm_c1_ind);
% subsys = string(subsys);
% subsys_uni = unique(subsys); % "" empty cell, rxns that do not belong to a pathway (yet to determine)
% [subsys_count, subsys_names] = histcounts(categorical(subsys));
% sub_diff = cellstr(setdiff(csm_unisub,subsys_names));
% subsys_names = [subsys_names,sub_diff'];
% subsys_count = [subsys_count, zeros(1,numel(sub_diff))];
% table_subsys = table(subsys_names', subsys_count');
% table_subsys.Properties.VariableNames = {'subsys_names','subsys_count'};
% table_subsysC5 = sortrows(table_subsys,{'subsys_count'},{'descend'});


% determine the reaction count in recon3d, divide csm with recon3d
r3d_paths = string(model.subSystems);
r3d_uni_paths = unique(r3d_paths);
[r3d_pathcount, r3d_pathnames] = histcounts(categorical(r3d_paths));
table_r3d_paths = table(r3d_pathnames',r3d_pathcount');
table_r3d_paths.Properties.VariableNames = {'subsys_names','paths_count'};

table_int1 = join(table_subsys1,table_r3d_paths);
table_int2 = join(table_subsys2,table_r3d_paths);
table_int3 = join(table_subsys3,table_r3d_paths);
table_int4 = join(table_subsys4,table_r3d_paths);
table_int5 = join(table_subsys5,table_r3d_paths);

table_frac1 = [table_int1(:,1),array2table(table_int1.subsys_count./table_int1.paths_count)];
table_frac1.Properties.VariableNames = {'subsys_names','Epithelial_1'};
table_frac2 = [table_int2(:,1),array2table(table_int2.subsys_count./table_int2.paths_count)];
table_frac2.Properties.VariableNames = {'subsys_names','Epithelial_2'};
table_frac3 = [table_int3(:,1),array2table(table_int3.subsys_count./table_int3.paths_count)];
table_frac3.Properties.VariableNames = {'subsys_names','Epithelial_3'};
table_frac4 = [table_int4(:,1),array2table(table_int4.subsys_count./table_int4.paths_count)];
table_frac4.Properties.VariableNames = {'subsys_names','T-cell'};
table_frac5 = [table_int5(:,1),array2table(table_int5.subsys_count./table_int5.paths_count)];
table_frac5.Properties.VariableNames = {'subsys_names','B-cell'};


% Bar graph of subsystem count
[~,I] = sort(table_subsys.subsys_count,'descend'); %sort by descreasing subsystems counts
figure
bar(table_subsys.subsys_count(I));
xlim([1 numel(table_subsys.subsys_names)]);
% xticklabels(T.pathway_names(I)) %enable this to see the names, the figure will be very small
xticks(1:numel(table_subsys.subsys_names));
xticklabels(table_subsys.subsys_names);
xtickangle(60);
title('Number of reactions per subsystem');

% Heatmap
table_frac1 = sortrows(table_frac1);
% table_subsys1.Properties.VariableNames = {'subsys_names','T-cell'};
heat_table = join(table_frac1, table_frac2);
heat_table = join(heat_table, table_frac3);
heat_table = join(heat_table, table_frac4);
heat_table = join(heat_table, table_frac5);
heat_table = sortrows(heat_table,'Epithelial_1','descend');
heat_array = table2array(heat_table(:,2:end));

heat_table_cut = heat_table(1:15,:);
heat_array_cut = table2array(heat_table_cut(:,2:end));

figure('WindowState','maximized')
% lngth = 20;
% Blue = [0, 0, 1];
% white = [1, 1, 1];
% colors_p = [linspace(white(1),Blue(1),lngth)', linspace(white(2),Blue(2),lngth)', linspace(white(3),Blue(3),lngth)'];
I = imagesc(heat_array);
caxis([0 max(max(heat_array))])
title('Model''s Cluster integrated Pathways based on Recon3D')
ylabel('Pathways','FontSize',16)
xlabel('Clusters','FontSize',16)
xticks([1:5]); 
xticklabels(strrep(heat_table.Properties.VariableNames(2:end),'_',' '))
xtickangle(0)
yticks([1:98])
yticklabels(table2cell(heat_table(:,1)))
ytickangle(30)
colorbar;
ax = gca;
c = ax.FontSize;
ax.FontSize = 11;

figure('WindowState','maximized')
I = imagesc(heat_array_cut);
caxis([0 max(max(heat_array_cut))])
title('Model''s Cluster integrated Pathways based on Recon3D')
ylabel('Pathways','FontSize',16)
xlabel('Clusters','FontSize',16)
xticks([1:5]); 
xticklabels(strrep(heat_table_cut.Properties.VariableNames(2:end),'_',' '))
xtickangle(0)
yticks([1:98])
yticklabels(table2cell(heat_table_cut(:,1)))
ytickangle(30)
colorbar;
ax = gca;
c = ax.FontSize;
ax.FontSize = 16;

cgo_J = clustergram(heat_array',...
    'ColumnLabels', table2cell(heat_table(:,1)),...
    'RowLabels', Headers(2:end),...
    'ColumnLabelsRotate',340, ...
    'Cluster', 'all', ...
    'Annotate', 'false',...
    'symmetric','False',...
    'AnnotColor','k')
addTitle(cgo_J,{'Model''s Cluster integrated Pathways based on Recon3D'});
plot(cgo_J);
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
fig_gcf = findall(0,'type','figure', 'tag', 'Clustergram');
set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)


cgo_J1 = clustergram(heat_array_cut',...
    'ColumnLabels', table2cell(heat_table_cut(:,1)),...
    'RowLabels', Headers(2:end),...
    'ColumnLabelsRotate',340, ...
    'Cluster', 'all', ...
    'Annotate', 'false',...
    'symmetric','False',...
    'AnnotColor','k')
addTitle(cgo_J1,{'Model''s Cluster integrated Pathways based on Recon3D'});
plot(cgo_J1);
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
fig_gcf = findall(0,'type','figure', 'tag', 'Clustergram');
set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)

%% Single gene deletion analysis
Context_Specific_model2 = Context_Specific_model;
Context_Specific_model2.genes = extractBefore(Context_Specific_model2.genes,'_');

Context_Specific_model2.c(find(ismember(Context_Specific_model2.rxns, 'biomassreaction_1')))=1;
Context_Specific_model2.c(find(ismember(Context_Specific_model2.rxns, 'biomassreaction_2')))=0;
Context_Specific_model2.c(find(ismember(Context_Specific_model2.rxns, 'biomassreaction_3')))=0;
Context_Specific_model2.c(find(ismember(Context_Specific_model2.rxns, 'biomassreaction_4')))=0;
Context_Specific_model2.c(find(ismember(Context_Specific_model2.rxns, 'biomassreaction_5')))=0;
cind = find(Context_Specific_model2.c);
Context_Specific_model2.c(cind) = clust_ratio;

unique(Context_Specific_model2.c)
Context_Specific_model2.rxns(find(Context_Specific_model2.c))
checkObjective(Context_Specific_model2);

tic
[grRatio,grRateKO,grRateWT5,hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(Context_Specific_model2,'FBA',[],0,1);
toc
table_sgd5 = [grRatio,grRateKO];

unique(table_sgd5(:,1))
ess_ind = find(table_sgd5(:,1)<0.55);
ess_genes5 = geneList(ess_ind);


%% Single reaction deletion
csm = Context_Specific_model2;
rxn_count = count(Context_Specific_model2.rxns(1:18825),'_');
rxn_list = Context_Specific_model.rxns(1:18825);
csm.rxns(1:18825) = extractBefore(rxn_list,'_');
csm.c(find(ismember(Context_Specific_model2.rxns, 'biomassreaction')))=1;
cind = find(csm.c);
csm.c(cind) = clust_ratio;

find(ismember(csm.rxns,'biomassreaction'))
csm.c(3448)=0;
csm.c(7218)=0;
csm.c(10730)=0;
csm.c(13642)=0;
csm.c(16185)=0;

unique(csm.c)
csm.rxns(find(csm.c))
checkObjective(csm);

tic
[grRatior, grRateKOr, grRateWTr, hasEffectr, delRxnr, fluxSolutionr] = singleRxnDeletion_MP(csm, 'FBA',csm.rxns);
toc
table_srd = [grRatior,grRateKOr];
unique(table_srd(:,1))
essr_ind = find(table_srd(:,1)<0.55);
ess_rxns = csm.rxns(essr_ind);


csm.grRules(find(contains(csm.rxns,ess_rxns)))


%% Manual KO
csm = Context_Specific_model2;
csm.lb(find(ismember(csm.rxns,'TRDR')))=0
csm.ub(find(ismember(csm.rxns,'TRDR')))=0
optimizeCbModel(csm)
optimizeCbModel(Context_Specific_model2)
7296
[r,c]=find(csm.rxnGeneMat(:,find(ismember(csm.genes,'7296'))))
csm.rxns(r)

[agrRatio,agrRateKO,agrRateWT,~, ~, ~, ageneList] = singleGeneDeletion_MISB(csm,'FBA',{'7296'});
[agrRatio, agrRateKO, agrRateWT, ahasEffect, adelRxns, afluxSolution, ageneList] = singleGeneDeletion_MISB(csm,'FBA', model.genes, 0, 1)
unique(csm.c)

%manual gene ko
tempcsm = Context_Specific_model2;
gr_ind = find(contains(tempcsm.grRules,'6566'));
tempcsm.lb(gr_ind) = 0;
tempcsm.ub(gr_ind) = 0;
tempcsm.c(find(ismember(tempcsm.rxns, 'biomassreaction_4')))=1;
checkObjective(tempcsm)
optimizeCbModel(tempcsm)
