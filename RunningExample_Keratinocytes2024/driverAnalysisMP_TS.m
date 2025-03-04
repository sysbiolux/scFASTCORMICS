% driver analysis, TS 07/23
%% scFASTCORMICS results
clearvars -except solverOK, clc, close all force

%% Load model/define the clusters. conditions and the medium
load('Results09-Jul-2024.mat')
Results_keep=Results
Results_keep.multi_cell_population=Results_keep.multi_cell_population_model;

load('simpleRev_recon3.mat')
Recon=model;

[NUM,TXT,RAW]=xlsread('medium_withconc.xlsx');
medium=[RAW(:,2); 'hdl_hs[e]'; 'nh4[e]'; 'so3[e]'; 'tsul[e]'; 'tag_hs(e)']
T=table(RAW(:,1),RAW(:,2),NUM,'VariableNames',{'medium_rxns_keep','medium_composition','met_Conc_mM'})

names_col={'timepoint1'}

ClusterNames={'H','M','P','TD1','TD2','F1','F2'}'

CellNumberPerCluster=table([897 1941 1279 1470 686 156 58]','RowNames',ClusterNames,'VariableNames',names_col)

biomass_reaction='biomassreaction' %sc model: biomassreaction_1 etc

changeCobraSolver('ibm_cplex')

%% model sizes
Modelsizes=zeros(size(Results_keep,1),1);
for counter=1:numel(Results_keep)
    Modelsizes(counter)=numel(Results_keep(counter).multi_cell_population_model.rxns);
end
disp('Model sizes:')
disp(names_col')
disp(Modelsizes)

%% model sizes per cluster
clusterSize=zeros(size(Results_keep,2),numel(ClusterNames));
for counter=1:size(clusterSize,1)
    %     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population_model;
    for counter2=1:size(clusterSize,2)
        idx = strfind(model.rxns, ['_' num2str(counter2)]);
        idx = find(not(cellfun('isempty', idx)));
        clusterSize(counter,counter2)=numel(idx);
    end
end
disp('Model sizes per clusters (cols):')
clusterSize=array2table(clusterSize);
clusterSize.Properties.VariableNames=ClusterNames;
clusterSize.Row=names_col;
disp(clusterSize)

%% Jaccard plot Whole Models
if numel(Results_keep)<2
    warning('Not enough models for model comparison')
else
    J=nan(numel(Results_keep));
    for counter=1:numel(Results_keep)
        for counter2=1:numel(Results_keep)
            A1=find(Results_keep(counter).A);
            A2=find(Results_keep(counter2).A);
            J(counter,counter2)=numel(intersect(A1,A2))/numel(union(A1,A2));
        end
    end
    disp('Jaccard similarity Whole Models:')
    
    altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
    cgo_J = clustergram(J,...
        'RowLabels', names_col,...
        'ColumnLabels', names_col,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all', ...
        'Annotate', 'true',...
        'symmetric','False',...
        'AnnotColor','k',...
        'Colormap', altcolor);
    addTitle(cgo_J,{'Jaccard similarity Whole Models'});
    plot(cgo_J);
end

%% Jaccard plot per Cluster over all models

%%% Why are reactions in each cluster in generic model slightly different?
if numel(Results_keep)<2
    warning('Not enough models for cluster comparison over all models')
else
    for counter=1:numel(clusterSize.Properties.VariableNames) %per cluster
        J=nan(numel(clusterSize.Properties.VariableNames));
        % identify from generic model rxns in specific cluster
        model=Results_keep(1).ExpandedInputModel;
        idx = strfind(model.rxns, ['_' num2str(counter)]);
        idx = find(not(cellfun('isempty', idx)));
        
        for counter2=1:numel(Results_keep) %per model
            for counter3=1:numel(Results_keep) %per model
                A1=find(Results_keep(counter2).A);
                A1=intersect(A1,idx);
                A2=find(Results_keep(counter3).A);
                A2=intersect(A2,idx);
                J(counter2,counter3)=numel(intersect(A1,A2))/numel(union(A1,A2));
            end
        end
        disp(['Jaccard similarity Whole Cluster:' num2str(counter)])
        
        altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
            255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
        cgo_J = clustergram(J,...
            'RowLabels', ClusterNames,...
            'ColumnLabels', ClusterNames,...
            'ColumnLabelsRotate',340, ...
            'Cluster', 'all', ...
            'Annotate', 'true',...
            'symmetric','False',...
            'AnnotColor','k',...
            'Colormap', altcolor);
        addTitle(cgo_J,{['Jaccard similarity Cluster: ' num2str(counter)]});
        plot(cgo_J);
    end
end

%% Jaccard plot per all Clusters within one model
for counter=1:numel(Results_keep) %per model
    J=nan(numel(clusterSize.Properties.VariableNames));
    %     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population_model;
    
    for counter2=1:numel(clusterSize.Properties.VariableNames) %per cluster
        idx1 = strfind(model.rxns, ['_' num2str(counter2)]);
        idx1 = find(not(cellfun('isempty', idx1)));
        A1 = cellfun(@(S) S(1:end-2), model.rxns(idx1), 'Uniform', 0);
        for counter3=1:numel(clusterSize.Properties.VariableNames) %per cluster
            idx2 = strfind(model.rxns, ['_' num2str(counter3)]);
            idx2 = find(not(cellfun('isempty', idx2)));
            A2 = cellfun(@(S) S(1:end-2), model.rxns(idx2), 'Uniform', 0);
            J(counter2,counter3)=numel(intersect(A1,A2))/numel(union(A1,A2));
        end
    end
    disp(['Jaccard similarity Model:' num2str(counter)])
    
    altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
    cgo_J = clustergram(J,...
        'RowLabels', clusterSize.Properties.VariableNames,...
        'ColumnLabels', clusterSize.Properties.VariableNames,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all', ...
        'Annotate', 'true',...
        'symmetric','False',...
        'AnnotColor','k',...
        'Colormap', altcolor);
    addTitle(cgo_J,{['Jaccard similarity Model: ' num2str(counter)]});
    plot(cgo_J);
end

%% Pathway anaysis model centred
multiCellModel=Results_keep(1).multi_cell_population_model;
ExpandedInputModel=Results_keep(1).ExpandedInputModel;
clear counter counter2 counter3 idx1 idx2 J A1 A2

Pathways_keep=struct();
for i=1:numel(names_col)
    multiCellModel=Results_keep(i).multi_cell_population_model;
    [Pathways]=PathwayAnalysisMultiCellModel(multiCellModel,ExpandedInputModel,  clusterSize.Properties.VariableNames);
    
    Pathways_keep(i).Analysis=Pathways;
    Pathways_keep(i).multiCellModel=multiCellModel;
    %Making a plot that show interesting pathways differences
    %between the different clusters (Tumor)
    
    %Finding the size of the top 30 most differentially present pathways
    Diff=zeros(size(Pathways,1),1);
    %Finding top thirty pathways differential
    for ii = 1:height(Pathways)
        maxValue = max(Pathways{ii,:});
        minValue = min(Pathways{ii,:});
        Diff(ii,1) = (maxValue - minValue);
    end
    Pathways2=[Pathways, table(Diff)];
    index=find(~isnan(Pathways2.Diff));
    
    ClusterRatioShort=Pathways2(index(1:30),:);
    
    cgo_Path = clustergram(table2array(ClusterRatioShort(:,1:end-1)),...
        'RowLabels',ClusterRatioShort.Row ,...
        'ColumnLabels',ClusterRatioShort.Properties.VariableNames(1:end-1),...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all',...
        'LabelsWithMarkers',true,...
        'Colormap',redbluecmap,...
        'DisplayRange',1000,...
        'DisplayRatio', [0.2 0.2]);
    
    addTitle(cgo_Path, cellstr(strcat('Pathway Analysis: ',' ', names_col{i})));
    
    graph = plot(cgo_Path);
    set(graph, 'Clim', [0,1])
    colormap(graph,"copper")
    
    figureHandle = gcf;
    %# make all text in the figure to wanted size and bold
    set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)
    
end


%% Pathway anaysis cluster centred
[ClusterRatioKeep]=PathwayAnalysisMultiCellModel4Cluster(Results_keep,ExpandedInputModel, clusterSize.Properties.VariableNames, names_col);
for ii=1:numel(ClusterRatioKeep)
    % Finding the size of the top 30 most differentially present pathways
    ClusterRatio=ClusterRatioKeep(ii).cluster;
    Diff=zeros(size(ClusterRatio,1),1);
    
    for i=1:size(ClusterRatio,1)
        %Finding top thirty pathways differential
        maxValue = max(ClusterRatio{i,:});
        minValue = min(ClusterRatio{i,:});
        Diff(i,1) = (maxValue - minValue);
    end
    ClusterRatio2=[ClusterRatio, table(Diff)];
    index=find(~isnan(ClusterRatio2.Diff));
    
    ClusterRatioShort=ClusterRatio2(index(1:30),:);
    
    if size(ClusterRatioShort(:,1:end-1),2)<2
        warning('Not enough models for cluster comparison over all models')
    else
        cgo_Path = clustergram(table2array(ClusterRatioShort(:,1:end-1)),...
            'RowLabels',ClusterRatioShort.Row ,...
            'ColumnLabels',ClusterRatioShort.Properties.VariableNames(1:end-1),...
            'ColumnLabelsRotate',340, ...
            'Cluster', 'all',...
            'LabelsWithMarkers',true,...
            'Colormap',redbluecmap,...
            'DisplayRange',1000,...
            'DisplayRatio', [0.2 0.2]);
        
        addTitle(cgo_Path, cellstr(strcat('Pathway Analysis: ',' ', clusterSize.Properties.VariableNames{ii})));
        
        graph = plot(cgo_Path);
        set(graph, 'Clim', [0,1])
        colormap(graph,"copper")
        
        figureHandle = gcf;
        %# make all text in the figure to wanted size and bold
        set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
        set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)
    end
end

%% run this for all models ... (change index in script!)
% FBA
% sortFlag=1; %1:max-min; 2:max; 3:min

[AnalysisKeep,ExchangesKeep]=ExchangesPrediction_weighted(Results_keep,Recon,T,CellNumberPerCluster,biomass_reaction);
%%
plotCutoff=0.1
plotMax=0.5

for i=1:numel(ExchangesKeep)
    ClusterExchanges=ExchangesKeep(i).Exchanges;
    sum_of_matrix_rows = sum(abs(table2array(ClusterExchanges)),2)>plotCutoff;
    ClusterExchangesShort=ClusterExchanges(sum_of_matrix_rows,:);
    
    ClusterExchangesShort=sortrows(ClusterExchangesShort,size(ClusterExchangesShort,2),'descend');
    
    temp=table2array(ClusterExchangesShort);
    temp(temp>plotMax)=plotMax;
    temp(temp<-plotMax)=-plotMax;
    
    cgo_Path = clustergram(temp,...
        'RowLabels',ClusterExchangesShort.Row ,...
        'ColumnLabels', ClusterExchangesShort.Properties.VariableNames,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all',...
        'LabelsWithMarkers',true,...
        'Colormap',redbluecmap,...
        'DisplayRange',10,...
        'DisplayRatio', [0.2 0.2]);
    
    %     addTitle(cgo_Path, cellstr(strcat('FBA of exchanged metabolites','Model:', names_col(i), ', sortflag: ',num2str(sortFlag))))
    addTitle(cgo_Path, cellstr(strcat('FBA of exchanged metabolites','Model:', names_col(i))))
    
    plot(cgo_Path);
    
    figureHandle = gcf;
    %# make all text in the figure to size 14 and bold
    set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)
    
    % red (positive values) indicates the metabolite leaves the corresponding cluster
    % blue (negative values) indicates the metabolite enters the corresponding cluster
    fileName=['data' num2str(i) '.xlsx'];
    if exist(fileName, 'file') == 2
        
        delete(fileName)
    end
    % % %     writetable(ExchangesKeep(i).Exchanges,fileName,'WriteRowNames',1)
    writetable(ClusterExchangesShort,fileName,'WriteRowNames',1)
end

%save(['FBAmodel' num2str(modelNr)],'FBAweighted','multiCellModel')
%%
if numel(Results_keep)<2
    warning('Not enough models for FBA comparison')
else
    [ComparisonFlux2Models,res4t, FluxPerCluster,FluxPerPathway, FluxPerPathwayPerCluster]=compareFBA(AnalysisKeep, Results_keep);
    k=0;
    for i=1:numel(FluxPerCluster)
        disp('reactions carrying a flux per cluster')
        disp(FluxPerCluster(i).FluxPerCluster)
    end
    for i=1:numel(FluxPerPathway)
        
        disp('reactions carrying a flux per pathways')
        disp(FluxPerPathway)
    end
    for counterA=1:numel(Results_keep)
        for counterB=1:numel(Results_keep)
            if counterA<counterB
                
                k=k+1;
                
                disp('relative flux difference per subSystem and cluster')
                res4= res4t(k).rest;
                res4m=table2array(res4);
                temp=sort(max(res4m,[],2),'descend');
                show=find(max(res4m,[],2)>=temp(10)); %show top10
                resTop10=res4(show,:);
                
                
                fileName=['compareFBA_' num2str(counterA) '_vs' num2str(counterB) '.xlsx'];
                if exist(fileName, 'file') == 2
                    
                    delete(fileName)
                end
                writetable(res4,fileName,'WriteRowNames',1)
            end
        end
    end
end

%% extract indivudal GEMs from community model
clearvars -except solverOK, clc, close all force

load('Results09-Jul-2024.mat')
Results_keep=Results
Results_keep.multi_cell_population=Results_keep.multi_cell_population_model;

model=Results_keep(1).multi_cell_population_model

modelToExtract=1 %which cluster to extract

% Extract individual GEM
% This will only include the reactions of the desired cluster, 
% but the resulting model will not be consistent.
% For consistency reactions of the other clusters will be necessary.
% See below. 
idx = strfind(model.rxns, ['_' num2str(modelToExtract)]);
idxToKeep = find(not(cellfun('isempty', idx)));
numel(idxToKeep)
modelInd   = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),idxToKeep)));
modelInd.description='Reactions Cluster 1'

% Extract consistent individual GEM
% This will include some reactions of other submodels to ensure consistency
A      = fastcore_4_rfastcormics(idxToKeep, model, 1e-4, []);
modelIndCons   = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A)));
modelIndCons.description='Model Cluster 1'

A2 = fastcc_4_rfastcormics(modelIndCons, 1e-4, 1);
numel(A2)
