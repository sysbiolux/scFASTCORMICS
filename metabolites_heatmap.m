
%load('P:\Jimmy Bis\workspace_V7.mat')
table_flux = [cmodel_rxnstrans printRxnFormula(Context_Specific_model, Context_Specific_model.rxns(rxns_need2)) num2cell(full(minFlux)) num2cell(full(maxFlux)) Context_Specific_model.subSystems(rxns_need2)];
%% keep only higher fluxes
clear range
FLUXON=zeros(size(table_flux,1),1)
for i=1:size(table_flux,1)
    FLUXON(i,1)=table_flux{i,4}> 1e-4 | table_flux{i,3}< -1e-4 ;
end
%find significant flux
Significant_flux=table_flux(find(FLUXON),:);

% both direction?
direction=zeros(size(Significant_flux,1),1);
for i=1:size(Significant_flux,1)
    
    if (Significant_flux{i,3})< -1 &&  (Significant_flux{i,4})> 1
        direction(i,1)=2;
    elseif (Significant_flux{i,3}) <-1
        direction(i,1)=-1;
    elseif (Significant_flux{i,4}) >1
        direction(i,1)=1;
    else
        direction(i,1)=0;
    end
end


% metabolites implicated in the  interaction

[m,~]=find(Context_Specific_model.S(:, ismember(Context_Specific_model.rxns,Significant_flux(:,1))));
mets=Context_Specific_model.mets(m);
mets=mets(contains(mets,'[u]'));

mets=unique(mets);

mets_matrix=zeros(numel(mets),5);
for i=1:numel(mets)
    
    [~,r]=find(Context_Specific_model.S( ismember(Context_Specific_model.mets, mets(i)),:));
    [names_reaction, IA, IB]=intersect(Context_Specific_model.rxns(r),Significant_flux(:,1));
    Context_Specific_model.S(ismember(Context_Specific_model.mets, mets(i)),r(IA));
    for j=1:numel(IA)
        if Context_Specific_model.S(ismember(Context_Specific_model.mets, mets(i)),r(IA(j)))>=1 && direction(IB(j)) ==1
           
            clu=extractAfter(Context_Specific_model.rxns(r(IA(j))),'_');
            clu=str2double(clu{1});
            if mets_matrix(i,clu)==0
            mets_matrix(i,clu)=-1;
            else
                if mets_matrix(i,clu)==2 ||   mets_matrix(i,clu)==1
                    
                mets_matrix(i,clu)=2
                elseif mets_matrix(i,clu)==- 1 ||mets_matrix(i,clu)==0
                    
                end
            end
            
        elseif Context_Specific_model.S(ismember(Context_Specific_model.mets, mets(i)),r(IA(j)))<=-1 && direction(IB(j)) ==1
           
            clu=extractAfter(Context_Specific_model.rxns(r(IA(j))),'_');
            clu=str2double(clu{1});
            if mets_matrix(i,clu)==0 | mets_matrix(i,clu)==1
            mets_matrix(i,clu)=1;
            elseif mets_matrix(i,clu)==-1 | mets_matrix(i,clu)==2
                 mets_matrix(i,clu)=2
            else
                bbb
                
            end
        elseif Context_Specific_model.S(ismember(Context_Specific_model.mets, mets(i)),r(IA(j)))<=-1 && direction(IB(j)) ==-1
           
            clu=extractAfter(Context_Specific_model.rxns(r(IA(j))),'_');
            clu=str2double(clu{1});
            if mets_matrix(i,clu)==0 || mets_matrix(i,clu)==-1
                
            mets_matrix(i,clu)=-1;
            elseif  mets_matrix(i,clu)==2
                
            elseif mets_matrix(i,clu)==1
                mets_matrix(i,clu)=2
            else
                sss
            end
        elseif direction(IB(j))==0
           
            clu=extractAfter(Context_Specific_model.rxns(r(IA(j))),'_');
            clu=str2double(clu{1});
            if mets_matrix(i,clu)==0
            mets_matrix(i,clu)=0;
            else
                
                % if there was sth previously then do not change it
            end
        elseif direction(IB(j))==2
            clu=extractAfter(Context_Specific_model.rxns(r(IA(j))),'_');
            clu=str2double(clu{1});
            
            mets_matrix(i,clu)=2;
            
        elseif Context_Specific_model.S(ismember(Context_Specific_model.mets, mets(i)),r(IA(j)))>=1 && direction(IB(j))==-1
            
            clu=extractAfter(Context_Specific_model.rxns(r(IA(j))),'_');
            clu=str2double(clu{1});
            if mets_matrix(i,clu)==0 | mets_matrix(i,clu)==1
            mets_matrix(i,clu)=1;
            elseif mets_matrix(i,clu)==2 | mets_matrix(i,clu)==1
                mets_matrix(i,clu)==2
            else
                
            end
        else
            ssss
                  
            
        end
    end
end

T=table(mets, mets_matrix)

match=sum(mets_matrix~=0,2)>1
T2=T(match,:)

T2_values = table2array(T2(:,2));

headers = {'Epithelial 1','Epithelial 2','Epithelial 3','T-cell','B-cell'}';
mets_uni_table = extractBefore(string(table2cell(T2(:,1))),'[');
mets_uni_table = append(mets_uni_table,'[e]');
for i=1:size(T2,1)
mets_uni_table(i,2) = model.metNames(find(ismember(model.mets,mets_uni_table(i,1))));
end
T2_fullnames = [table(mets_uni_table(:,2)),T2(:,2)];
metsNames = mets_uni_table(:,2);
metsNames(20) = 'Met 1';
metsNames(35) = 'Met 2';
metsNames(40) = 'Met 3';

figure('WindowState','maximized')
I = clustergram(T2(:,2:end));
%caxis([-1 1])
title('Metabolite inter-cluster flux')
ylabel('Metabolites','FontSize',16)
xlabel('Clusters','FontSize',16)
xticks([1:5]); 
xticklabels({'Epithelial 1','Epithelial 2','Epithelial 3','T-cell','B-cell'})
xtickangle(0)
yticks([1:50])
yticklabels(T2(:,1))
ytickangle(30)
colorbar;
ax = gca;
c = ax.FontSize;
ax.FontSize = 11;

% defaultPosition = get(0,'DefaultFigurePosition');
% screensize = get(groot, 'Screensize');
% screensize = [1,1,800,800];
% set(0, 'DefaultFigurePosition', screensize);
              
cgo_J = clustergram(T2_values,...
    'RowLabels', mets_uni_table(:,2),...
    'RowLabelsRotate', 360,...
    'ColumnLabels', headers,...
    'ColumnLabelsRotate', 360,...
    'Cluster', 'Column', ...
    'symmetric','False',...
    'Colormap', redbluecmap)
% addTitle(cgo_J,{'Metabolite inter-cluster flux'});
set(0,'ShowHiddenHandles','on')
allhnds = get(0,'Children');
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 8)
plot(cgo_J)

cgo_J = clustergram(T2_values',...
    'RowLabels', headers,...
    'RowLabelsRotate', 360,...
    'ColumnLabels', metsNames,...
    'ColumnLabelsRotate', 90,...
    'Cluster', 'Column', ...
    'symmetric','False',...
    'Colormap', redbluecmap)
% addTitle(cgo_J,{'Metabolite inter-cluster flux'});
set(0,'ShowHiddenHandles','on')
allhnds = get(0,'Children');
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 8)
plot(cgo_J)