function[AnalysisKeep,ExchangesKeep]=ExchangesPrediction_weighted(Results_keep,model,T,CellNumberPerCluster,biomass_reaction)
AnalysisKeep=struct();
ExchangesKeep=struct();

for clusterrouned=1:numel(Results_keep)
modelNr=clusterrouned;

multiCellModel=Results_keep(modelNr).multi_cell_population;
LastTwo = cellfun(@(S) S(end-1:end), multiCellModel.rxns, 'Uniform', 0);
clusters=cell(size(LastTwo));
match=contains(LastTwo,'_');
clusters(match) = cellfun(@(S) S(end), LastTwo(match), 'Uniform', 0);

LastThree = cellfun(@(S) S(end-2:end), multiCellModel.rxns, 'Uniform', 0);
match=~contains(LastTwo,'_') & ~contains(LastTwo,'[e]');
clusters(match) = cellfun(@(S) S(end-1), LastThree(match), 'Uniform', 0);

% set medium uptake rates to concentrations
[temp,IA,IB]=intersect(multiCellModel.rxns,T.medium_rxns_keep);
tempS=nonzeros(multiCellModel.S(:,IA));
for counter=1:numel(IA)
    if tempS(counter)>0
        multiCellModel.ub(IA(counter))=T.met_Conc_mM(IB(counter));
    else
        multiCellModel.lb(IA(counter))=-T.met_Conc_mM(IB(counter));
    end
end

%setting a 10% lower bounds for each cluster's biomass reaction (10% of the
%max biomass number obtained when we fully optimise for the specific
%biomass reaction)
% % % biomass=multiCellModel.rxns(contains(multiCellModel.rxns,'biomass'));
%for all metabolites that contains [u], R are the reactions that involve
%that metabolite -->

% % % for i=1: size(CellNumberPerCluster,2)
% % %     ordre = zeros(1,size(CellNumberPerCluster,2));
% % %     if i==i
% % %         ordre(i) = 1;
% % %     end
% % %     multiCellModel.c(contains(multiCellModel.rxns,'biomassmaintenance')) = ordre;
% % %     %Flux Balance Analysis (FBA) is used to find the max optimised value if the
% % %     %model was optimising only for this specific biomass reaction
% % %     FBAsolMax1cluster=optimizeCbModel(multiCellModel);
% % %     tenpercent = 0.1*FBAsolMax1cluster.f;
% % %     %set lower bound
% % %     multiCellModel.lb(contains(multiCellModel.rxns,biomass(i))) = tenpercent;
% % % end
% % % clear i ordre
% % % TotalCells=sum(table2array(CellNumberPerCluster(clusterrouned,:)));
% % % BiomassFix=table2array(CellNumberPerCluster(clusterrouned,:))./repmat(TotalCells,size(CellNumberPerCluster(clusterrouned,:),1),1);
% % % multiCellModel.c(contains(multiCellModel.rxns,'biomassmaintenance')) = BiomassFix;
% % % FBAweighted = optimizeCbModel(multiCellModel,'max','zero');
% % % indx_non_zero_flux=FBAweighted.v~=0;
% % % clear BiomassFix biomass FBAsolMax1cluster tenpercent TotalCells

%% NEW idea: add new overall biomass reaction, with biomass fix as stoichiometric coefficients
temp=table2array(CellNumberPerCluster(:,modelNr));
TotalCells=sum(temp);
BiomassFix=temp./repmat(TotalCells,size(temp,1),1);
sum(BiomassFix);

mm=multiCellModel;
overallBiomassRxn=[];
clear biomass
for i=1:size(CellNumberPerCluster,1)
    biomass(i) = cellstr([biomass_reaction '_' num2str(i)]);
end 
for counter=1:numel(biomass)
    temp = printRxnFormula(mm,'rxnAbbrList',cellstr(biomass(counter)));
    newRxn=[temp{:} ['+ BIOMASS' num2str(counter)]];
    mm = removeRxns(mm,biomass{counter});
    mm = addReaction(mm,biomass{counter},'reactionFormula',newRxn);
    printRxnFormula(mm,'rxnAbbrList',biomass{counter});
    if counter < numel(biomass)
        overallBiomassRxn=[overallBiomassRxn, [num2str(BiomassFix(counter,modelNr)) ' BIOMASS' num2str(counter) ' + ']];
    else
        overallBiomassRxn=[overallBiomassRxn, [num2str(BiomassFix(counter,modelNr)) ' BIOMASS' num2str(counter) ' -> ']];
    end
end
mm = addReaction(mm,'biomassreaction_All','reactionFormula',overallBiomassRxn);
mm.rxns(contains(mm.rxns,'biomass'))

mm.c=0*mm.c;
mm.c(find(contains(mm.rxns,'biomassreaction_All')))=1; %_All

multiCellModel=mm;

FBAweighted = optimizeCbModel(multiCellModel,'max','zero')
FBAweighted.x(end-11:end)
indx_non_zero_flux=FBAweighted.v~=0;

%%
ReactionsNonZerFluxFBA = multiCellModel.rxns(indx_non_zero_flux);
Flux_non_zero_FBA = FBAweighted.v(indx_non_zero_flux);
TableFlux=table(ReactionsNonZerFluxFBA, Flux_non_zero_FBA);
clear Flux_non_zero_FBA indx_non_zero_flux
%--> finding the U reactions and searching for the index so that we  can find
%the metabolites present in it with the S matrix
UMetsID=find(contains(multiCellModel.mets,'[u]'));

%%
ClusterExchanges = zeros(length(UMetsID),size(CellNumberPerCluster,1)+1);
% ClusterExchanges.Properties.RowNames =  list_of_U_metabolites_that_have_flux(:,1);
% ClusterExchanges.Properties.VariableNames = list_U_reactions_FBAflux(:,1);

for i = 1:length(UMetsID)
    
    [~,Rindx]=find(multiCellModel.S(UMetsID(i),:));
    
    tmp=zeros(numel(Rindx),size(CellNumberPerCluster,1)+1);
    %to note down the flux of the reactions for that specific metabolite
    %decide if the metabolite is consumed or produced depending on which
    %side of the reaction it is (S matrix value) and depending on its flux
    % --> THIS  CODE HAS BEEN CORRECTED THROUGH THE TROUBLE SHOOTING PART
    % (SEE DOWN BELOW)
    for v=1:length(Rindx)
        %right and positive flux --> produced
        [~,~,IB]=intersect(multiCellModel.rxns(Rindx(v)),TableFlux.ReactionsNonZerFluxFBA);
        if ~isempty(IB)
            
            c=(clusters(Rindx(v)));
            if ~strcmp(c,'e')
                c=str2double(c{1});
            else
                c=size(CellNumberPerCluster,1)+1;
            end
            if multiCellModel.S(UMetsID(i),Rindx(v)) > 0 && TableFlux.Flux_non_zero_FBA(IB) > 0
                tmp(v,c)= table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
                %left and positive flux --> consumed
            elseif multiCellModel.S(UMetsID(i),Rindx(v)) < 0 && TableFlux.Flux_non_zero_FBA(IB) > 0
                tmp(v,c)= -table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
                %right and negative flux --> consumed
            elseif multiCellModel.S(UMetsID(i),Rindx(v)) > 0 && TableFlux.Flux_non_zero_FBA(IB) < 0
                tmp(v,c)= table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
                %left and negative flux --> produced
            elseif multiCellModel.S(UMetsID(i),Rindx(v)) < 0 && TableFlux.Flux_non_zero_FBA(IB) < 0
                tmp(v,c)= -table2array(TableFlux(ismember(TableFlux.ReactionsNonZerFluxFBA,multiCellModel.rxns(Rindx(v))),2));
            end
            
            
            ClusterExchanges(i,:)=sum(tmp,1);
        end
    end
end
ClusterExchanges=array2table(ClusterExchanges);
ClusterExchanges.Properties.VariableNames=[CellNumberPerCluster.Row;cellstr('e')];
tmp=strrep(multiCellModel.mets(UMetsID),'[u]','[e]');
ClusterExchanges.Row=tmp;

[~,Ia,Ib]=(intersect(model.mets,tmp));
% % % ClusterExchanges.Row(Ib)=model.metNames(Ia) ;
ClusterExchanges.Row(Ib)=model.mets(Ia) ;

Analysis=table(multiCellModel.rxns,FBAweighted.v, FBAweighted.x);
AnalysisKeep(clusterrouned).Analysis=Analysis;
ExchangesKeep(clusterrouned).Exchanges=ClusterExchanges;

end
end

