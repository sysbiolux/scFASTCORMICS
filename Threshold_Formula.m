function[best]=Threshold_Formula(model_composite,path, T, plot_formula, Cover_range, REI_range, Optimization_global)
%% Define the threshold/formula elements - loading the information


%dis = [path,'\Discretization_Step'];

mkdir ([path,'\Formula']);

% Number of reactions that should or should not be there according to bulk
% data
% bulk_rxns_should = T.rxns_should_bulk_data;
% bulk_rxns_should_not =T.rxns_should_not_bulk_data;

% Number of reactions which should or not be there in composite model (Be int N & Bu int N reactions)

model_composite_reactions = (cellfun(@(S) S(1:end-2), model_composite.rxns, 'Uniform', 0));


% Metabolic genes in generic model according to bulk data (Be int N & Bu int N)
%model_composite_genes = str2double(cellfun(@(S) S(1:end-2), model_composite.genes, 'Uniform', 0));


% Core genes/reactions (C)

% % Reactions (bulk data) intersect with context specific models
Reactions_Should_be = T.Table_Reaction_Expressed;

% Genes/Reactions in context specific models (A)
A = sum(Optimization_global(1).A,1)';
C = sum(Optimization_global(1).C_keep,1)';

BE_inter_N=numel(intersect(model_composite_reactions,T(1).Rxns_Present));
BU_inter_N=numel(intersect(model_composite_reactions,T(1).Rxns_Not_Present));
BE_inter_A= table2array(T.Table_Reaction_Expressed(:,3));
BU_inter_A= table2array(T.Table_Reaction_Never_Expressed(:,3));


%% Formula calculation% Reaction Formula

Reaction_Formula_F1 =(A-C)./A; % Completness
Reaction_Formula_F2 = (BE_inter_N-(BE_inter_A))./BE_inter_N;
Reaction_Formula_F3 = BU_inter_A./BU_inter_N;
Whole_Reaction_Formula = (Reaction_Formula_F1 + Reaction_Formula_F2 + Reaction_Formula_F3)/3 ;

if any(~isfinite(Reaction_Formula_F1))
    ind = find(~isfinite(Reaction_Formula_F1));
    Whole_Reaction_Formula(ind) = 1;
    Reaction_Formula_F1(ind) = 1;
end

% Find the best threshold setting (Reaction_Formula)

Index_Min_Reaction_Formula = find(Whole_Reaction_Formula == min(Whole_Reaction_Formula));
Best_Cover_Threshold_Reaction_Formula = Reactions_Should_be.Cover(Index_Min_Reaction_Formula);
Best_REI_Threshold_Reaction_Formula = Reactions_Should_be.Percentile(Index_Min_Reaction_Formula);
best.Best_REI_Threshold_Reaction_Formula=Best_REI_Threshold_Reaction_Formula;
best.Best_Cover_Threshold_Reaction_Formula=Best_Cover_Threshold_Reaction_Formula;
best.Whole_Reaction_Formula=Whole_Reaction_Formula;


%writecell(best,[path,'\Formula\Table_Best_Thresholds.txt']);
%% Get an overview over the whole table
% Create a table with three column --> Cover, Percentile, Whole_Formula
% Saving the results --> nice data representation

%% Display of heatmaps of gene formula and components
if plot_formula == 1
    % As an input give a variable like : Whole_Gene_Formula ,Gene_Formula_F1 , Gene_Formula_F2 , Gene_Formula_F3
    % Preparing the input for heatmap

    
    
    %% Display of heatmaps of reaction formula and components
    
   A_Reactions = A;
     C_Reactions = C;
%     
     Reaction_Formula = {Whole_Reaction_Formula,Reaction_Formula_F1,Reaction_Formula_F2,Reaction_Formula_F3, A_Reactions, C_Reactions, BE_inter_A, BU_inter_A, BE_inter_N, BU_inter_A};
     Reaction_Titles = {'Whole Reaction Formula','Reaction Formula F1','Reaction Formula F2','Reaction Formula F3','A','C','Be int A', 'Bu int A', 'Be int N' ,'Bu int N'};
%     
    figure('WindowState','maximized')
    lngth = 20;
    Blue = [0, 0, 1];
    white = [1, 1, 1];
    colors_p = [linspace(white(1),Blue(1),lngth)', linspace(white(2),Blue(2),lngth)', linspace(white(3),Blue(3),lngth)'];
    
    for f=1:length(Reaction_Formula)
        Formula_Cover_Percentage = Reaction_Formula{f};
        if f < 5
           Table_Formula_1=reshape(Formula_Cover_Percentage, numel(Cover_range), numel(REI_range));
           Table_Formula_2 = transpose(Table_Formula_1);
           
            subplot(3,3,f);
            imagesc(Table_Formula_2);
            caxis([0 1])
            title(['Heatmap of dataset, ', Reaction_Titles{f}])
            ylabel('Coverage') %Cover
            xlabel('Percentile') %Percentage
%             xmin = 5 ;
%             xmax = 95 ;
%             ymin = 0.1 ;
%             ymax = 0.9 ;
            xticks([1 2 3 4 5 6 7 8 9 10 11]); % 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
            xticklabels({'5','10','20','30','40','50','60','70','80','90'}) %'0.1/70','0.1/80','0.1/90','0.1/95','0.2/70','0.2/80','0.2/90','0.2/95','0.3/70','0.3/80','0.3/90','0.3/95','0.4/70','0.4/80','0.4/90','0.4/95','0.5/70','0.5/80','0.5/90','0.5/95','0.6/70','0.6/80','0.6/90','0.6/95','0.7/70','0.7/80','0.7/90','0.7/95','0.8/70','0.8/80','0.8/90','0.8/95','0.9/70','0.9/80','0.9/90','0.9/95'})
            xtickangle(45)
            yticks([1 2 3 4 5 6 7 8 9 10])
            yticklabels({'0.0025','0.005','0.0075','0.01','0.0125','0.015','0.02','0.05','0.1','0.5'})
            colorbar;
            colormap(colors_p)
        end
        
        if f > 4 && f < 9
            
            % Build a matrix where the rows are the covers and the columns are the percentiles
            
            Table_Formula_1=reshape(Formula_Cover_Percentage, numel(Cover_range), numel(REI_range));
            Table_Formula_2 = transpose(Table_Formula_1);
            
            subplot(3,3,f)
            imagesc(Table_Formula_2);
            caxis([0 max(Reaction_Formula{f})])
            ylabel('Coverage') %Cover
            xlabel('REI') %Percentage
            title(['Heatmap of dataset, ', Reaction_Titles{f}])

%             xmin = 5 ;
%             xmax = 95 ;
%             ymin = 0.1 ;
%             ymax = 0.9 ;
            xticks([1 2 3 4 5 6 7 8 9 10 11]); % 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
            xticklabels({'5','10','20','30','40','50','60','70','80','90'}) %'0.1/70','0.1/80','0.1/90','0.1/95','0.2/70','0.2/80','0.2/90','0.2/95','0.3/70','0.3/80','0.3/90','0.3/95','0.4/70','0.4/80','0.4/90','0.4/95','0.5/70','0.5/80','0.5/90','0.5/95','0.6/70','0.6/80','0.6/90','0.6/95','0.7/70','0.7/80','0.7/90','0.7/95','0.8/70','0.8/80','0.8/90','0.8/95','0.9/70','0.9/80','0.9/90','0.9/95'})
            xtickangle(45)
            yticks([1 2 3 4 5 6 7 8 9 10])
            yticklabels({'0.0025','0.005','0.0075','0.01','0.0125','0.015','0.02','0.05','0.1','0.5'})
            colorbar;
            colormap(colors_p)
            
        end
        
        if f == 9
            subplot(3,3,f)
            xlim([0 1])
            ylim([0 1])
            text(0.2, 0.65, ['Be int N: ',num2str(BE_inter_N)],'FontSize',24)
            text(0.2, 0.45, ['Bu int N: ',num2str(BU_inter_N)],'FontSize',24)
            axis off
        end
    end
    
    saveas(gcf,[path,'\Formula\Heatmap_Reaction_Formula_Recap.png']);
end
