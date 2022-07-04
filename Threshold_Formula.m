%% Define the threshold/formula elements - loading the information

%dis = [path,'\Discretization_Step'];

mkdir ([path,'\Formula']);

% Number of reactions that should or should not be there according to bulk
% data
bulk_rxns_should = size(readtable([dis,'\Table_Reaction_Should_Bulk.txt']),1);
bulk_rxns_should_not = size(readtable([dis,'\Table_Reaction_Should_not_Bulk.txt']),1);

% Number of reactions which should or not be there in composite model (Be int N & Bu int N reactions)
num_reactions_table = readcell([dis,'\Table_Reactions_Presence_Count.txt']);
num_reactions = num_reactions_table{2,1};
num_reactions_not = num_reactions_table{2,2};

model_composite_reactions = str2double(cellfun(@(S) S(1:end-2), model_composite.rxns, 'Uniform', 0));

% Number of genes which should or not be there, is taken from
% Metabolic_Gene_ID_Expressed/Never (Be & Bu genes)
Metabolic_Expressed_Genes = readmatrix([dis,'\Gene_ID_Expressed_Bulk_Data.txt']);
Metabolic_Never_Expressed_Genes = readmatrix([dis,'\Gene_ID_Never_Expressed_Bulk_Data.txt']);

% Metabolic genes in generic model according to bulk data (Be int N & Bu int N)
model_composite_genes = str2double(cellfun(@(S) S(1:end-2), model_composite.genes, 'Uniform', 0));
num_genes = numel(intersect(Metabolic_Expressed_Genes, model_composite_genes));
num_genes_not = numel(intersect(Metabolic_Never_Expressed_Genes, model_composite_genes));

% Core genes/reactions (C)
Core_Genes = readtable([dis,'\Table_Unique_Core_Genes.txt']);
Core_Reactions = readtable([dis,'\Table_Unique_Core_Reactions.txt']);

% % Genes (bulk data) intersect with context specific models
% Genes_Should_be = readtable([dis,'\Table_Number_Genes_Should_be_there.txt']);
% Genes_Should_not_be = readtable([dis,'\Table_Number_Genes_Should_not_be_there.txt']);
% 
% % Reactions (bulk data) intersect with context specific models
Reactions_Should_be = readtable([dis,'\Table_Number_of_reactions_should_be_there.txt']);
Reactions_Should_not_be = readtable([dis,'\Table_Number_of_reactions_should_not_be_there.txt']);

% Genes/Reactions in context specific models (A)
Gene_Associated = readtable([dis,'\Table_Unique_Gene_Asso.txt']);
Reactions_in_Context_Specific_model = readtable([dis,'\Table_Unique_Reactions_in_Context_Specific_Model.txt']);

% Genes/reactions intersect with composite model (Be int A & Bu int A)
Genes_int_present = readtable([dis,'\Table_intersect_genes_present.txt']);
Genes_int_not_present = readtable([dis,'\Table_intersect_genes_not_present.txt']);
Reactions_int_present = readtable([dis,'\Table_intersect_rxns_present.txt']);
Reactions_int_not_present = readtable([dis,'\Table_intersect_rxns_not_present.txt']);

%% Formula calculation
% Gene Formula

Gene_Formula_F1 = (Gene_Associated.Unique_Gene_Associated - Core_Genes.Unique_Core_Genes)./Gene_Associated.Unique_Gene_Associated ;
Gene_Formula_F2 = (num_genes - Genes_int_present.Var3)/num_genes ;
Gene_Formula_F3 = Genes_int_not_present.Var3/num_genes_not;
Whole_Gene_Formula = (Gene_Formula_F1 + Gene_Formula_F2 + Gene_Formula_F3)/3 ;

if any(~isfinite(Gene_Formula_F1)) == 1
    ind = find(~isfinite(Gene_Formula_F1));
    Whole_Gene_Formula(ind) = 1;
    Gene_Formula_F1(ind) = 1;
end

% Find the best threshold setting (Gene_Formula)

Index_Min_Gene_Formula = find(Whole_Gene_Formula == min(Whole_Gene_Formula));
Best_Cover_Threshold_Gene_Formula = Genes_int_present.Var1(Index_Min_Gene_Formula);
Best_Percentile_Threshold_Gene_Formula = Genes_int_present.Var2(Index_Min_Gene_Formula);

% Reaction Formula

Reaction_Formula_F1 = (Reactions_in_Context_Specific_model.Unique_Reactions_in_Context_Specific_Model - Core_Reactions.Unique_Core_Reactions)./Reactions_in_Context_Specific_model.Unique_Reactions_in_Context_Specific_Model; % Completness
Reaction_Formula_F2 = (num_reactions - Reactions_int_present.Var3)/num_reactions ; % Expressed rate
Reaction_Formula_F3 = Reactions_int_not_present.Var3/num_reactions_not ; % Not expressed rate
Whole_Reaction_Formula = (Reaction_Formula_F1 + Reaction_Formula_F2 + Reaction_Formula_F3)/3 ;

if (any(~isfinite(Reaction_Formula_F1)) == 1)
    ind1 = find(~isfinite(Reaction_Formula_F1));
    Whole_Reaction_Formula(ind) = 1;
    Reaction_Formula_F1(ind) = 1;
end

% Find the best threshold setting (Reaction_Formula)

Index_Min_Reaction_Formula = find(Whole_Reaction_Formula == min(Whole_Reaction_Formula));
Best_Cover_Threshold_Reaction_Formula = Reactions_Should_be.Cover(Index_Min_Reaction_Formula);
Best_Percentile_Threshold_Reaction_Formula = Reactions_Should_be.Percentile(Index_Min_Reaction_Formula);

best = {'Type','Percentile','Cover';
    'Genes',Best_Percentile_Threshold_Gene_Formula,Best_Cover_Threshold_Gene_Formula;
    'Reactions',Best_Percentile_Threshold_Reaction_Formula,Best_Cover_Threshold_Reaction_Formula};
writecell(best,[path,'\Formula\Table_Best_Thresholds.txt']);
%% Get an overview over the whole table
% Create a table with three column --> Cover, Percentile, Whole_Formula
% Saving the results --> nice data representation

Reaction_Formula_Table = array2table(Whole_Reaction_Formula);
Gene_Formula_Table = array2table (Whole_Gene_Formula);

Cover_table = array2table(Core_Genes.Cover);
Header_Cover = {'Cover'};
Cover_table.Properties.VariableNames = Header_Cover;

Percentile_table = array2table(Core_Genes.Percentile);
Header_Percentile = {'Percentile'};
Percentile_table.Properties.VariableNames = Header_Percentile;

Reaction_Formula_Table = [Cover_table,Percentile_table, Reaction_Formula_Table];
Gene_Formula_Table = [Cover_table,Percentile_table, Gene_Formula_Table];

writetable(Reaction_Formula_Table, [path,'\Formula\Table_Whole_Reaction_Formula.txt'])
writetable(Gene_Formula_Table,[path,'\Formula\Table_Whole_Gene_Formula.txt'])


%% Display of heatmaps of gene formula and components
if plot_formula == 1
    % As an input give a variable like : Whole_Gene_Formula ,Gene_Formula_F1 , Gene_Formula_F2 , Gene_Formula_F3
    % Preparing the input for heatmap
    
    % Start with a one column with 100 rows , each row is the result of one threshold setting
    % First Cover 5 and Percentile 5 , Last Cover 90 and Percentile 90
    
    Be_int_A_Genes = Genes_int_present.Var3;
    Be_Bu_int_N_Genes = [num_genes, num_genes_not];
    Bu_int_A_Genes = Genes_int_not_present.Var3;
    A_Genes = Gene_Associated.Unique_Gene_Associated;
    C_Genes = Core_Genes.Unique_Core_Genes;
    
    Gene_Formula = {Whole_Gene_Formula,Gene_Formula_F1,Gene_Formula_F2,Gene_Formula_F3, A_Genes, C_Genes, Be_int_A_Genes, Bu_int_A_Genes, Be_Bu_int_N_Genes};
    Gene_Titles = {'Whole Gene Formula','Gene Formula F1','Gene Formula F2','Gene Formula F3','A','C','Be int A', 'Bu int A', 'Be int N & Bu int N'};
    
    figure('WindowState','maximized')
    lngth = 20;
    Blue = [0, 0, 1];
    white = [1, 1, 1];
    colors_p = [linspace(white(1),Blue(1),lngth)', linspace(white(2),Blue(2),lngth)', linspace(white(3),Blue(3),lngth)'];
    
    for f=1:length(Gene_Formula)
        Formula_Cover_Percentage = Gene_Formula{f};
        
        if f < 5
            Cover_5  = Formula_Cover_Percentage(1:10,:); % Cover 5 and all the possible percentile values : 5,10,20,30...90
            Cover_10 = Formula_Cover_Percentage(11:20,:);
            Cover_20 = Formula_Cover_Percentage(21:30,:);
            Cover_30 = Formula_Cover_Percentage(31:40,:);
            Cover_40 = Formula_Cover_Percentage(41:50,:);
            Cover_50 = Formula_Cover_Percentage(51:60,:);
            Cover_60 = Formula_Cover_Percentage(61:70,:);
            Cover_70 = Formula_Cover_Percentage(71:80,:);
            Cover_80 = Formula_Cover_Percentage(81:90,:);
            Cover_90 = Formula_Cover_Percentage(91:100,:);
            %         Cover_05  = Formula_Cover_Percentage(1:10,:); % Cover 5 and all the possible percentile values : 5,10,20,30...90
            %         Cover_1 = Formula_Cover_Percentage(11:20,:);
            %         Cover_2 = Formula_Cover_Percentage(21:30,:);
            %         Cover_3 = Formula_Cover_Percentage(31:40,:);
            %         Cover_4 = Formula_Cover_Percentage(41:50,:);
            % Build a matrix where the rows are the covers and the columns are the percentiles
            %         Table_Formula = table(Cover_05, Cover_1, Cover_2, Cover_3, Cover_4);
            
            Table_Formula = table(Cover_5, Cover_10, Cover_20, Cover_30, Cover_40, Cover_50, Cover_60, Cover_70, Cover_80, Cover_90);
            Table_Formula_1 = table2array(Table_Formula);
            Table_Formula_2 = transpose(Table_Formula_1);
            
            subplot(3,3,f)
            I = imagesc(Table_Formula_2);
            caxis([0 1])
            title(['Heatmap of dataset, ', Gene_Titles{f}])
            ylabel('Coverage') %Cover
            xlabel('Percentile') %Percentage
            xmin = 5 ;
            xmax = 95 ;
            ymin = 0.1 ;
            ymax = 0.9 ;
            xticks([1 2 3 4 5 6 7 8 9 10 11]); % 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
            xticklabels({'5','10','20','30','40','50','60','70','80','90'}) %'0.1/70','0.1/80','0.1/90','0.1/95','0.2/70','0.2/80','0.2/90','0.2/95','0.3/70','0.3/80','0.3/90','0.3/95','0.4/70','0.4/80','0.4/90','0.4/95','0.5/70','0.5/80','0.5/90','0.5/95','0.6/70','0.6/80','0.6/90','0.6/95','0.7/70','0.7/80','0.7/90','0.7/95','0.8/70','0.8/80','0.8/90','0.8/95','0.9/70','0.9/80','0.9/90','0.9/95'})
            xtickangle(45)
            yticks([1 2 3 4 5 6 7 8 9 10])
            yticklabels({'0.0025','0.005','0.0075','0.01','0.0125','0.015','0.02','0.05','0.1','0.5'})
            colorbar;
            colormap(colors_p)
        end
        
        if f>4 && f<9
            Cover_5  = Formula_Cover_Percentage(1:10,:); % Cover 5 and all the possible percentile values : 5,10,20,30...90
            Cover_10 = Formula_Cover_Percentage(11:20,:);
            Cover_20 = Formula_Cover_Percentage(21:30,:);
            Cover_30 = Formula_Cover_Percentage(31:40,:);
            Cover_40 = Formula_Cover_Percentage(41:50,:);
            Cover_50 = Formula_Cover_Percentage(51:60,:);
            Cover_60 = Formula_Cover_Percentage(61:70,:);
            Cover_70 = Formula_Cover_Percentage(71:80,:);
            Cover_80 = Formula_Cover_Percentage(81:90,:);
            Cover_90 = Formula_Cover_Percentage(91:100,:);
            
            % Build a matrix where the rows are the covers and the columns are the percentiles
            
            Table_Formula = table(Cover_5, Cover_10, Cover_20, Cover_30, Cover_40, Cover_50, Cover_60, Cover_70, Cover_80, Cover_90);
            Table_Formula_1 = table2array(Table_Formula);
            Table_Formula_2 = transpose(Table_Formula_1);
            
            subplot(3,3,f)
            I = imagesc(Table_Formula_2);
            caxis([0 max(Gene_Formula{f})])
            title(['Heatmap of dataset, ', Gene_Titles{f}])
            ylabel('Coverage') %Cover
            xlabel('Percentile') %Percentage
            xmin = 5 ;
            xmax = 95 ;
            ymin = 0.1 ;
            ymax = 0.9 ;
            xticks([1 2 3 4 5 6 7 8 9 10 11]); % 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
            xticklabels({'5','10','20','30','40','50','60','70','80','90'}) %'0.1/70','0.1/80','0.1/90','0.1/95','0.2/70','0.2/80','0.2/90','0.2/95','0.3/70','0.3/80','0.3/90','0.3/95','0.4/70','0.4/80','0.4/90','0.4/95','0.5/70','0.5/80','0.5/90','0.5/95','0.6/70','0.6/80','0.6/90','0.6/95','0.7/70','0.7/80','0.7/90','0.7/95','0.8/70','0.8/80','0.8/90','0.8/95','0.9/70','0.9/80','0.9/90','0.9/95'})
            xtickangle(45)
            yticks([1 2 3 4 5 6 7 8 9 10])
            yticklabels({'0.0025','0.005','0.0075','0.01','0.0125','0.015','0.02','0.05','0.1','0.5'})
            colorbar;
            colormap(colors_p)
        end
        
        if f==9
            subplot(3,3,f)
            xlim([0 1])
            ylim([0 1])
            text(0.2, 0.65, ['Be int N: ',num2str(Be_Bu_int_N_Genes(1,1))],'FontSize',24)
            text(0.2, 0.45, ['Bu int N: ',num2str(Be_Bu_int_N_Genes(1,2))],'FontSize',24)
            axis off
        end
    end
    
    saveas(gcf,[path,'\Formula\Heatmap_Gene_Formula_Recap.png']);
    
    
    %% Display of heatmaps of reaction formula and components
    
    Be_int_A_Reactions = Reactions_int_present.Var3;
    Be_Bu_int_N_Reactions = [num_reactions, num_reactions_not];
    Bu_int_A_Reactions = Reactions_int_not_present.Var3;
    A_Reactions = Reactions_in_Context_Specific_model.Unique_Reactions_in_Context_Specific_Model;
    C_Reactions = Core_Reactions.Unique_Core_Reactions;
    
    Reaction_Formula = {Whole_Reaction_Formula,Reaction_Formula_F1,Reaction_Formula_F2,Reaction_Formula_F3, A_Reactions, C_Reactions, Be_int_A_Reactions, Bu_int_A_Reactions, Be_Bu_int_N_Reactions};
    Reaction_Titles = {'Whole Reaction Formula','Reaction Formula F1','Reaction Formula F2','Reaction Formula F3','A','C','Be int A', 'Bu int A', 'Be int N & Bu int N'};
    
    figure('WindowState','maximized')
    lngth = 20;
    Blue = [0, 0, 1];
    white = [1, 1, 1];
    colors_p = [linspace(white(1),Blue(1),lngth)', linspace(white(2),Blue(2),lngth)', linspace(white(3),Blue(3),lngth)'];
    
    for f=1:length(Reaction_Formula)
        Formula_Cover_Percentage = Reaction_Formula{f};
        
        if f < 5
            Cover_5  = Formula_Cover_Percentage(1:10,:); % Cover 5 and all the possible percentile values : 5,10,20,30...90
            Cover_10 = Formula_Cover_Percentage(11:20,:);
            Cover_20 = Formula_Cover_Percentage(21:30,:);
            Cover_30 = Formula_Cover_Percentage(31:40,:);
            Cover_40 = Formula_Cover_Percentage(41:50,:);
            Cover_50 = Formula_Cover_Percentage(51:60,:);
            Cover_60 = Formula_Cover_Percentage(61:70,:);
            Cover_70 = Formula_Cover_Percentage(71:80,:);
            Cover_80 = Formula_Cover_Percentage(81:90,:);
            Cover_90 = Formula_Cover_Percentage(91:100,:);
            
            % Build a matrix where the rows are the covers and the columns are the percentiles
            
            Table_Formula = table(Cover_5, Cover_10, Cover_20, Cover_30, Cover_40, Cover_50, Cover_60, Cover_70, Cover_80, Cover_90);
            Table_Formula_1 = table2array(Table_Formula);
            Table_Formula_2 = transpose(Table_Formula_1);
            
            subplot(3,3,f)
            I = imagesc(Table_Formula_2);
            caxis([0 1])
            title(['Heatmap of dataset, ', Reaction_Titles{f}])
            ylabel('Coverage') %Cover
            xlabel('Percentile') %Percentage
            xmin = 5 ;
            xmax = 95 ;
            ymin = 0.1 ;
            ymax = 0.9 ;
            xticks([1 2 3 4 5 6 7 8 9 10 11]); % 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40])
            xticklabels({'5','10','20','30','40','50','60','70','80','90'}) %'0.1/70','0.1/80','0.1/90','0.1/95','0.2/70','0.2/80','0.2/90','0.2/95','0.3/70','0.3/80','0.3/90','0.3/95','0.4/70','0.4/80','0.4/90','0.4/95','0.5/70','0.5/80','0.5/90','0.5/95','0.6/70','0.6/80','0.6/90','0.6/95','0.7/70','0.7/80','0.7/90','0.7/95','0.8/70','0.8/80','0.8/90','0.8/95','0.9/70','0.9/80','0.9/90','0.9/95'})
            xtickangle(45)
            yticks([1 2 3 4 5 6 7 8 9 10])
            yticklabels({'0.0025','0.005','0.0075','0.01','0.0125','0.015','0.02','0.05','0.1','0.5'})
            colorbar;
            colormap(colors_p)
        end
        
        if f > 4 && f < 9
            Cover_5  = Formula_Cover_Percentage(1:10,:); % Cover 5 and all the possible percentile values : 5,10,20,30...90
            Cover_10 = Formula_Cover_Percentage(11:20,:);
            Cover_20 = Formula_Cover_Percentage(21:30,:);
            Cover_30 = Formula_Cover_Percentage(31:40,:);
            Cover_40 = Formula_Cover_Percentage(41:50,:);
            Cover_50 = Formula_Cover_Percentage(51:60,:);
            Cover_60 = Formula_Cover_Percentage(61:70,:);
            Cover_70 = Formula_Cover_Percentage(71:80,:);
            Cover_80 = Formula_Cover_Percentage(81:90,:);
            Cover_90 = Formula_Cover_Percentage(91:100,:);
            
            % Build a matrix where the rows are the covers and the columns are the percentiles
            
            Table_Formula = table(Cover_5, Cover_10, Cover_20, Cover_30, Cover_40, Cover_50, Cover_60, Cover_70, Cover_80, Cover_90);
            Table_Formula_1 = table2array(Table_Formula);
            Table_Formula_2 = transpose(Table_Formula_1);
            
            subplot(3,3,f)
            I = imagesc(Table_Formula_2);
            caxis([0 max(Reaction_Formula{f})])
            title(['Heatmap of dataset, ', Gene_Titles{f}])
            ylabel('Coverage') %Cover
            xlabel('Percentile') %Percentage
            xmin = 5 ;
            xmax = 95 ;
            ymin = 0.1 ;
            ymax = 0.9 ;
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
            text(0.2, 0.65, ['Be int N: ',num2str(Be_Bu_int_N_Reactions(1,1))],'FontSize',24)
            text(0.2, 0.45, ['Bu int N: ',num2str(Be_Bu_int_N_Reactions(1,2))],'FontSize',24)
            axis off
        end
    end
    
    saveas(gcf,[path,'\Formula\Heatmap_Reaction_Formula_Recap.png']);
end
