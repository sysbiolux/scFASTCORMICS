%% Creating Jaccard similarity plot of reaction count between clusters of each threshold condition

mkdir ([extract,'\Jaccard_Similarity']);
cb_keep = [0.9200, 0.1100, 0.0267, 0.8150];
for c=1:numel(Cover)
    for p=1:numel(percent)
        
        if isfile([extract,'\Reaction_Matrix\Reaction_Matrix_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.txt']) == 1
            Reactions = readtable([extract,'\Reaction_Matrix\Reaction_Matrix_Cover_',num2str(Cover(c)),'_Percentile_',num2str(percent(p)),'.txt']);
            Reactions.Properties.VariableNames = Headers;
            Reaction_Matrix = table2array(Reactions(:,2:end));
            
            J = squareform(pdist(Reaction_Matrix','jaccard'));
            if numel(unique(J)) > 1
                if sum(sum(isnan(J)))>0
                    altcolor = [0 0 255; 255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
                        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
                else
                    altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
                        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
                end
                
                J(isnan(J)) = 1.1;
                
                %get and save the default size
                defaultPosition = get(0,'DefaultFigurePosition');
                %get the current screen size
                screensize = get(groot, 'Screensize');
                %screensize = get(0, 'Screensize'); %for earlier Matlab versions (e.g. Matlab 2010)
                %set default figure position to full screen
                screensize = [1,1,800,800];
                set(0, 'DefaultFigurePosition', screensize);
                
                data_values_names = Reactions.Properties.VariableNames(2:end)';
                
                cgo_J = clustergram(1-J,...
                    'RowLabels', data_values_names,...
                    'ColumnLabels', data_values_names,...
                    'ColumnLabelsRotate',340, ...
                    'Cluster', 'all', ...
                    'Annotate', 'true',...
                    'symmetric','False',...
                    'AnnotColor','k',...
                    'Colormap', altcolor)
                %addTitle(cgo_J,{'Model similarity using Jaccard index',' based on the reconstructed models'' reactions', ['Cover ', num2str(Cover(c)) ,' and Percentage ', num2str(percent(p))]});
                addTitle(cgo_J,{'Cluster similarity using Jaccard index',' based on the reconstructed models'' reactions'});
                                
                plot(cgo_J);
                % cb = colorbar;
                %cb.Position = cb_keep;

                figureHandle = gcf;
                %# make all text in the figure to size 14 and bold
                fig_gcf = findall(0,'type','figure', 'tag', 'Clustergram');
                set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
                set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)
                saveas(gcf,[extract,'\Jaccard_Similarity\1Jaccard_C_',num2str(Cover(c)),'_P_',num2str(percent(p)),'.png']);
            end
            J = 1-J;
            rowNames = Headers(2:end);
            colNames = Headers(2:end);
            Jac_table = array2table(J,'RowNames',rowNames,'VariableNames',colNames);
            
            writetable(Jac_table,[extract,'\Jaccard_Similarity\1Jaccard_Similarity_Index_C_',num2str(Cover(c)),'_P_',num2str(percent(p)),'.txt'], 'WriteRowNames', true);
        end
    end
end
close all force
%'ColumnLabels', data_values_names,...