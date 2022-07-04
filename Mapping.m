%% Mapping the genes from the sc-RNA data to the composite model

mkdir ([path,'\Mapped']);

for i=1:files_length
        
    table_numbered = readtable([path,'\Numbered\Numbered_',scdataset(i).name]);
        
    gene_numeration = cellfun(@convertStringsToChars,table2cell(table_numbered(:,2)),'UniformOutput',false);
    
    intersect_genes = intersect(gene_numeration, model_composite.genes);
    intersect_indice = find(ismember(cellfun(@convertStringsToChars,table2cell(table_numbered(:,2)),'UniformOutput',false), intersect_genes));
    table_final = table_numbered(intersect_indice,2:end);
    
    name = [path,'\Mapped\Mapped_',scdataset(i).name];
    writetable(table_final, name);
    
end


