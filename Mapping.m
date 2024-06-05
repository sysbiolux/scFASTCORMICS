function[input_data]=Mapping(input_data, model_composite, path, printLevel, scdataset)
%% Mapping the genes from the sc-RNA seq data to the composite model
if printLevel==1
    %Creates a folder named Mapped
    mkdir ([path,'/Mapped']);
end

files_length = length(input_data);


for i=1:files_length
    
    table_numbered= input_data(i).input_data;
    gene_numeration = cellfun(@convertStringsToChars,table2cell(table_numbered(:,2)),'UniformOutput',false);
    
    intersect_genes = intersect(gene_numeration, model_composite.genes);
    intersect_indice = (ismember(cellfun(@convertStringsToChars,table2cell(table_numbered(:,2)),'UniformOutput',false), intersect_genes));
    table_final = table_numbered(intersect_indice,2:end);
    
    name = [path,'/Mapped/Mapped_',scdataset(i).name];
    input_data(i).table_final=table_final;
    if printLevel==1
        writetable(table_final, name);
    end    
end


