%% Loading/Importing the data

mkdir ([path,'\Numbered'])

for i=1:files_length
    
    %import data
    data = readtable([scdataset(i).folder,'\', scdataset(i).name]);
    
    %add 'Gene Numeration' column to the right position, after Gene ID
    gene_numeration = cellstr(convertStringsToChars(strcat(num2str(data.GeneID) ,'_',num2str(i))));
    gene_numeration = strrep(gene_numeration,' ','');
    
    %build the table
    table_data = [data(:,2),gene_numeration, data(:,4:end)];
    table_data.Properties.VariableNames(2) = {'Gene_Numeration'};
    
    %save the table
    name = [path,'\Numbered\Numbered_', scdataset(i).name];
    writetable(table_data,name);
end