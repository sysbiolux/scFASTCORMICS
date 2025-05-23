function[input_data]=Numbering(scdataset, path,printLevel, cell_index )
 
 
%% Loading/Importing the data
if printLevel==1
    mkdir ([path,'/Numbered'])
end
files_length = length(scdataset);

input_data=struct();

    
for i=1:files_length    
    %import data
    data = readtable([scdataset(i).folder,'/', scdataset(i).name]);
    if ~isnumeric(data.GeneID)
        data.GeneID=str2double(data.GeneID);
    end
    
    %add 'Gene Numeration' column to the right position, after Gene ID
    gene_numeration = cellstr(convertStringsToChars(strcat(num2str(data.GeneID) ,'_',num2str(i))));
    gene_numeration = strrep(gene_numeration,' ','');
    
    %build the table
    table_data = [data(:,2),gene_numeration, data(:,4:end)];
    table_data.Properties.VariableNames(2) = {'Gene_Numeration'};
    if nargin==4
      table_data=table_data(:,cell_index)

    else

    end
    
    %save the table
    if printLevel==1
        name = [path,'/Numbered/Numbered_', scdataset(i).name];
        writetable(table_data,name);
    end
    input_data(i).input_data=table_data;
    
end