function[selectedData]=readtableMaria(filename)
% Open the file
fileID = fopen(filename, 'r');

% Read the header line to determine column indices
headerLine = fgetl(fileID);
headers = strsplit(headerLine, ',');

% Specify the columns you want to read
columnsToRead = headers(1:1000);
colIndices = cellfun(@(col) find(strcmp(headers, col)), columnsToRead);

% Initialize a cell array to store the data
data = cell(1, numel(columnsToRead));

% Read the file line by line
line = fgetl(fileID);
while ischar(line)
    % Split the line into columns
    cols = strsplit(line, ',');
    
    % Extract the specified columns
    for i = 1:numel(colIndices)
        data{i} = [data{i}; cols{colIndices(i)}]; %#ok<AGROW>
    end
    
    % Read the next line
    line = fgetl(fileID);
end

% Close the file
fclose(fileID);

% Convert cell arrays to numeric arrays if necessary
for i = 1:numel(data)
    data{i} = str2double(data{i});
end

% Combine the columns into a matrix if needed
selectedData = [data{:}];

% Display the data
