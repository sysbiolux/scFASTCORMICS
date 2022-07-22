% Driver script for HPC containing all 20 data sets
%starts 22/07 12h49
clear all
changeCobraSolver('ibm_cplex')
user_path = 'C:\Users\maria.pacheco\Desktop\scFASTCORMICS';

% Data_1
biodbnet = readtable([user_path,'\dico_biodbnet.txt']); % The dictionary should map the model genes id to the data ids
Best_keep=zeros(20,2);
load('simpleRev_recon3.mat','model');
TIME=zeros(20,1);

for ii=1:20
tic
if ii==1
    Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_tumor.txt']);
elseif ii==2
    Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_NM.txt']);
elseif ii> 3 && ii<15
    Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
elseif ii==15 || ii==16
    Discretization_Table = readtable([user_path,'\Discretization_Table_PACA.txt']);
elseif ii== 17
    Discretization_Table = readtable([user_path,'\Discretization_Table_PBMC.txt']);

elseif ii==18 || ii==19
    Discretization_Table = readtable([user_path,'\Discretization_Table_Liver_Cancer.txt']);
elseif ii>19
    Discretization_Table_Breast_Cancer
end
scdataset = dir([user_path,'\Datasets\Data_', num2str(ii),'\*.txt']);
set_name = strcat('Data_',num2str(ii)','_model_orig');
  
[best]=Driver_script_M(Discretization_Table, set_name, scdataset, biodbnet,user_path, model);
Best_keep(ii,1)=best.Best_REI_Threshold_Reaction_Formula;
Best_keep(ii,2)=best.Best_Cover_Threshold_Reaction_Formula;
TIME(ii)=toc;
end