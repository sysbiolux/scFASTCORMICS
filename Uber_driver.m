% Driver script for HPC containing all 20 data sets
feature astheightlimit 2000

changeCobraSolver('ibm_cplex')
% Data_1
tic
feature astheightlimit 2000
user_path = 'C:\Users\0170676825\Desktop\sc_transfer';
scdataset = dir([user_path,'\Datasets\Data_1\*.txt']);
set_name = 'Data_1_model_orig';
Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_tumor.txt']);
Driver_script;
toc
clear; clc;

% Data_2
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_2\*.txt']);
set_name = 'Data_2_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_CRC_NM.txt']);
Driver_script;
clear; clc;

% Data_3
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_3\*.txt']);
set_name = 'Data_3_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_4
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_4\*.txt']);
set_name = 'Data_4_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_5
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_5\*.txt']);
set_name = 'Data_5_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_6
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_6\*.txt']);
set_name = 'Data_6_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_7
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_7\*.txt']);
set_name = 'Data_7_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_8
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_8\*.txt']);
set_name = 'Data_8_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_9
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_9\*.txt']);
set_name = 'Data_9_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_10
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_10\*.txt']);
set_name = 'Data_10_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_11
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_11\*.txt']);
set_name = 'Data_11_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_12
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_12\*.txt']);
set_name = 'Data_12_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_13
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_13\*.txt']);
set_name = 'Data_13_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_14
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_14\*.txt']);
set_name = 'Data_14_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_hiPancreas.txt']);
Driver_script;
clear; clc;

% Data_15
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_15\*.txt']);
set_name = 'Data_15_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_PACA.txt']);
Driver_script;
clear; clc;

% Data_16
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_16\*.txt']);
set_name = 'Data_16_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_PACA.txt']);
Driver_script;
clear; clc;

% Data_17
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_17\*.txt']);
set_name = 'Data_17_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_PBMC.txt']);
Driver_script;
clear; clc;

% Data_18
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_18\*.txt']);
set_name = 'Data_18_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_Liver_Cancer.txt']);
Driver_script;
clear; clc;

% Data_19
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_19\*.txt']);
set_name = 'Data_19_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_Liver_Cancer.txt']);
Driver_script;
clear; clc;

% Data_20
feature astheightlimit 2000
user_path = 'D:\scFASTCORMICS';
scdataset = dir([user_path,'\Datasets\Data_20\*.txt']);
set_name = 'Data_20_model';
Discretization_Table = readtable([user_path,'\Discretization_Table_Breast_Cancer.txt']);
Driver_script;
clear; clc;
