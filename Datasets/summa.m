clear all
s=zeros(20,17) 

for i=1:20
    i
    name=strcat('Data_',num2str(i));
 scdataset = dir(['Data_', num2str(i),'\*.txt']);


for ii=1:size(scdataset,1)
  data = readtable([scdataset(ii).folder,'\', scdataset(ii).name]);
  s(i,ii)=size(data (:,4:end),2);
end
 end