function[expanded_input_model]=compute_expanded_input_model(model,expanded_input_model,internal_mets,internal_rxns,intercluster_transport,mets_connected_to_ext_metabolite,umets,exmets,exRxns,number_of_cluster,fields)
% (c) Maria Pires Pacheco 2020 -University of Luxembourg

%Reconstruction of the S matrix for internal internal_rxns and internal
%metabolites of every cluster: internal_S_matrix (size (internal mets*
%numbers of clusters x  internal_rxns* number of clusters)
%The script starts computing the S internal matrix for 1 cluster: internal_S_matrix_for_1_cluster
% and places at the right places in the internal_S_matrix

% reconstructs the matrix for the internal metabolites moving between to the u compartment and the clusters : intercluster_trans_Smat
%(size (internal mets* numbers of clusters x   intercluster_transport * numbers of clusters)

% reconstructs the matrix for the u metabolite moving between to the u compartment and the clusters: intercluster_trans_Smat2
%(size (umets *numbers of clusters x   intercluster_transport* numbers of clusters)


% initializes an empty sparse matrix that has the correct size
internal_S_matrix=sparse(numel(internal_mets)*number_of_cluster,numel(internal_rxns)*number_of_cluster);
intercluster_trans_Smat=sparse(numel(internal_mets)*number_of_cluster,numel(intercluster_transport)*number_of_cluster);
intercluster_trans_Smat2=sparse(numel(umets),numel(intercluster_transport)*number_of_cluster);

% finds the indices of the internal metabolites and reactions
[~,internal_mets_index, ~]=intersect(model.mets,internal_mets);
[~,internal_rxns_index, ~]=intersect(model.rxns,internal_rxns);
[~,intercluster_transportID, ~]=intersect(model.rxns,intercluster_transport);

% builds a matrix with the internal reactions and metabolites of a sub-model
internal_S_matrix_for_1_cluster=model.S(internal_mets_index,internal_rxns_index);
internal_rxns_array=cell(numel(internal_rxns)*number_of_cluster,numel(fields)+1);
intercluster_trans_array=cell(numel(intercluster_transport)*number_of_cluster,numel(fields)+1);
internal_rxns_array_n=sparse(numel(internal_rxns)*number_of_cluster,numel(fields)+1);
intercluster_trans_array_n=sparse(numel(intercluster_transport)*number_of_cluster,numel(fields)+1);

[~, IB, IA]=intersect(model.mets,mets_connected_to_ext_metabolite);
[~, IC, ID]=intersect(model.rxns,intercluster_transport);
tmp(IA,ID)=model.S(IB,IC);
[~, IB, IA]=intersect(internal_mets,mets_connected_to_ext_metabolite);
intercluster_trans_Smat_4_1_cluster=sparse(numel(internal_mets),numel(intercluster_transport));

intercluster_trans_Smat_4_1_cluster(IB,:)=tmp(IA,:);
[~, IB, ~]=intersect(model.mets,exmets);
[~, IC, ~]=intersect(model.rxns,intercluster_transport);
intercluster_trans_Smat_4_1_cluster2=model.S(IB,IC);


missing_ex_mets= setdiff(model.mets(find(contains(model.mets,'[e]'))),exmets);
match=ismember(internal_mets,missing_ex_mets);
internal_mets(match)=strrep(internal_mets(match),'[e]','[u]');

% begin indices in the respective matrix
ind=1; %internal_S_matrix begin_dim1
ind2=size(internal_S_matrix_for_1_cluster,1);%internal_S_matrix end dim1
ind3=1;%internal_S_matrix begin dim2
ind4=size(internal_S_matrix_for_1_cluster,2);%internal_S_matrix end dim2
ind5=1; %intercluster_trans_Smat_4_1_cluster begin dim1 (ind6 == end)
ind7=1; %intercluster_trans_Smat_4_1_cluster begin dim2 (ind8 == end)
ind9=1;%intercluster_trans_Smat_4_1_cluster2 begin dim1 (ind10 == end)
%dim2 == ind7 and ind8


for i=1:number_of_cluster
    expanded_input_model.genes=[expanded_input_model.genes;strcat(model.genes,'_', num2str(i))];
    
    if i>1
        ind=ind2+1;
        ind3=ind4+1;
        ind5=ind6+1;
        ind7=ind8+1;
    end
    ind2=ind+size(internal_S_matrix_for_1_cluster,1)-1;
    ind4= ind3+ size(internal_S_matrix_for_1_cluster,2)-1;
    ind6=ind5+size(intercluster_trans_Smat_4_1_cluster,1)-1;
    ind8=ind7+size(intercluster_trans_Smat_4_1_cluster,2)-1;
    ind10=ind9+size(intercluster_trans_Smat_4_1_cluster2,1)-1;
    internal_S_matrix(ind: ind2,ind3:ind4)=internal_S_matrix_for_1_cluster;
    intercluster_trans_Smat(ind5:ind6,ind7:ind8)=intercluster_trans_Smat_4_1_cluster;
    intercluster_trans_Smat2(ind9:ind10,ind7:ind8)=intercluster_trans_Smat_4_1_cluster2;
    internal_rxns_array(ind3:ind4,1)=strcat(internal_rxns,'_', num2str(i));
    mets(ind: ind2)=strcat(internal_mets,'_', num2str(i));
    intercluster_trans_array(ind7:ind8,1)=(strcat(intercluster_transport,'_', num2str(i)));
   % do first metabolites
    for ii=1:numel(fields)
        if size( model.(fields{ii}),1)==size( model.mets,1)&& size( model.(fields{ii}),2)~=size( model.rxns,1)
            temp=model.(fields{ii});
           
            expanded_input_model.(fields{ii})= [expanded_input_model.(fields{ii});temp(internal_mets_index)];
        elseif size( model.(fields{ii}),1)==size( model.rxns,1) && size( model.(fields{ii}),1)~=size( model.mets,1)
            temp=model.(fields{ii});
            if strcmp(fields{ii},'rules') && i>1
                %Update the indices of the genes in the rules. The number of total genes is added to the indices for
                %each cluster (i.e 22 + numel(model.genes) * number of
                %clusters. This is done in two steps, the expression is
                %changed. Plus (+) appear in the expressions then it is
                %evaluated to remove the (+)
                
                temp=strrep(temp, 'x(', strcat( 'x(' ,num2str(numel(model.genes)*(i-1)),'+'));
                for counter1=1:numel(temp)
                    temp2=temp{counter1};
                    plus_pos=strfind(temp2,'+');
                    if ~isempty(plus_pos)
                        bracket1_pos=strfind(temp2,'(');
                        bracket2_pos=strfind(temp2,')');
                        if numel(plus_pos)==1 && numel(bracket2_pos)==1 && numel(bracket1_pos)==1
                            temp3=temp2(bracket1_pos+1: bracket2_pos-1);
                            temp4=eval(temp3);
                            temp{counter1}=strrep(temp2,temp3, num2str(temp4));
                        else
                            while ~isempty(plus_pos)
                                u22= plus_pos(1)-bracket1_pos;
                                u22=bracket1_pos(u22== min(u22(u22>0)));
                                u33= bracket2_pos-plus_pos(1);
                                u33=bracket2_pos(u33== min(u33(u33>0)));
                                temp3=temp2(u22+1: u33-1);
                                temp4=eval(temp3);
                                temp3a=strcat('(',temp3,')');
                                temp4a=strcat('(',num2str(temp4),')');
                                temp2=strrep(temp2,temp3a,temp4a);
                                plus_pos=strfind(temp2,'+');
                                bracket1_pos=strfind(temp2,'(');
                                bracket2_pos=strfind(temp2,')');
                            end
                            temp{counter1}=temp2;
                        end
                    end
                end
            else
            end
            if iscell(temp)
                internal_rxns_array(ind3:ind4,ii+1)=temp(internal_rxns_index);
                intercluster_trans_array(ind7:ind8,ii+1)=temp(intercluster_transportID);
            else
                internal_rxns_array_n(ind3:ind4,ii+1)=temp(internal_rxns_index);
                if strcmp(fields{ii}, 'lb')
                    
                    intercluster_trans_array_n(ind7:ind8,ii+1)=temp(intercluster_transportID);
                else
                    intercluster_trans_array_n(ind7:ind8,ii+1)=temp(intercluster_transportID);
                end
                
            end
            
            expanded_input_model.(fields{ii})= [expanded_input_model.(fields{ii});temp(internal_rxns_index)];
            
        end
    end
    
end


varname = regexp(exmets, '\[(.*?)\]', 'match', 'once');
expanded_input_model.mets=[mets';sort(strrep(exmets,varname,'[u]'));sort(strrep(exmets,varname,'[e]'))];
umets=strrep(exmets,varname,'[u]');
S=[internal_S_matrix,intercluster_trans_Smat];

S2=[sparse(size(umets,1),size(internal_S_matrix,2)),intercluster_trans_Smat2];
S=[S;S2];
S7=sparse(size(intercluster_trans_Smat,1), numel(exRxns));
S8=[S7;eye(size(intercluster_trans_Smat2,1))*-1];
S2=[S;sparse(size(intercluster_trans_Smat2,1),size(S,2))];
S=[S2, [S8;eye(size(intercluster_trans_Smat2,1))]];

S7=sparse(size(S,1), numel(exRxns));
S7(end+1-size(eye(size(intercluster_trans_Smat2,1)),1):end,:)=eye(size(intercluster_trans_Smat2,1));
S=[S, S7];


expanded_input_model.S=S;
trans=strcat('Trans_',sort(exmets));
Ex=strcat('Ex_',exmets);
dico_EX=cell(numel(Ex),2);
dico_EX(:,1)=Ex;
for i=1:numel(exmets)
    r=find(model.S(ismember(model.mets,exmets(i)),:));
    r=intersect(model.rxns(r),exRxns);
    dico_EX(i,2)=r;
end
expanded_input_model.rxns=[internal_rxns_array(:,1);intercluster_trans_array(:,1);trans;Ex];

for ii=1:numel(fields)
    if size( model.(fields{ii}),1)==size( model.rxns,1) && size( model.(fields{ii}),1)~=size( model.mets,1)
        if iscell(expanded_input_model.(fields{ii}))
            
            rul=[internal_rxns_array(:,ii+1);intercluster_trans_array(:,ii+1);cell(numel(trans),1);cell(numel(Ex),1)];
            rul(cellfun('isempty',rul))=cellstr('');
            expanded_input_model.(fields{ii})=rul;
            
        else
            if strcmp(fields{ii}, 'lb')
                
                %% lower bounds set -1000 by defaut for trans and Ex
                expanded_input_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*-1000; ones(numel(Ex),1)*-1000];
            else
                expanded_input_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*median(internal_rxns_array_n(:,ii+1)); ones(numel(Ex),1)*median(internal_rxns_array_n(:,ii+1))];
            end
        end
    else
        if iscell(expanded_input_model.(fields{ii}))
            [~,index]=ismember(exmets,model.mets);
            temp2=cell(numel(exmets),1);
            temp=model.(fields{ii});
            temp2(index>0)=temp(index(index>0));
            expanded_input_model.(fields{ii})=[expanded_input_model.(fields{ii});temp2;temp2];
        else
            expanded_input_model.(fields{ii})=[expanded_input_model.(fields{ii});ones(numel(umets),1)*mean(expanded_input_model.(fields{ii}));ones(numel(umets),1)*mean(expanded_input_model.(fields{ii}))];
        end
    end
end
'Hell0';
%% correct reversibility
lb=zeros(numel(exmets),1);
ub=zeros(numel(exmets),1);
grRules=cell(numel(exmets),1);
R_ex=cell(numel(exmets),1);
for i=1:numel(exmets)
    [~,r]=find(model.S(ismember(model.mets,exmets(i)),:));
    r=intersect(model.rxns(r),exRxns);
    R_ex(i,1)=r;
    
    sign=model.S(ismember(model.mets,exmets(i)),ismember(model.rxns,r));
    if full(sign) >0
        % then flip
        model.S(ismember(model.mets,exmets(i)),ismember(model.rxns,r))=-sign;
        tmp=model.ub(ismember(model.rxns,r));
        model.ub(ismember(model.rxns,r))=-model.lb(ismember(model.rxns,r));
        model.lb(ismember(model.rxns,r))=-tmp;
    end
end
[~,index]=ismember(R_ex,model.rxns);
lb(index>0)=model.lb(index(index>0));
ub(index>0)=model.ub(index(index>0));


Table_bound=table(R_ex,lb, ub, exmets,strcat('Ex_',exmets));
[~,index]=ismember(expanded_input_model.rxns,table2array(Table_bound(:,5)));
expanded_input_model.lb(index>0)=Table_bound.lb(index(index>0));
expanded_input_model.ub(index>0)=Table_bound.ub(index(index>0));

%%
expanded_input_model.rev= sparse(numel(expanded_input_model.rxns),1);
expanded_input_model.rev(expanded_input_model.lb<0)=1;

expanded_input_model = buildRxnGeneMat(expanded_input_model);
expanded_input_model = creategrRulesField(expanded_input_model);
expanded_input_model=fixIrr_rFASTCORMICS(expanded_input_model);

%% correct rules for exchanges
[~,index]=ismember(R_ex,model.rxns);
grRules(index>0)=model.grRules(index(index>0));
Table_bound=[Table_bound,table(grRules)];
rules_keep=cell(size(Table_bound,1),1);
for counter2=1:size(Table_bound,1)
    if ~cellfun('isempty',table2array(Table_bound(counter2,6)))
        
        
        for counter=1:number_of_cluster
            if counter==1
                rules=strcat(table2array(Table_bound(counter2,6)),'_', num2str(counter));
            else
                tmp= strcat(table2array(Table_bound(counter2,6)),'_',num2str(counter));
                rules=strcat(string(rules), ' or ','$ ', string(tmp));
            end
        end
        rules=strrep(rules,'$',' ');
        rules_keep(counter2,1)=cellstr(rules);
    end
    [~,index]=ismember(expanded_input_model.rxns,table2array(Table_bound(:,5)));
    

    expanded_input_model.grRules(index>0)=rules_keep(index(index>0),1);
end
    expanded_input_model.grRules(cellfun('isempty',expanded_input_model.grRules))=cellstr('');
for i=1:numel(expanded_input_model.grRules)
expanded_input_model.grRules(i)=strrep(expanded_input_model.grRules(i),'or',' or ');
     expanded_input_model.grRules(i)=strrep(expanded_input_model.grRules(i),'and',' and ');
end
    expanded_input_model = generateRules(expanded_input_model);


end

