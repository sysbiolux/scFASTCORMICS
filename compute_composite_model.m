function[composite_model]=compute_composite_model(model,composite_model,internal_mets,internal_rxns,intercluster_transport,internal_2_intersticiel_mets,intersticiel_mets,exmets,exRxns,number_of_cluster,fields)
% (c) Maria Pires Pacheco 2020 -University of Luxembourg


%reconstruct the Smatrix for internal internal_rxns and intercluster_transport
internalS_matrix=sparse(numel(internal_mets)*number_of_cluster,numel(internal_rxns)*number_of_cluster);
intercluster_trans_Smat=sparse(numel(internal_mets)*number_of_cluster,numel(intercluster_transport)*number_of_cluster);
intercluster_trans_Smat2=sparse(numel(intersticiel_mets),numel(intercluster_transport)*number_of_cluster);
[~,internal_metsID, ~]=intersect(model.mets,internal_mets);
[~,internal_rxnsID, ~]=intersect(model.rxns,internal_rxns);
[~,intercluster_transportID, ~]=intersect(model.rxns,intercluster_transport);

internalS_matrix_4_1_cluster=model.S(internal_metsID,internal_rxnsID);
internal_rxns_array=cell(numel(internal_rxns)*number_of_cluster,numel(fields)+1);
intercluster_trans_array=cell(numel(intercluster_transport)*number_of_cluster,numel(fields)+1);
internal_rxns_array_n=sparse(numel(internal_rxns)*number_of_cluster,numel(fields)+1);
intercluster_trans_array_n=sparse(numel(intercluster_transport)*number_of_cluster,numel(fields)+1);

[~, IB, IA]=intersect(model.mets,internal_2_intersticiel_mets);
[~, IC, ID]=intersect(model.rxns,intercluster_transport);
tmp(IA,ID)=model.S(IB,IC);
[~, IB, IA]=intersect(internal_mets,internal_2_intersticiel_mets);
intercluster_trans_Smat_4_1_cluster=sparse(numel(internal_mets),numel(intercluster_transport));

intercluster_trans_Smat_4_1_cluster(IB,:)=tmp(IA,:);
[~, IB, ~]=intersect(model.mets,exmets);
[~, IC, ~]=intersect(model.rxns,intercluster_transport);
intercluster_trans_Smat_4_1_cluster2=model.S(IB,IC);

 
    missing_ex_mets= setdiff(model.mets(find(contains(model.mets,'[e]'))),exmets);
    match=ismember(internal_mets,missing_ex_mets);
    internal_mets(match)=strrep(internal_mets(match),'[e]','[u]');

x=1;
y=1;
xx=1;
xxx=1;

yy=1;
x2=size(internalS_matrix_4_1_cluster,1);
y2=size(internalS_matrix_4_1_cluster,2);


for i=1: (number_of_cluster)
    %%
    composite_model.genes=[composite_model.genes;strcat(model.genes,'_', num2str(i))];
    
    if i>1
        x=x2+1;
        y=y2+1;
        xx=xx2+1;
        yy=yy2+1;
        
        
    end
    x2=x+size(internalS_matrix_4_1_cluster,1)-1;
    y2= y+ size(internalS_matrix_4_1_cluster,2)-1;
    xx2=xx+size(intercluster_trans_Smat_4_1_cluster,1)-1;
    yy2=yy+size(intercluster_trans_Smat_4_1_cluster,2)-1;
    xxx2=xxx+size(intercluster_trans_Smat_4_1_cluster2,1)-1;
    
    internalS_matrix(x: x2,y:y2)=internalS_matrix_4_1_cluster;
    intercluster_trans_Smat(xx:xx2,yy:yy2)=intercluster_trans_Smat_4_1_cluster;
    intercluster_trans_Smat2(xxx:xxx2,yy:yy2)=intercluster_trans_Smat_4_1_cluster2;
    
    internal_rxns_array(y:y2,1)=strcat(internal_rxns,'_', num2str(i));
    
    
    mets(x: x2)=strcat(internal_mets,'_', num2str(i));
   
 
    
    intercluster_trans_array(yy:yy2,1)=(strcat(intercluster_transport,'_', num2str(i)));
    for ii=1:numel(fields)
      
        
        if size( model.(fields{ii}),1)==size( model.mets,1)&& size( model.(fields{ii}),2)~=size( model.rxns,1)
            
            temp=model.(fields{ii});
            
            composite_model.(fields{ii})= [composite_model.(fields{ii});temp(internal_metsID)];
            
            
        elseif size( model.(fields{ii}),1)==size( model.rxns,1) && size( model.(fields{ii}),1)~=size( model.mets,1)
            temp=model.(fields{ii});
            if strcmp(fields{ii},'rules') && i>1
               
                temp=strrep(temp, 'x(', strcat( 'x(' ,num2str(numel(model.genes)*(i-1)),'+'));
                %%
                for ooo=1:numel(temp)
                    temp2=temp{ooo};
                    
                    
                    
                    u=strfind(temp2,'+');
                    if ~isempty(u)
                        
                        u2=strfind(temp2,'(');
                        u3=strfind(temp2,')');
                        if numel(u)==1 & numel(u3)==1 & numel(u2)==1
                            temp3=temp2(u2+1: u3-1);
                            temp4=eval(temp3);
                            temp{ooo}=strrep(temp2,temp3, num2str(temp4));
                        else
                            while ~isempty(u)
                                
                                io=1;
                                u22= u(io)-u2;
                                u22=u2(u22== min(u22(u22>0)));
                                u33= u3-u(io);
                                u33=u3(u33== min(u33(u33>0)));
                                
                                temp3=temp2(u22+1: u33-1);
                                temp4=eval(temp3);
                                temp3a=strcat('(',temp3,')');
                                temp4a=strcat('(',num2str(temp4),')');

                                temp2=strrep(temp2,temp3a,temp4a);
                                u=strfind(temp2,'+');
                                u2=strfind(temp2,'(');
                                u3=strfind(temp2,')');
                            end
                            temp{ooo}=temp2;
                            
                        end
                    end
                end
                
                
            else
            end
            if iscell(temp)
                
                internal_rxns_array(y:y2,ii+1)=temp(internal_rxnsID);
                intercluster_trans_array(yy:yy2,ii+1)=temp(intercluster_transportID);
            else
                internal_rxns_array_n(y:y2,ii+1)=temp(internal_rxnsID);
                if strcmp(fields{ii}, 'lb')
                    
                    intercluster_trans_array_n(yy:yy2,ii+1)=temp(intercluster_transportID);
                else
                    intercluster_trans_array_n(yy:yy2,ii+1)=temp(intercluster_transportID);
                end
                
            end
            
            composite_model.(fields{ii})= [composite_model.(fields{ii});temp(internal_rxnsID)];
            
        end
    end
    
end


varname = regexp(exmets, '\[(.*?)\]', 'match', 'once');
composite_model.mets=[mets';sort(strrep(exmets,varname,'[u]'));sort(strrep(exmets,varname,'[e]'))];
a=strrep(exmets,varname,'[u]');
S=[internalS_matrix,intercluster_trans_Smat];

S2=[sparse(size(intersticiel_mets,1),size(internalS_matrix,2)),intercluster_trans_Smat2];
S=[S;S2];
S7=sparse(size(intercluster_trans_Smat,1), numel(exRxns));
S8=[S7;eye(size(intercluster_trans_Smat2,1))*-1];
S2=[S;sparse(size(intercluster_trans_Smat2,1),size(S,2))];
S=[S2, [S8;eye(size(intercluster_trans_Smat2,1))]];

S7=sparse(size(S,1), numel(exRxns));
S7(end+1-size(eye(size(intercluster_trans_Smat2,1)),1):end,:)=eye(size(intercluster_trans_Smat2,1));
S=[S, S7];


composite_model.S=S;
trans=strcat('Trans_',sort(exmets));
Ex=strcat('Ex_',exmets);
dico_EX=cell(numel(Ex),2);
dico_EX(:,1)=Ex;
for i=1:numel(exmets)
    r=find(model.S(ismember(model.mets,exmets(i)),:));
    r=intersect(model.rxns(r),exRxns);
    dico_EX(i,2)=r;
end
composite_model.rxns=[internal_rxns_array(:,1);intercluster_trans_array(:,1);trans;Ex];
for ii=1:numel(fields)
    if size( model.(fields{ii}),1)==size( model.rxns,1) && size( model.(fields{ii}),1)~=size( model.mets,1)
        if iscell(composite_model.(fields{ii}))
            rul=[internal_rxns_array(:,ii+1);intercluster_trans_array(:,ii+1);cell(numel(trans),1);cell(numel(Ex),1)];
            rul(cellfun('isempty',rul))=cellstr('');
            composite_model.(fields{ii})=rul;
            exmissing=model.(fields{ii})(ismember(model.rxns, dico_EX(:,2)));
            composite_model.(fields{ii})(ismember(composite_model.rxns, Ex))=exmissing;
        else
            if strcmp(fields{ii}, 'lb')
                composite_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*-1000; ones(numel(Ex),1)*-1000];
            else
                composite_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*median(internal_rxns_array_n(:,ii+1)); ones(numel(Ex),1)*median(internal_rxns_array_n(:,ii+1))];
            end
        end
    else
        if iscell(composite_model.(fields{ii}))
            composite_model.(fields{ii})=[composite_model.(fields{ii});cell(numel(a),1);cell(numel(a),1)];
        else
            composite_model.(fields{ii})=[composite_model.(fields{ii});ones(numel(a),1)*mean(composite_model.(fields{ii}));ones(numel(a),1)*mean(composite_model.(fields{ii}))];
        end
    end
end
composite_model.rev= sparse(numel(composite_model.rxns),1);
composite_model.rev(composite_model.lb<0)=1; 

composite_model = buildRxnGeneMat(composite_model);
composite_model = creategrRulesField(composite_model);


end


