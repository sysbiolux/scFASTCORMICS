

    %% %%%%%% LOOPING OVER THE NUMBER OF CLUSTERS/CELLPOPULATIONS & THE FIELDS INSIDE THE MODEL OBJECT
    % this has multiple outputs: 
    % - the S matrix for the internal_rxns and intercluster transport is filled
    % - S matrix: it is assembled then after the loop
    % - intercluster_trans_array, internal_rxns_array entail the slots of the
    % expanded model as columns, they are filled during the loop and then later
    % assembled and put into the right slot of the expanded model
    % - and the genes slot is expanded - just the prefix according to which
    % cluster it belongs - then the rules need to be adjusted accordingly 

    % summary:  3 main task in this for loop
    % - smatrix and slots of the model are expanded for the internal and
    % intercluster rxns (the rest is added after the loop) 
    % - the genes are copied for every cluster + added postfix 
    % - since the gene names change by adding the postfix the rules also need
    % to be adjusted accordingly 


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
        % loops through every cell population 
        % the previous section  created the empty S matrix, this section does
        % loop through all cell population through every slot of the model
        % object to expand the other slots of the model (not only S needs to be
        % expanded

        % genes are copied with a postfix for all clusters, all genes can be
        % copied! no further action needed
        postfix = strcat( '_', num2str(i));
        expanded_input_model.genes=[expanded_input_model.genes;strcat(model.genes,postfix)];

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
        % now the 3 empty matrices created in the first part of the function
        % are filled with values
        internal_S_matrix(ind: ind2,ind3:ind4)=internal_S_matrix_for_1_cluster;
        intercluster_trans_Smat(ind5:ind6,ind7:ind8)=intercluster_trans_Smat_4_1_cluster;
        intercluster_trans_Smat2(ind9:ind10,ind7:ind8)=intercluster_trans_Smat_4_1_cluster2;
        % the internal rxns and mets are copied per cluster as well as the
        % metabolites which are intracellular part of the transport reactions
        internal_rxns_array(ind3:ind4,1)=strcat(internal_rxns,postfix);
        mets(ind: ind2)=strcat(internal_mets,postfix);
        intercluster_trans_array(ind7:ind8,1)=(strcat(intercluster_transport,postfix));

        for ii=1:numel(fields)
            % looping through every slot of the model

            % either we go into the slots defining the characteristics of the
            % metabolites - slot b 
            if size( model.(fields{ii}),1)==size( model.mets,1)&& size( model.(fields{ii}),2)~=size( model.rxns,1)
                temp=model.(fields{ii});

                expanded_input_model.(fields{ii})= [expanded_input_model.(fields{ii});temp(internal_mets_index)];

            % or into the slots defining the characteristics of the rxns -
            % c,lb,ub,rev,rules
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
                                % in case there is only one gene for the rxns,
                                % then only one gene must be changed
                                temp3=temp2(bracket1_pos+1: bracket2_pos-1);
                                temp4=eval(temp3);
                                temp{counter1}=strrep(temp2,temp3, num2str(temp4));
                            else
                                while ~isempty(plus_pos)
                                    % if there are multipel genes for one rxns,
                                    % then there is more then one plus,
                                    % therefore this while loop eradicates all
                                    % the pluses until there are none left 
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
                % some of the slots in the model are strings some are numeric 
                % therefore there are two arrays for the internal and transport
                % rxns - one with the strings and one with postfix n (numeric)
                if iscell(temp)
                    %disp(fields{ii})
                    internal_rxns_array(ind3:ind4,ii+1)=temp(internal_rxns_index);
                    intercluster_trans_array(ind7:ind8,ii+1)=temp(intercluster_transportID);
                else
                    internal_rxns_array_n(ind3:ind4,ii+1)=temp(internal_rxns_index);
                    %disp(fields{ii})
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
    
    
    
%% Finalizing the S-matrix 
% - the expansion is done for the S matrix part belonging to the internal
% and intercluster transport rxns 
% - they therefore need to be assembled and the exchange and transport rxns
% need to be added (they are just identity matrices which can be created
% using the arrays in which the names of these rxns/metabolites are stored)

varname = regexp(exmets, '\[(.*?)\]', 'match', 'once');
expanded_input_model.mets=[mets';sort(strrep(exmets,varname,'[u]'))];
umets=strrep(exmets,varname,'[u]');
S=[internal_S_matrix,intercluster_trans_Smat];

S2=[sparse(size(umets,1),size(internal_S_matrix,2)),intercluster_trans_Smat2];
S=[S;S2];
S7=sparse(size(intercluster_trans_Smat,1), numel(exRxns));
S8=[S7;eye(size(intercluster_trans_Smat2,1))*-1];
% S2=[S;sparse(size(intercluster_trans_Smat2,1),size(S,2))];
S=[S, S8]; 


expanded_input_model.S=S;




%% adding all the rxn names and add to .rxns slot 
trans=strcat('Trans_',sort(exmets));
expanded_input_model.rxns=[internal_rxns_array(:,1);intercluster_trans_array(:,1);trans];


%% doing what we just did for the rxns - for all the other fields in the model object 
% - pasting together the information we get from the arrays created in the
% for internal and intercluster rxns and adding the information for the
% exchange and transport rxns (which we know because they are exchange and
% transport reactions) so it is clear that for example the lower bound is
% set to -1000, and can be done for all exchange and transport rxns

%%% QUESTIONS: 
% - why do we here take the median and the mean for the rev,b and c slot
% for the exchange and transport rxns ? 

for ii=1:numel(fields)
    if size( model.(fields{ii}),1)==size( model.rxns,1) && size( model.(fields{ii}),1)~=size( model.mets,1)
        if iscell(expanded_input_model.(fields{ii}))
            
            rul=[internal_rxns_array(:,ii+1);intercluster_trans_array(:,ii+1);cell(numel(trans),1)];
            rul(cellfun('isempty',rul))=cellstr('');
            expanded_input_model.(fields{ii})=rul;
            
        else
            if strcmp(fields{ii}, 'lb')
                
                %% lower bounds set -1000 by defaut for trans and Ex
                expanded_input_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*-1000];
            else
                % c & rev  - why median ? 
                expanded_input_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*median(internal_rxns_array_n(:,ii+1))];
            end
        end
    else
        if iscell(expanded_input_model.(fields{ii}))
            [~,index]=ismember(exmets,model.mets);
            temp2=cell(numel(exmets),1);
            temp=model.(fields{ii});
            temp2(index>0)=temp(index(index>0));
            expanded_input_model.(fields{ii})=[expanded_input_model.(fields{ii});temp2];
        else
            % b slot - why the mean ? 
            expanded_input_model.(fields{ii})=[expanded_input_model.(fields{ii});ones(numel(umets),1)*mean(expanded_input_model.(fields{ii}));ones(numel(umets),1)*mean(expanded_input_model.(fields{ii}))];
        end
    end
end











