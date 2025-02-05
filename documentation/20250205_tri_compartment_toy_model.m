% this code is from the build_expanded_input_model function & compute
% expanded input model


%%
%%%%%%%%%%%%%% CLASSIFYING THE RXNS AND METABOLITES 
% - ARE THE RXNS INTRA CELLULAR/TRANSPORT OR EXCHANGE RXNS ?
% - ARE THE METS A PLAYER IN AN INTRA/TRANSPORT OR EXCHANG RXN? IN WHICH
% COMPARTMENT ARE THE METABOLITES -> [U]/ [E] OR INTRACELLULAR ? 
% mets: extracellular metabolites, u compartment metabolites, intracellular
% metabolites 
% rxns: intercluster_rxns(intracellular), intercluster_transport, transport
% into compartment -> trans, exchang reactions -> transport into the model


%% get metabolites in [e] and the corresponding exchange Rxns

model.genes=regexprep(model.genes,'\.[0-9]+$','');
% get external rxns - searching the S matrix for columns with only one
% entry - 1 entry 
exRxns=find(sum(abs(model.S),1)==1);
% find the corrsponding metabolite -> through S matrix
[exmets,~]=find(model.S(:, exRxns));
% from all the metabolites which are participating in the exchange rxns,
% only the ones should be choosen which are in the e compartment, not the
% metabolites within the cell -> therefore filtering for [e]
exmets=exmets(contains(model.mets(exmets),'[e]'));
exmets=unique(model.mets(exmets));
% find the index within the model.mets
[~,exmetsID,~]=intersect(model.mets, exmets);
% find again the corresponding rxns
[~,r]=find(model.S(exmetsID,:));
% make sure that the found rxns are ex rxns found in the beginning
exRxns=intersect(exRxns,r);
% find the names from the indices
exRxns=model.rxns(exRxns);

u_mets=unique(regexprep(exmets,'\[e]','[u]'));



% find all the metabolites which are not in the e mets
% there are two types of rxns which are not exchange, intracellular
% (internal) and intercluster transport rxns - internal_2_u_mets are
% therefore the metabolites intracellular which take part in a transport
% from intracellular to [u] compartment
internal_mets=setdiff(model.mets,exmets);
% get the indices
[~,internal_metsID,~]=intersect(model.mets, internal_mets);
% find the rxns where the internal mets are taking part in
[~,internal_rxnsID]=find(model.S(internal_metsID, :));
% get the names of the rxns where internal mets take part in
internal_rxns=model.rxns(internal_rxnsID);
internal_rxns=unique(internal_rxns);
% find the ids of the exchange reactions -> this will be the new
% intercluster transport rxns -> transporting from the cells to the u
% compartment 
[~,intercluster_transport_ID]=find(model.S(exmetsID, :));
% find the names of the intercluster transport rxns
intercluster_transport=model.rxns(intercluster_transport_ID);
% get rid of the ex reactions in the internal_rxns - internal rxns was
% determined by getting all the rxns the internal mets take part in, but
% internal mets also take part in ex rxns, therefore we need to get rid of
% those
internal_rxns=setdiff(internal_rxns,intercluster_transport);


% find inter(intra?)cluster_transport(internal_rxns), and metabolites that are transported between the clusters and the intersticial  (internal_2_u_mets),
intercluster_transport=setdiff(intercluster_transport,exRxns); % these are the ones which are exRxns but are not associated with a [e] metabolite
[~,intercluster_transport_ID,~]=intersect(model.rxns, intercluster_transport);
[internal_2_u_mets_ID,~]=find(model.S(:, intercluster_transport_ID));
internal_2_u_mets=model.mets(internal_2_u_mets_ID);
internal_2_u_mets=intersect(internal_2_u_mets,internal_mets);



%% create a empty expanded_input_model
[expanded_input_model, fields]=create_empty_expanded_input_model(model, fields_to_keep);

% exclude the fields you will not be looping over in the for loops, for
% which the input will be set manually, not via for loop
fields=setdiff(fields,'genes');
fields=setdiff(fields,'rxnGeneMat');
fields=setdiff(fields,'grRules');
fields=setdiff(fields,'mets');
fields=setdiff(fields,'rxns');
fields=setdiff(fields,'S');


% 
mets_connected_to_ext_metabolite = internal_2_u_mets;
umets = u_mets; 
%

%% CONSTRUCTION OF SMATRIX FROM SMALLER MATRICES 

% initializes an EMPTY SPARSE MATRICES that has the correct size

% mets/rxns which are intracellular/ not connected to u compartment
internal_S_matrix=sparse(numel(internal_mets)*number_of_cluster,numel(internal_rxns)*number_of_cluster);
% rxns which transport to/from u compartment to intracellular x metabolites
% which are the intracellular part of these transport rxns (therefore need
% to be copied for every cluster/population)
intercluster_trans_Smat=sparse(numel(internal_mets)*number_of_cluster,numel(intercluster_transport)*number_of_cluster);
% rxns which transport to/from u compartment to intracellular (same as
% above) but x metabolites which are laying in the [u] compartment of the
% transport rxns - therefore 1 u compartment means these mets do not have to
% be copied
intercluster_trans_Smat2=sparse(numel(umets),numel(intercluster_transport)*number_of_cluster);

% finds the indices of the internal metabolites and reactions
[~,internal_mets_index, ~]=intersect(model.mets,internal_mets);
[~,internal_rxns_index, ~]=intersect(model.rxns,internal_rxns);
[~,intercluster_transportID, ~]=intersect(model.rxns,intercluster_transport);


% builds a matrix with the internal reactions and metabolites of a sub-model
internal_S_matrix_for_1_cluster=model.S(internal_mets_index,internal_rxns_index);

% initializes an EMPTY SPARSE ARRAYS FOR THE RXNS 
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
    % now the 3 empty matrices created in the first part of the function
    % are filled with values
    internal_S_matrix(ind: ind2,ind3:ind4)=internal_S_matrix_for_1_cluster;
    intercluster_trans_Smat(ind5:ind6,ind7:ind8)=intercluster_trans_Smat_4_1_cluster;
    intercluster_trans_Smat2(ind9:ind10,ind7:ind8)=intercluster_trans_Smat_4_1_cluster2;
    % the internal rxns and mets are copied per cluster as well as the
    % metabolites which are intracellular part of the transport reactions
    internal_rxns_array(ind3:ind4,1)=strcat(internal_rxns,'_', num2str(i));
    mets(ind: ind2)=strcat(internal_mets,'_', num2str(i));
    intercluster_trans_array(ind7:ind8,1)=(strcat(intercluster_transport,'_', num2str(i)));
   
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

%% adding all the rxn names and add to .rxns slot 

% creating a Trans and Ex rxn for every ex met in u compartment
trans=strcat('Trans_',sort(exmets));
Ex=strcat('Ex_',exmets);
% dico_EX=cell(numel(Ex),2);
% dico_EX(:,1)=Ex;
% % find the rxns for the ex mets 
% for i=1:numel(exmets)
%     r=find(model.S(ismember(model.mets,exmets(i)),:));
%     r=intersect(model.rxns(r),exRxns);
%     dico_EX(i,2)=r;
% end
% add all the rxns to the rxns slot of the model
expanded_input_model.rxns=[internal_rxns_array(:,1);intercluster_trans_array(:,1);trans;Ex];

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
            
            rul=[internal_rxns_array(:,ii+1);intercluster_trans_array(:,ii+1);cell(numel(trans),1);cell(numel(Ex),1)];
            rul(cellfun('isempty',rul))=cellstr('');
            expanded_input_model.(fields{ii})=rul;
            
        else
            if strcmp(fields{ii}, 'lb')
                
                %% lower bounds set -1000 by defaut for trans and Ex
                expanded_input_model.(fields{ii})=[internal_rxns_array_n(:,ii+1);intercluster_trans_array_n(:,ii+1);ones(numel(trans),1)*-1000; ones(numel(Ex),1)*-1000];
            else
                % c & rev  - why median ? 
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
            % b slot - why the mean ? 
            expanded_input_model.(fields{ii})=[expanded_input_model.(fields{ii});ones(numel(umets),1)*mean(expanded_input_model.(fields{ii}));ones(numel(umets),1)*mean(expanded_input_model.(fields{ii}))];
        end
    end
end




%% Correct reversibility of the rxns

lb=zeros(numel(exmets),1);
ub=zeros(numel(exmets),1);
grRules=cell(numel(exmets),1);
R_ex=cell(numel(exmets),1);

% since the ex rxns are still defined as 1 in the s matrix the sign needs
% to be flipped - toy model convention ? Recon convention ???
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
% this part I do not get... - 
% so we create lb and ub with all zeros, then we manipulate this 
% we get the indices of the R_ex in the model.rxns 
[~,index]=ismember(R_ex,model.rxns);
% then we filter out all the R_ex which are not in the model.rxns by
% excluding the zero entrie in the indices vector
lb(index>0)=model.lb(index(index>0));
ub(index>0)=model.ub(index(index>0));
% transfer the lb & ub from the consensus model to the expanded model, the
% indices of the ex rxns change between those two models, therefore this
% transfer is needed

Table_bound=table(R_ex,lb, ub, exmets,strcat('Ex_',exmets));
[~,index]=ismember(expanded_input_model.rxns,table2array(Table_bound(:,5)));
expanded_input_model.lb(index>0)=Table_bound.lb(index(index>0));
expanded_input_model.ub(index>0)=Table_bound.ub(index(index>0));



%%
% now the reversibility slot is fixed, for all the rxns the lb were set,
% so for the rxns which have a negative lb we can assume that they are
% reversible 
expanded_input_model.rev= sparse(numel(expanded_input_model.rxns),1);
expanded_input_model.rev(expanded_input_model.lb<0)=1;


%%
% building the RxnGene matrix from grRules slot
expanded_input_model = buildRxnGeneMat(expanded_input_model);
% creating the .rule slot from grRules
expanded_input_model = creategrRulesField(expanded_input_model);
% fixing Irr rxns ? what does that mean? 
expanded_input_model=fixIrr_rFASTCORMICS(expanded_input_model);

% % this is the fixIrr function code: 
% the code check if the rxn is really reversible, looking at the lower
% upper bound, if the upper bound is lower then 0 or the lower bound
% bigger then 0  then the rxn is not reversible rev = 0 
% this works for the rxns which have a lower bound bigger than 0, since
% the rxns is defined A + B <=> C but if the ub is <= 0 then the rxns
% needs to be turned around 
% model.rev = zeros(numel(model.rxns),1);
% % the  
% model.rev(model.lb <0 & model.ub> 0) = 1;
% Irr=(model.lb >=0 & model.ub>0| model.ub<=0 & model.lb<0);
% model.rev(Irr) = 0;
% FakeIrr= model.ub <=0 & model.lb<0;
% model.S(:, FakeIrr) = -model.S(:,FakeIrr);
% model.ub(FakeIrr) = -model.lb(FakeIrr);
% model.lb(FakeIrr) = zeros(sum(FakeIrr),1);

%% 
% correct rules for exchanges
% question: the exchange rxns are empty bevore and after -> 
% first of all, why do this code in the frist place and why should exhange
% reactions have a gprule ? 

[~,index]=ismember(R_ex,model.rxns);
grRules(index>0)=model.grRules(index(index>0));
% getting the rules for the exchange rxns from the consistent generic model
% then adding to the table for the exchange reactions 
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

% update the rules field after altering the grRules slot
expanded_input_model = generateRules(expanded_input_model);






