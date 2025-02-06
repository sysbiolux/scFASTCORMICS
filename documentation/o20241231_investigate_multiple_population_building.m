

%% 
% this code creates compartments with 3 cellclusters and the second time it 
% is expanded it creates different compartments 

% expand the single celltype model
number_of_cluster = 2;
o20250206_expand_S_matrix
transport_intercompartment = intercluster_transport;
model = expanded_input_model;
% expand the single compartment model
number_of_cluster = 3;
o20250206_expand_S_matrix

if small_model 
    graphObj=createMetIntrcNetwork(expanded_input_model,expanded_input_model.mets);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor doesend
end
%% 
% 

model = expanded_input_model
%%

%%
if small_model 
    graphObj=createMetIntrcNetwork(expanded_input_model,expanded_input_model.mets);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does
end
%%

% now I need to add the exchange and transport rxns & metabolites
%% 
% * first the exchang rxns 
% * from which of the compartments do I add the exchange rxns ? -> 1 & 3

% find transport reactions 


%% e 1


transport = find((contains(model.rxns,"Trans_")));

transport_to_ex = find((contains(model.rxns,"Trans_")& (contains(model.rxns,"_1"))));
e_mets=unique(regexprep(model.rxns(transport_to_ex),'Trans_',''));

num_ex_mets = length(transport)/number_of_cluster;

new_S = sparse(num_ex_mets, size(model.S,2));
new_S(:,transport_to_ex) = [eye(num_ex_mets)];

%%


model.S = [model.S;new_S];
model.mets = [model.mets;e_mets];%%


%% e 3


transport = find((contains(model.rxns,"Trans_")));

transport_to_ex = find((contains(model.rxns,"Trans_")& (contains(model.rxns,"_3"))));
e_mets=unique(regexprep(model.rxns(transport_to_ex),'Trans_',''));

num_ex_mets = length(transport)/number_of_cluster;

new_S = sparse(num_ex_mets, size(model.S,2));
new_S(:,transport_to_ex) = eye(num_ex_mets);

model.S = [model.S;new_S];
model.mets = [model.mets;e_mets];%%

model_with_2_e = model;
%%
if small_model 
    graphObj=createMetIntrcNetwork(model_with_2_e,model_with_2_e.mets);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does
end


% next the exchange between the compartments & ex: 

% this will not change the model in the figure, sine it only adds exchange
% reactions to supply with e metabolites
number_of_ex_envs = 2;
model.S = [model.S,[sparse(size(model.S,1)-( length(transport_to_ex)* number_of_ex_envs),length(transport_to_ex)* number_of_ex_envs);eye(length(transport_to_ex) * number_of_ex_envs)]];
model.rxns = [model.rxns;strcat(model.rxns(transport_to_ex),"_to_1");strcat(model.rxns(transport_to_ex),"_to_2")];
% next the exchange between compartments
% 


ele_transport_matrix = [sparse(((size(model.S,1) - number_of_ex_envs * length(transport_to_ex))/3) - length(transport_intercompartment),length(transport_intercompartment)); ...
                        eye(length(transport_intercompartment))];
intercompartment_rxns = [[[ele_transport_matrix;ele_transport_matrix *-1;sparse(size(ele_transport_matrix,1),size(ele_transport_matrix,2))],[sparse(size(ele_transport_matrix,1),size(ele_transport_matrix,2));ele_transport_matrix;ele_transport_matrix*-1]];sparse(length(transport_to_ex) * number_of_ex_envs,length(transport_to_ex) * number_of_ex_envs)];

full_model = model;
full_model.S = [model.S, intercompartment_rxns];
full_model.rxns = [model.rxns; [strcat(transport_intercompartment, "_1");strcat(transport_intercompartment, "_2")]];



%%
size(model.S)
%%
if small_model 
    graphObj=createMetIntrcNetwork(full_model,full_model.mets);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does
end

%% 
% 


