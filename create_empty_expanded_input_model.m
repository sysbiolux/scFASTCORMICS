function[composite_model, fields]=create_empty_expanded_input_model(model, fields_to_keep)
%(c)Maria Pires Pacheco
composite_model=struct();
% retrieve the fields from the input model
fields = fieldnames(model);%retrieve the fields of the model;

if ~isempty(fields_to_keep)
    fields= intersect(fields,fields_to_keep);
end
for i=1:numel(fields)
    if iscell(model.(fields{i}))
        composite_model.(fields{i})=cell(0,1);
    else
        composite_model.(fields{i})=zeros(0,1);
    end    
end