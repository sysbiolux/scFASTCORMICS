function[model]=removeRxns_maria(model,rxns)
fields = fieldnames(model);
fields=setdiff(fields,'genes');
all_rxns=getfield(model, 'rxns');
all_mets=getfield(model, 'mets');
if numel(all_mets)== numel(all_rxns)
    disp('you need identify which field is a rxns_realted')
else
    
    remove=ismember (all_rxns, rxns);
    for i=1: numel(fields)
        items=getfield(model, fields{i});
        
        if size(items,1)== numel(all_rxns)
            
            items(remove,:)=[];
            model = rmfield(model,fields(i));
            model=setfield(model,fields{i},items);
            
        end
        if  size(items,2)== numel(all_rxns)
            items(:,remove)=[];
            model = rmfield(model,fields(i));
            model=setfield(model,fields{i},items);
            
        end
        
    end
    Struct_matrix=model.S~=0;
    remove=model.mets(sum(Struct_matrix,2)==0);
    [model]=removeMets_maria(model,remove);
    
end
