function[model]=removeMets_maria(model,remove) 
remove=find(ismember(model.mets, remove));
fields = fieldnames(model);
 all_rxns=getfield(model, 'rxns');
  all_mets=getfield(model, 'mets');
  if numel(all_mets)== numel(all_rxns)
      disp('you need identify which field is a rxns_realted')
  else
      for i=1: numel(fields)
           items=getfield(model, fields{i});
           
           if size(items,1)== numel(all_mets)
               
               items(remove,:)=[];
               model = rmfield(model,fields(i));
               model=setfield(model,fields{i},items);
               
           end
           if  size(items,2)== numel(all_mets)
                items(:,remove)=[];
               model = rmfield(model,fields(i));
               model=setfield(model,fields{i},items);
               
      end
      
  end

end
