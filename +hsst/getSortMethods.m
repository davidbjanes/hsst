
function [ output_obj, output_string ] = getSortMethods()
    
    obj = what('hsst/sortMethod');
    output_string = obj(1).classes;
    
    sortMethodObjList = {};
    for n = 1:length(output_string),
        [~,sortMethodObjName,~] = fileparts(output_string{n});
        
        try 
            sortMethodObjList{end+1} = eval(['hsst.sortMethod.' sortMethodObjName '.empty(1,0)']);
        catch error_message
            fprintf('WARNING: Class ''hsst.sortMethod.%s'' is improperly defined. ''hsst.sortMethod.%s'' will be ignored. \n', sortMethodObjName, sortMethodObjName);
        end
    end
    
    output_obj = sortMethodObjList;  %[sortMethodObjList{:}];
            
end

