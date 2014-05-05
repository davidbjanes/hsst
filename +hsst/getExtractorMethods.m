
function [ output_obj, output_string ] = getExtractorMethods()
    
    obj = what('hsst/extractorMethod');
    output_string = obj(1).classes;
    
    extractorMethodObjList = {};
    for n = 1:length(output_string),
        [~,extractorMethodObjName,~] = fileparts(output_string{n});
        
        try 
            extractorMethodObjList{end+1} = eval(['hsst.extractorMethod.' extractorMethodObjName '.empty(1,0)']);
        catch error_message
            fprintf('WARNING: Class ''hsst.extractorMethod.%s'' is improperly defined. ''hsst.extractorMethod.%s'' will be ignored. \n', extractorMethodObjName, extractorMethodObjName);
            fprintf('ERROR MESSAGE: %s', error_message.message);
        end
    end
    
    output_obj = extractorMethodObjList;  %[extractorMethodObjList{:}];
            
end

