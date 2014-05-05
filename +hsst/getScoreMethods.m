
function [ output_obj, output_string ] = getScoreMethods()
    
    obj = what('hsst/scoreMethod');
    output_string = obj(1).classes;
    
    scoreMethodObjList = {};
    for n = 1:length(output_string),
        [~,scoreMethodObjName,~] = fileparts(output_string{n});
        
        try 
            scoreMethodObjList{end+1} = eval(['hsst.scoreMethod.' scoreMethodObjName '.empty(1,0)']);
        catch error_message
            fprintf('WARNING: Class ''hsst.sortMethod.%s'' is improperly defined. ''hsst.sortMethod.%s'' will be ignored. \n', scoreMethodObjName, scoreMethodObjName);
        end
    end
    
    output_obj = scoreMethodObjList; %[scoreMethodObjList{:}];
            
end
