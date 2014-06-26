
function [ property ] = runExampleCode()


    %% GENERATE Fake inputs 
    rawWaveform = randn([1, 1000000]);
    wf = randn([32, 1000]);
    ts = unifrnd([0:999],[1:1000]);
    fs = 1000;
    noiseEstimate = mean(std(wf,[],2));
    
    
    %% Get Methods of Sorting/Scoring
    sortMethods    = hsst.getSortMethods();
    sortMethodObj  = sortMethods{6};
    
    sortParameters = [1:3];  
    
    
    %% Run HSST on fake snippets
    [ sort_IDs, parameter, property ] = hsst.sortSnips(wf, ts, fs, noiseEstimate, ...
                                                       sortMethodObj, ...
                                                       sortParameters);
    fprintf('The optimal parameter for sorting with ''%s'', is: %d. \n', class(sortMethodObj), parameter);   
    

    %% Run HSST on raw waveforms
    [ sort_IDs, parameter, property ] = hsst.sortRaw(rawWaveform, fs, ...
                                                     sortMethodObj, ...
                                                     sortParameters);
    fprintf('The optimal parameter for sorting with ''%s'', is: %d. \n', class(sortMethodObj), parameter);   
    
    
    %% Get Unit Information
    param_ind = cellfun(@isequal, repmat({parameter}, [1 property.num_sortCodes]), property.sortParameterList);
    property.sortScoreObjList{param_ind}.noise_units
    property.sortScoreObjList{param_ind}.undersorted_units                       
    property.sortScoreObjList{param_ind}.oversorted_units
        
end

