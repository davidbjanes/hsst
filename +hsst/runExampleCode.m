
function [ property ] = runExampleCode()


    %% GENERATE Fake inputs 
    rawWaveform     = randn([1, 1000000]);               % Raw waveform
    spikeWaveform   = randn([32, 1000]);                 % Waveform Snippets
    spikeTimeStamps = unifrnd([0:999],[1:1000]);         % time stamps of snippets
    sampleFrequency = 1000;                              % sampling frequency (Hz)
    noiseEstimate   = mean(std(spikeWaveform,[],2));     % noise estimate
    
    
    %% Load Example Neural Data
    load('NeuralExampleDataset.mat');    
    
    
    %% Get Methods of Sorting/Scoring
    sortMethods    = hsst.getSortMethods();
    sortMethodObj  = sortMethods{6};
    
    % OR ALTERATIVELY
    sortMethodObj  = hsst.sortMethod.MATLAB_GMM;
    
    % select range of input parameters (unique values for each sorting
    % method)
    sortParameters = [1:6];  
    
    
    %% Run HSST on fake snippets
    [ sort_IDs, parameter, property ] = hsst.sortSnips(spikeWaveform, ...
                                                       spikeTimeStamps, ...
                                                       sampleFrequency, ...
                                                       noiseEstimate, ...
                                                       sortMethodObj, ...
                                                       sortParameters);
    fprintf('The optimal parameter for sorting with ''%s'', is: %d. \n', class(sortMethodObj), parameter);   
    

    %% Run HSST on raw waveforms
    [ sort_IDs, parameter, property ] = hsst.sortRaw(rawWaveform, ...
                                                     sampleFrequency, ...
                                                     sortMethodObj, ...
                                                     sortParameters);
    fprintf('The optimal parameter for sorting with ''%s'', is: %d. \n', class(sortMethodObj), parameter);   
    
    
    %% Get Unit Information
    param_ind = cellfun(@isequal, repmat({parameter}, [1 property.num_sortCodes]), property.sortParameterList);
    property.sortScoreObjList{param_ind}.noise_units
    property.sortScoreObjList{param_ind}.undersorted_units                       
    property.sortScoreObjList{param_ind}.oversorted_units
        
    hsst.gui(property);
    
end

