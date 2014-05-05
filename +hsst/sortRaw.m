
function [ sort_IDs, parameter, property ] = sortRaw( rawWaveform, sample_freq, sortMethodObj, sortParameters )                                                    
% sortRaw - sorts raw waveform 
% Written by David Bjanes (dab194@pitt.edu)
%  
% Algorithm
%   1 : Extracts Snippets from raw waveform using the given extraction
%       method class
%   2 : Sorts snippets M times, each time with a different parameter
%   3 : Uses HSST to determing optimal parameter setting and returns output
%
% INPUTS
%   rawWaveform         (1xN) : raw waveform 
%   sample_freq         (1x1) : sampling frequency (Hz)
%   sortMethodObj       (1x1) : subclass obj of superclass "hsst.sorter"
%   sortParameters      (1xM) : parameters for sorting method class
%
% INPUTS: NOT IMPLEMENTED (due to only a single option so far)
%   scoreMethodObj      (1x1) : subclass obj of superclass "hsst.scorer"
%   extractorMethodObj  (1x1) : subclass obj of superclass "hsst.scorer"
%
% Outputs
%   sort_IDs            (1xN) : vector of cluster id integers 
%   parameter           (1x1) : ideal parameter value from input "sortParameters"
%   property            (1x1) : dataObject Class instance containing sort 
%                               information from each parameter listed 
%                               in the input "sortParameters" 
%
    %% extractorMethodObj definition (Could be an input variable,
    % but there is only one option so far)
    extractorMethodList = hsst.getExtractorMethods();
    extractorMethodObj = extractorMethodList{1};
    
    
    %% Check extractorMethodObj is correct class type (if it was an input)
    [ boolean_output, message ] = hsst.isObjectExtractorMethod( extractorMethodObj );
    if ~boolean_output, error(message); end
      
    
    %% Extract Snippets
    threshold = [];
    [snippets, time_stamps, property] = extractorMethodObj.getSnippets( rawWaveform, ...
                                                                        sample_freq, ...
                                                                        threshold );
    noise_estimate  = property.noise_estimate;
    threshold       = property.threshold;
       
    
    %% Call "sortSnips"
    [ optimalSortCode, optimalParam, dataObj ] = hsst.sortSnips( snippets, ...
                                                                 time_stamps, ...
                                                                 sample_freq, ...
                                                                 noise_estimate, ...
                                                                 sortMethodObj, ...
                                                                 sortParameters );
    
    %% Update dataObj
    updateDataObj(dataObj, 'wf_continuous', rawWaveform);
    updateDataObj(dataObj, 'threshold_value', threshold);
    updateDataObj(dataObj, 'alignment_sample_num', property.window_align);

    
    %% Return Result
    sort_IDs  = optimalSortCode;
    parameter = optimalParam;
    property  = dataObj;    
    
    
end

