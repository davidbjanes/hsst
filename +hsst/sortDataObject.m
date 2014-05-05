
function [ dataObj ] = sortDataObject( dataObj, sortMethodObj, sortParameters )                                                    
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
    %% Check dataObj
    if ~strcmp('hsst.dataObject', class(dataObj)),
        error('Invalid hsst.dataObject Class Object');
        dataObj = [];
        return
    end
    
    
    %% Get Values
    wf_snippets = dataObj.wf_snippets;
    time_stamps = dataObj.time_stamps;
    sampling_freq = dataObj.sampling_freq;
    noise_estimate = dataObj.noise_estimate;
    
    
    %% Call hsst.sortSnips
    [ ~, ~, new_dataObj ] = hsst.sortSnips( wf_snippets, ...
                                            time_stamps, ...
                                            sampling_freq, ...
                                            noise_estimate, ...
                                            sortMethodObj, ...
                                            sortParameters );
                         
    %% Update dataObj
    for n = 1:new_dataObj.num_sortCodes,
        dataObj.addSortCode(new_dataObj.sortCodeList{n}, ...
                            new_dataObj.sortMethodNameList{n},...
                            new_dataObj.sortParameterList{n},...
                            new_dataObj.sortScoreObjList{n});
    end
    
    %% Return Result
    dataObj = dataObj;    
    
    
end

