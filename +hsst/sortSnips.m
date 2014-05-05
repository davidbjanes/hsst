
function [ sort_IDs, parameter, property ] = sortSnips( snippets, time_stamps, sample_freq, noise_estimate, sortMethodObj, sortParameters ) 
% sortSnips - sorts snippet waveforms already extracted
% Written by David Bjanes (dab194@pitt.edu)
%  
% Algorithm
%   1 : Sorts snippets M times, each time with a different parameter
%   2 : Uses HSST to determing optimal parameter setting and returns output
%
% INPUTS
%   snippets        (NxD) : waveform snippets
%   time_stamps     (1xN) : vector of time stamps 
%   sample_freq     (1x1) : sampling frequency (Hz)
%   noise_estimate  (1x1) : estimate of noise floor (same units as snippets)
%   sortMethodObj   (1x1) : subclass obj of superclass "hsst.sorter"
%   sortParameters	(1xM) : parameters for sorting method class
% NOT IMPLEMENTED YET (due to only a single option so far)
%   scoreMethodObj	(1x1) : subclass obj of superclass "hsst.scorer"
%
% Outputs
%   sort_IDs        (1xN) : vector of cluster id integers 
%   parameter       (1x1) : ideal parameter value from input "sortParameters"
%   property        (1x1) : dataObject Class instance containing sort 
%                           information from each parameter listed 
%                           in the input "sortParameters" 
%
    %% scoreMethodObj definition (Could be an input variable,
    % but there is only one option so far)
    scoreMethodList = hsst.getScoreMethods();
    scoreMethodObj = scoreMethodList{2};

    
    %% Check SortMethodObj is correct class type
    [ boolean_output, message ] = hsst.isObjectSortMethod( sortMethodObj );
    if ~boolean_output, error(message); end
    
      
    %% Check ScoreMethodObj is correct class type (if it was an input)
    [ boolean_output, message ] = hsst.isObjectScoreMethod( scoreMethodObj );
    if ~boolean_output, error(message); end
    
    
    %% Build DataObject
    dataObj = hsst.dataObject([], snippets, time_stamps, sample_freq, [], noise_estimate);
            
    
    %% Run Sort Algorithm and Optimizer for Parameter Selection
    optimizer = hsst.optimizer( dataObj, ...
                                sortMethodObj, ...
                                scoreMethodObj, ...
                                sortParameters);
    optimalSortCode = optimizer.returnOptimalSortCode();
    optimalParam    = optimizer.returnOptimalParameter();
    
    
    %% Return Result
    sort_IDs  = optimalSortCode;
    parameter = optimalParam;
    property  = dataObj;
    
    
end

