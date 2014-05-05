
% extraction method class for HSST
classdef thresholding < hsst.extractor
    
    % Public Viewable Properties
    properties (Constant)         
        extractorMethodLabel = 'Thresholding';
    end 
    
    % Abstract Methods
    methods (Static)
        
        [snippets, time_stamps, property] = getSnippets( rawWaveform, sample_freq, threshold );
        
        [noise,info] = getNoiseEstimate(data,fs,varargin)
        
    end
       
end

