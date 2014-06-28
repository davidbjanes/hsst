
% extraction method class for HSST
classdef DummyExtactorMethod < hsst.extractor
    
    % Public Viewable Properties
    properties (Constant)         
        extractorMethodLabel = 'Thresholding';
    end 
    
    % Abstract Methods
    methods (Static)
        
        function [snippets, time_stamps, property] = getSnippets( rawWaveform, sample_freq, threshold )
        
            %% Spike Detection Algorithm Here
            
            %
            %
            %
            %
            
            %% Outputs
            snippets = 1;                   % Waveform Snippets
            time_stamps = 1;                % Time Stamps
            
            property = struct;
            property.noise_estimate = 1;    % Estimate of noise floor
            property.threshold = 1;         % Spike Detection Threshold
    
        end
        
    end
       
end

