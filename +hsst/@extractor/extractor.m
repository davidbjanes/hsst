
% Abstract extraction class for HSST
classdef extractor < handle
    
    % Public Viewable Properties
    properties (Abstract, Constant)         
        extractorMethodLabel         % I.e. 'Thresholding'
    end
    
    % Public Settable Properties
    properties (Access = public)
        verbose = false;
    end  
    
    % Abstract Methods
    methods (Abstract, Static)
        
        [snippets, time_stamps, property] = getSnippets( rawWaveform, sample_freq, threshold );
        
    end
    
    % Public Methods
    methods
        
        function obj = extractor(verbose)
            
            if nargin > 1,
                obj.verbose = verbose;
            end
            
            if obj.verbose,
                fprintf('Sort Object Created: %s \n', obj.extractorMethodLabel)
            end

        end
                
    end

        
end

