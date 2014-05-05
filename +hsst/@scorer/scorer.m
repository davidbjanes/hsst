
% Abstract scoring class 
classdef scorer < handle
    
    % Public Viewable Properties
    properties (Abstract, Constant)         
        scoreMethodLabel         % I.e. 'SQM'
    end
        
    % Public Settable Properties
    properties (Access = public)
        verbose = false;
    end      
    
    % Abstract Methods
    methods (Abstract)
              
        % INPUTS ----------------------------------------------------------
        % none
        % OUTPUTS ---------------------------------------------------------
        % score     - [0 to 1] value 
        score = getScore(this);       
        
    end
    
    methods (Access = public)
        
        function obj = scorer(sortCode, waveforms, timeStampData, noiseEstimate)
            
            % DO NOTHING with these variables.  Just ensure they are given
%             obj.wf              = waveforms;
%             obj.ts              = timeStampData;
%             obj.noiseEstimate   = noiseEstimate;
%             obj.sort_IDs        = sortCode;
            
        end
        
    end
          
end