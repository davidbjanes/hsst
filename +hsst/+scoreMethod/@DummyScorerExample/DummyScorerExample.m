
classdef DummyScorerExample < hsst.scorer
    
    % Public Viewable Properties
    properties (Constant)         
        scoreMethodLabel = 'Dummy Example Scorer';
    end
    
    % Visible (but nonEditable) Properties
    properties (SetAccess = private)
        score    
    end
    
    methods
              
        % INPUTS ----------------------------------------------------------
        % none
        % OUTPUTS ---------------------------------------------------------
        % score     - [0 to 1] value 
        function score = getScore(this)
            score = this.score;
        end
        
    end
    
    methods
        
        function obj = DummyScorerExample(sortCode, waveforms, timeStampData, ...
                                          noiseEstimate, varargin)
            
        	obj@hsst.scorer(sortCode, waveforms, timeStampData, noiseEstimate)
            
            %% Add scoring method(s) here
            
            % Do Something Here
            % 
            %
            %
            %
            
            
            %% Output
            
            obj.score = 1;
            
        end
        
    end
    
end

