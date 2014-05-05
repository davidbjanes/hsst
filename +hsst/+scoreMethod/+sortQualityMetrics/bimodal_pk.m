
% sortQualityMetric of SNR values for sortQualityClass
classdef bimodal_pk < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'Bimodal Pk';
        verbose_label = 'Bimodal Peak'
    end
    
    methods (Static)

        function [score raw_val] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val
        
            import hsst.scoreMethod.sortQualityMetrics.*

            %% VARIABLES FOR METRIC
            NOISE_FLOOR = noiseEst;   
            SIGNIFICANT_AREA = 0.33;
            SIGNIFICANT_NORM_HEIGHT = 0.33;
                       
            
            %% METRIC ANALYSIS CODE           
            
            % Distribution of max amplitudes of the snips
            max_amplitude = max(abs(wf),[],2);
            
            % Conversion into SNR values
            max_amplitude = abs(max_amplitude ./ NOISE_FLOOR);
            
            bimodal = bimodalityDetector(max_amplitude, property_struct.showplots);
                
            score = ~bimodal;
            raw_val = [];
            
            
            %% PLOT Functionality
            if property_struct.showplots,

                % a/a
                
            end
            
                       
        end
        
    end
end

