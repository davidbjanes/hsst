
% sortQualityMetric called template_match for sortQualityClass
classdef threshold_slope < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'Threshold Slope';
        verbose_label = 'Threshold Slope Bimodality'
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
            PLOT_FLAG = property_struct.showplots;
            [~, ~, ALIGNMENT_POINT] = findThresholdCrossing(property_struct.wf);
            ALIGNMENT_POINT = ALIGNMENT_POINT + 1;

            
            %% METRIC ANALYSIS CODE
                    
            num_samples = size(wf,2); % Waveform Sample Length
            num_waveforms = size(wf,1); % Number of waveforms 

            if num_samples >= ALIGNMENT_POINT,
                crossing = wf(:,ALIGNMENT_POINT+[-1 0]);
                threshold_slopes = diff(crossing, [], 2);    

                bimodal = bimodalityDetector(threshold_slopes, PLOT_FLAG);
            else
                bimodal = false;
            end
            
            score = ~bimodal;
            raw_val = 0;            
        
        end
    end
end