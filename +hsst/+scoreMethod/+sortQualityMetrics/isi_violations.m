
% sortQualityMetric of SNR values for sortQualityClass
classdef isi_violations < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'ISI Violations';
        verbose_label = 'ISI - Percent of Violators'
    end
    
    methods (Static)
        
        function [score percent_violations] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val
        
            %% VARIABLES FOR METRIC
            MIN_ISI_TIME = 0.001;           % In Seconds   
            MAX_PERCENT_OF_VIOLATIONS = 5;	% In (Int) Percent
            
            
            %% METRIC ANALYSIS CODE
            
            ts_diff = diff(ts);    

            num_isi_violators = sum(ts_diff <= MIN_ISI_TIME);

            percent_violations = num_isi_violators/length(ts_diff) * 100;
            score = percent_violations <= MAX_PERCENT_OF_VIOLATIONS;
                           
            
            %% PLOT Functionality
            if property_struct.showplots,

                interval = MIN_ISI_TIME;
            
                ts_diff = ts_diff(ts_diff < MIN_ISI_TIME*10);
                X = 0:interval:MIN_ISI_TIME*10;
                histogram = histc(ts_diff, X(1:end));
                norm_histogram = histogram ./ sum(histogram);
                
                bar(X+MIN_ISI_TIME/2, norm_histogram);
                hold on
                grid on
                v = axis;
                hmax = v(4);
                plot(MIN_ISI_TIME * ones(1,10), linspace(0,hmax,10), '--', 'Color', 'r', 'LineWidth', 2);

                xlabel('ISI [sec]')
                xlim([0 MIN_ISI_TIME*10]);
                
            end
        end
        
    end
end

