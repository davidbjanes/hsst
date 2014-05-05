
% sortQualityMetric of SNR values for sortQualityClass
classdef isi_exp_fit < hsst.scoreMethod.sortQualityMetrics.metric 
    
    properties (Constant)
        label = 'ISI - Exp Fit';
        verbose_label = 'ISI - MSE of Exp Fit'
    end
    
    methods (Static)
        
        function [score mserror] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val
        
            %% VARIABLES FOR METRIC
            INTERVAL = 0.001;   % In Seconds
            MAX_ISI = 0.1;   % seconds
            MIN_MSE_PERCENT = 10;   % (int) percent
            MIN_ISI_TIME = 0.001;  % seconds
            
            
            %% METRIC ANALYSIS CODE
                                    
            num_wf = length(ts);
            num_wf = fix(num_wf/10^floor(log10(num_wf)))*10^floor(log10(num_wf));
            INTERVAL = 1/num_wf;
            
            ts_diff = diff(ts);
            ts_diff = ts_diff(ts_diff < MAX_ISI);
            
            % Histogram of interspike intervals of a 10us bins
            X = 0:INTERVAL:MAX_ISI;
            h = histc(ts_diff, X);

            % MSE of exponential fit to ISI histogram
            if length(ts_diff) > 10,
                
                % Normalize fitted curve by area of histogram
                h_area = sum(h.*INTERVAL);

                probdist_exp = fitdist(ts_diff', 'exponential');
                fitted_dist_exp = pdf(probdist_exp, 0:INTERVAL:MAX_ISI);
                fitted_dist_exp = fitted_dist_exp .* h_area;

                % Normalized root-mean-square-deviation 
                % http://en.wikipedia.org/wiki/Root-mean-square_deviation
                MSE_expFit_error = sqrt(mse(h - fitted_dist_exp))/(max(h) - min(h));
            
            else
                MSE_expFit_error = 1;
            end
           
            mserror = MSE_expFit_error*100;
            score = mserror >= MIN_MSE_PERCENT;
            
            
            %% PLOT Functionality
            if property_struct.showplots,

                bar(X, h);
                hold on
                grid on
                hmax = max(max(h), 2);
                plot(MIN_ISI_TIME * ones(1,hmax+1), 0:hmax, '--', 'Color', 'r', 'LineWidth', 2);

                axis tight
                v = axis;
                str2(1) = {sprintf('MSE - Exp. Fit = %2.1f%%', mserror)};
                text(v(2)*.95, v(4)*.925, str2, 'HorizontalAlignment', 'right');

                % Print fitted distribution  
                if length(ts_diff) > 10, 
                    plot(X, fitted_dist_exp, 'Color', 'c', 'LineWidth', 2)
                end

                xlabel('ISI [sec]')
                xlim([0 MAX_ISI]);
                
            end
            
        end
                
    end
end
