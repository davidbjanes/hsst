
% sortQualityMetric of cross correlation for sortQualityClass
classdef cross_correlation < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'Cross Correlation';
        verbose_label = 'Cross Correlation Values: ';
    end
    
    methods (Static)

        function [score, num_peaks] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val

            import hsst.scoreMethod.sortQualityMetrics.*

            %% METRIC ANALYSIS CODE
            
            warning ('off','all');

            % Number of waveforms
            min_sample_point_not_met = size(wf,1) < 50;
            if min_sample_point_not_met,
                num_peaks = 1;
                score = 1;
            
            else

                mean_wf = mean(wf);
                
                number_waveforms = size(wf,1);
                                
                xcorr_score = max( xcorr2( mean_wf, wf(:,:) ), [], 2 );
                
                bimodal = bimodalityDetector(xcorr_score, property_struct.showplots);
                
                score = ~bimodal;
                num_peaks = [];
                
                
                %% Attempt 2
%                 % Range of correlation values within 5sigma
%                 std_normalized = abs((xcorr_score-mean(xcorr_score))./std(xcorr_score));
%                 xcorr_score( std_normalized > 5) = [];

%                 % Build histogram
% %                 X = linspace(min(sample_points), max(sample_points), 100);
%                 X = linspace(min(xcorr_score), max(xcorr_score), sqrt(length(xcorr_score)));
%                 
%                 histogram = histc(xcorr_score, X);
%                 histogram = histogram./sum(histogram);
%                 
%                 invaluable_points = histogram > MIN_PERCENT_OF_HISTOGRAM_PER_BIN;
%                 start_ind = find(invaluable_points, 1, 'first');
%                 end_ind = find(invaluable_points, 1, 'last');
%                 invaluable_points(start_ind:end_ind) = 1;
%                 histogram(~invaluable_points) = [];
% 
%                 % Smooth Histogram
%                 N = length(histogram);
%             
%                 if N > 3,
%                     smoothedSlice = histogram;
% 
%                     % 3-Sample-Median Smoothing
%                     for i = 2:N-1,
%                         smoothedSlice(i) = median(histogram(i-1:i+1));
%                     end
%                     smoothedSlice([1 N]) = histogram([1 N]);
%                     histogram = smoothedSlice;            
% 
% 
%                     % Sliding Window Smoothing
%                     smoothedSlice = smooth(histogram,3);
%                     histogram = smoothedSlice';
%                 end
%                 
%                 
%                 
%                 plot_X = X(invaluable_points);
%                 X = 1:length(histogram);
% %                 X = plot_X;
%                     
%                 warning('off', 'all')                
%                 poly_val_n = [];
%                 squ_error = [];
%                 coeff = [6];
%                 poly_fit = polyfit(X,histogram,coeff);
%                 fitted_curve = polyval(poly_fit,X);
%                 warning('on', 'all')
%                     
%                 N = length(histogram);
%                 num_parameters = length(coeff);
%                 
%                 % Find all local maximum
%                 [pks, locs] = findpeaks(fitted_curve);
%                 
%                 % Find all local maximum
%                 [max_pks, max_locs] = findpeaks(fitted_curve);
%                 [min_pks, min_locs] = findpeaks(-fitted_curve);
%                 
%                 all_locs = [1 min_locs max_locs length(fitted_curve)];
%                 [all_locs, ~] = sort(all_locs);
%                 all_pks = [fitted_curve(all_locs)];
%                 
%                 significant_pk_height = 0;
%                 if length(max_pks) > 1,
%                     pk_height = [];
%                     for pk_ind = max_locs,
%                         all_pks_ind = find(all_locs == pk_ind);
%                         valley_min_pk_ind = all_pks_ind + [-1 1];
%                         
%                         pk_height_per_side = all_pks(all_pks_ind) - all_pks(valley_min_pk_ind);
%                         pk_height(pk_ind == max_locs) = max(pk_height_per_side);                        
%                     end
%                     norm_pk_height = (pk_height - 0)/(max(fitted_curve) - min(fitted_curve));
%                     significant_pk_height = all(norm_pk_height > 0.25);
%                 end
%                                   
%                 
%                 num_peaks = length(pks);
%                 score = ~(length(pks) > 1 & significant_pk_height);                

                %% Attempt 1
%                 mean_wf = mean(wf);
%                 
%                 number_waveforms = size(wf,1);
%                 
%                 xcorr_score = [];
%                 for n = 1:number_waveforms
%                     xcorr_score(n) = max( xcorr( mean_wf, wf(n,:) ) );
%                 end
%                 
%                 % Range of correlation values
%                 X = linspace(min(xcorr_score), max(xcorr_score), 100);
%                 interval = abs(X(1)-X(2));
%                 
%                 % Remove spurious activity
% %                 [hist_xcorr_score, bin_num] = histc(xcorr_score, X);
% %                 removal_hist_ind = find(hist_xcorr_score < sum(hist_xcorr_score)*0);
% %                 removal_xcorr_ind = ismember(bin_num, removal_hist_ind);
% %                 xcorr_score(removal_xcorr_ind) = [];
%                 
%                 % Number of clusters
%                 num_componments = 2; %length(unique(property_struct.sort_IDs));
%                 repeat_flag = true;
%                 
%                 gmm = {};
%                 N = 5;
%                 gmm{end+1} = gmdistribution.fit(xcorr_score', num_componments, ...
%                                              'regularize', interval/1000, ...
%                                              'covtype', 'diagonal', ...
%                                              'Replicates', N);
%                 num_componments = num_componments + 1;
%                 
%                 % Ensure each component of the GMM fit is greater than 5%
%                 % of total fit
%                 while repeat_flag,
%                     % Gaussian Mixture Model of correlation histogram (Non
%                     % Deterministic : run 10 times and pick lowest NlogL value
%                     gmm{end+1} = gmdistribution.fit(xcorr_score', num_componments, ...
%                                              'regularize', interval/1000, ...
%                                              'covtype', 'diagonal', ...
%                                              'Replicates', N);    
%                     % PDF of GMM
%                     gmm_pdf_fit = pdf(gmm{end-1}, X') * interval;
%                     
%                     % Find all local maximum
%                     [pks, locs] = findpeaks(gmm_pdf_fit);
%                     
%                     if gmm{end-1}.NlogL < gmm{end}.NlogL | length(pks) > 1 | num_componments > 4
%                         gmm = gmm{end-1};
%                         repeat_flag = false;
%                         num_componments = num_componments - 1;
%                         
%                     else
%                         num_componments = num_componments + 1;
%                         
%                     end
%                 end
%                 
%                 gmm_piecewise_fit = zeros([num_componments length(X)]);
%                 for i = 1:num_componments,
%                     gmm_piecewise_fit(i,:) = gmm.PComponents(i)*normpdf(X,gmm.mu(i),sqrt(gmm.Sigma(i))) * interval;
%                 end
%                 gmm_area = sum(gmm_piecewise_fit, 2);
%                 gmm_area_norm = gmm_area./sum(gmm_area);                  
% 
%                 % Find which Gaussian piece each peak belongs to
%                 [~, guassian_id] = max(gmm_piecewise_fit(:,locs));
%                 
%                 % Eliminate Peaks with less than 10% of total area
%                 significant_mean_pks_mask = gmm_area_norm(guassian_id) > 0.10;
%                 
%                 % Eliminate Peaks with standard deviation greater than 10%
%                 % of range
%                 norm_std = squeeze(gmm.Sigma(guassian_id))./(max(X)-min(X));
%                 bigger_gmm_component = gmm_area_norm(guassian_id) == max(gmm_area_norm(guassian_id));
%                 significant_std_pks_mask = norm_std < 0.10 | bigger_gmm_component; 
%                 
%                 % Significant Peaks
%                 significant_pks = pks(significant_mean_pks_mask & significant_std_pks_mask);
%            
%                 % Number of local maximum
%                 number_peaks = length(unique(significant_pks));
%                 
%                 % If there is only one local maximum, metric passes
%                 num_peaks = number_peaks;
%                 score = number_peaks == 1;
                
            end
            warning ('on','all');    
                
            %% PLOT Functionality
%             if property_struct.showplots,
%                 
%                 if min_sample_point_not_met,
%                     str1(1) = {sprintf('Min # of points not meet: %d', min_sample_point_not_met)}; 
%                     str1(2) = {sprintf('Score: %d', 1)}; 
%                     text(0.05, .8, str1, 'HorizontalAlignment', 'left');
%                     
%                 else
% %                     X = plot_X;
%                     bar(X, histogram);
%                     v = axis;
%                     ylim([0 v(4).*1.25])
%                     xlim([min(X), max(X)])
% 
%                     hold on
%                     grid on
%           
%                     plot(X, fitted_curve, 'r', 'LineWidth', 2)
% 
%                     % Plot Peaks
%                     plot(X(locs), pks, '*k', 'LineWidth', 10)
%                     plot(X(locs), pks, '*y', 'LineWidth', 7)    
% 
%                     v = axis; 
%                     str1(1) = {sprintf('Score: %d', score)}; 
%                     text(abs(v(1)-v(2))*.05 + v(1), v(4)*.9, str1, 'HorizontalAlignment', 'left');
%                 end                  
                

                %% Attempt 1
%                 hold on
%                 grid on
%                                                
%                 histogram = histc(xcorr_score, X);
%                 histogram = histogram./sum(histogram);
% 
%                 bar(X, histogram);
%                 
%                 for k = 1:num_componments,
%                     line(X, gmm_piecewise_fit(k,:), 'Color', [rand rand rand], 'LineWidth', 2);
%                 end
%                 
%                 plot(X, gmm_pdf_fit, 'r', 'LineWidth', 2)
%                 plot(X(locs), pks, '*g')     
%                      
%                 v = axis;
%                 str1(1) = {sprintf('Peaks with >10%% of Total Area: %d', number_peaks)}; %, 
%                 text(abs(v(1)-v(2))*.05 + v(1), v(4)*.925, str1, 'HorizontalAlignment', 'left');
                
%             end
        end
        
    end
end

