
% sortQualityMetric of SNR values for sortQualityClass
classdef snr < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'SNR';
        verbose_label = 'SNR Values'
    end
    
    methods (Static)

        function [score snr_val] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val
        
            %% VARIABLES FOR METRIC
            SNR_THRESHOLD = 1.5;
            NOISE_FLOOR = noiseEst;           
            
            
            %% METRIC ANALYSIS CODE

            warning ('off','all');
            
            % Distribution of max amplitudes of the snips
            max_amplitude = max(abs(wf),[],2);
            
            % Conversion into SNR values
            max_amplitude = abs(max_amplitude ./ NOISE_FLOOR);
            
            
            % Bin Histogram
            data = max_amplitude;
            length_X = round(sqrt(length(data)));  % 100;
            X = linspace(min(data), max(data), length_X);     
            histogram = histc(data, X)';
            
            N = length(histogram);
            if N > 5,
                smoothedHist = histogram;

                % 3-Sample-Median Smoothing
                for i = 2:N-1,
                    smoothedHist(i) = median(histogram(i-1:i+1));
                end
                smoothedHist([1 N]) = histogram([1 N]);         

                % Sliding Window Smoothing
                smoothedHist = smooth(smoothedHist,3)';
                
            else
                smoothedHist = histogram;
            end
            
            % Find all local maximum
            N = length(smoothedHist);
            max_locs = [];  
            min_locs = [];
            if N > 3,
                [~, max_locs] = findpeaks(smoothedHist);
                [~, min_locs] = findpeaks(-smoothedHist); 
            end    

            if ~isempty(max_locs) | ~isempty(min_locs)

                all_locs = sort([max_locs, min_locs]);

                if smoothedHist(1) < smoothedHist(all_locs(1)),
                    min_locs = [1 min_locs];
                else
                    max_locs = [1 max_locs];
                end

                if smoothedHist(N) < smoothedHist(all_locs(end)),
                    min_locs = [min_locs N];
                else
                    max_locs = [max_locs N];
                end        

            else
                if smoothedHist(1) < smoothedHist(N),
                    min_locs = 1;
                    max_locs = N;
                else
                    min_locs = N;
                    max_locs = 1;
                end

            end
            
            
            % Fit LOG GMM to data
            warning('off', 'all')    
            
            model_components = [];
            snr_val = [];
            score = 1;
                
            if length(data) > 10,
                K = length(max_locs);
                interval = X(2) - X(1);
                if K > 1,

                    % log GMM model
                    log_X = log(X);
                    log_data = log(data);

                    S.mu = log_X([max_locs])';
                    S.Sigma = ones(1,1,K) .* abs(log_X(3)-log_X(1));
                    S.PComponents = ones([K 1]).*1/K;
                    options = statset('TolFun', 1e-2); % 'MaxIter',10, 

                    % LOG space
                    gmm = gmdistribution.fit(log_data, K, ...
                                             'regularize', abs(min(diff(log_data)))/1000, ...
                                             'covtype', 'diagonal', ...
                                             'Start', S, ...
                                             'Options', options);
                    for k = 1:K,
                        gmmLOG = makedist('Lognormal', 'mu', gmm.mu(k), 'sigma', sqrt(gmm.Sigma(:,:,k)));
                        gmm_piecewise_fit = gmm.PComponents(k) * pdf(gmmLOG, X');
                        model_components(k,:) = gmm_piecewise_fit .* interval .* sum(histogram);
                    end

                    [~, sort_ind] = sort(gmm.mu);
                    model_components = model_components(sort_ind, :);

                elseif K == 1,
                    model_components = normpdf(X,mean(data),std(data)) .* interval .* sum(histogram);
                
                end
                
                model_pdf = sum(model_components, 1);
            
                [pks, locs] = findpeaks(model_pdf);
                snr_val = X(locs(pks == max(pks)));
                if isempty(snr_val), 
                    [~,ind] = max(model_pdf);
                    snr_val = X(ind);
                end
                score = (snr_val >= SNR_THRESHOLD);
            end
            
            
            %% Attempt 1
%             % Model the distribution as a two component GMM
%             interval = 0.01;
%             X = [min(max_amplitude):interval:max(max_amplitude)];
%             
%             if length(max_amplitude) > 10,
%                 % Number of clusters
%                 
%                 num_componments = 2; %length(unique(property_struct.sort_IDs));
%                 % Gaussian Mixture Model of correlation histogram (Non
%                 % Deterministic : run 10 times and pick lowest NlogL value
%                 N = 10;
%                 mu = [min(max_amplitude), mode(max_amplitude)]';
%                 sigma = cat(3,1,1);
%                 p = [0.5, 0.5]';
%                 S = gmdistribution(mu,sigma,p);
%                 gmm = gmdistribution.fit(max_amplitude, num_componments, ...
%                                          'regularize', interval/1000, ...
%                                          'covtype', 'diagonal', ...
%                                          'Replicates', N);
% 
%                 % PDF of GMM
%                 model_pdf = pdf(gmm, X') * interval;
% 
%                 gmm_piecewise_fit = zeros([num_componments length(X)]);
%                 for i = 1:num_componments,
%                     gmm_piecewise_fit(i,:) = gmm.PComponents(i)*normpdf(X,gmm.mu(i),sqrt(gmm.Sigma(i)));
%                 end
%                 gmm_area = sum(gmm_piecewise_fit, 2);
%                 gmm_area_norm = gmm_area./sum(gmm_area);                  
%                                        
%                 [~, ind] = max(max(gmm_piecewise_fit, [], 2));
% 
%                 % Find all local maximum
%                 [max_pks, max_locs] = findpeaks(model_pdf);
%                 [min_pks, min_locs] = findpeaks(-model_pdf);
%                 
%                 significant_pk_area = 0;
%                 significant_pk_height = 0;
%                 if length(max_pks) > 1,
%                     valley_found = max_locs(1) < min_locs & min_locs < max_locs(2);
%                     
%                     if valley_found,
%                         peak_height = max_pks + min_pks;
%                         norm_peak_height = (peak_height - 0)/(max(model_pdf) - min(model_pdf));
%                         
%                         significant_pk_area = all(gmm_area_norm > 0.25);
%                         significant_pk_height = all(norm_peak_height > 0.5);
%                     end
%                 end
% 
%                 model_pdf = model_pdf .* sum(histogram);
%                                     
%                 % Select the mean SNR value as the greater of the two peaks
%                 snr_val = gmm.mu(ind);
%                 
%             else
%                 % Mean Value of Distribution
%                 snr_val = mean( max_amplitude );
%                 significant_pk_area = 0;
%                 significant_pk_height = 0;    
%                 
%             end
%             
%             score = (snr_val >= SNR_THRESHOLD); % & ~(significant_pk_area | significant_pk_height);
%             
%             warning ('on','all');
            
            
            %% PLOT Functionality
            if property_struct.showplots,
                
%                 histogram = histc(max_amplitude, X);
                bar(X, histogram);
                hold on
                grid on
                
                axislim = axis;
                xlim([1 axislim(2)])
                
                plot(SNR_THRESHOLD, linspace(0, axislim(4), 100), 'r--', 'LineWidth', 2)
                plot(snr_val, linspace(0, axislim(4), 100), 'g:', 'LineWidth', 2)
                plot(X, model_pdf, 'c', 'LineWidth', 2)

                xlabel('SNR')
                ylabel('Histogram of |Peak Val|')
                
                set(gcf, 'Renderer', 'zbuffer');

            end
            
                       
        end
        
    end
end

