
% sortQualityMetric for statistically similar peaks for sortQualityClass
classdef mean_wf_similarity < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'DPrime Dissimilarity';
        verbose_label = 'DPrime Dissimilarity'
    end
    
    methods (Static)

        function [score dprime_raw] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % Min/Max Peak Distribution Statistical Difference
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val

            import hsst.scoreMethod.sortQualityMetrics.*
            import stdlib.*
            warning('off', 'all')
            
            %% VARIABLES FOR METRIC
            MIN_DPRIME_DISTANCE_THRESHOLD = 1;
            PCNT_TO_KEEP = 0.99;
            
            if isfield(property_struct, 'wf') && ...
               isfield(property_struct, 'sort_IDs') && ...
               isfield(property_struct, 'cur_ID'),
            	
                all_wf   = property_struct.wf;
                sort_IDs = property_struct.sort_IDs;
                cur_ID   = property_struct.cur_ID;
                
                unitID_list   = unique(sort_IDs);
            end
            
            %% METRIC ANALYSIS CODE
                
            num_sample_pnts = size(all_wf,2);
            num_wf = size(all_wf,1);
            num_units = length(unitID_list);
         
            if num_units < 2,
                score = 1;
                dprime_raw = [];
                
            else
            
                % Grab Unit Indices
                unitIndices = bsxfun(@eq, unitID_list', sort_IDs);
                cur_unit_ind = find(unitID_list == cur_ID);
                
                mean_waveform = mean(all_wf(unitIndices(cur_unit_ind,:),:));
                
                % Calculate Sum of Squared Error
                x = all_wf(:,:);
                y = mean_waveform;
                squared_error = abs(bsxfun(@minus, x, y)).^2;
                SSE = sum(squared_error, 2);
                                              
                % Get GammaFit pdf of each SSE(unit)
                X = linspace(min(SSE), max(SSE), sqrt(num_wf));
                norm_cumsum_hist = cumsum(hist(SSE, X))/num_wf;
                
                end_ind = find(norm_cumsum_hist > PCNT_TO_KEEP, 1, 'first');
                X = linspace(min(SSE), X(end_ind), 2*sqrt(num_wf));
                interval = abs(X(1)-X(2));
                
                gamma_PDF = [];
                mu_est = [];
                for unitX_ind = 1:num_units,
                    
                    data = [SSE(unitIndices(unitX_ind,:))];
                    mu_est(unitX_ind) = mean(data);
                    std_est(unitX_ind) = std(data);
                    
                    if property_struct.showplots,
                        min_data = min(data) - 0.001*interval;
                        shifted_data = data - min_data;

                        mean_subracted_shifted_data = abs(shifted_data - mean(shifted_data));
                        outliers_removed_shifted_data = shifted_data(mean_subracted_shifted_data < 7*std(shifted_data));

                        phat = gamfit(outliers_removed_shifted_data);
                        gamma_PDF(unitX_ind,:) = gampdf(X,phat(1),phat(2));

                        shift_ind = find(X > min_data, 1, 'first');
                        gamma_PDF(unitX_ind,:) = [zeros([1 shift_ind]) gamma_PDF(unitX_ind,[1:end-shift_ind])]; 

                        [~, max_ind] = max(gamma_PDF(unitX_ind,:));
                    end
                end
                
                
                % Check bimodality when combining each unit with curUnit
                bimodal = true;
                db = []; db(cur_unit_ind) = 0; 
                for unitX_ind = 1:num_units,

                    if unitX_ind ~= cur_unit_ind,
                        x1 = SSE(unitIndices(unitX_ind,:));
                        x2 = SSE(unitIndices(cur_unit_ind,:));
                             
                        dprime = abs(mu_est(cur_unit_ind) - mu_est(unitX_ind)) ./ ...
                                     sqrt(0.5 * (std_est(cur_unit_ind).^2 + std_est(unitX_ind).^2));

                        % Round
                        precision = 0.05;
                        round_dprime = round(dprime ./ precision) .* precision;
                        bimodal(unitX_ind) = round_dprime > 1;
                        
%                         db(unitX_ind) = bhattacharyya(x1, x2);
%                         bimodal(unitX_ind) = db(unitX_ind) > 0.25;

%                         data_combined = [ x1, x2 ];
%                         bimodal(unitX_ind) = bimodalityDetector(data_combined, false);
%                         
%                         if ~bimodal(unitX_ind),
%                             figure
%                             bimodalityDetector(data_combined, true);
%                         end
                    
                    else
                        bimodal(unitX_ind) = true;
                    end
                end
                
                % Output
                score = all(bimodal);
                dprime_raw = find(~bimodal);
                
                warning('on', 'all')
                
%% Attempt 1
%                 mean_val = [];
%                 std_val = [];
%                 for n = 1:num_units
%                     unit_indices    = sort_IDs == unitID_list(n);
%                     mean_val(n,:)   = mean(all_wf(unit_indices,:));
%                     std_val(n,:)    = std(all_wf(unit_indices,:));
%                 end
% 
% 
%                 % D'  (d-prime)
%                 dprime     = [];
%                 cur_unitID_ind = find(unitID_list == cur_ID);
%                 for n = 1:num_units   
%                     dprime(n,:) = abs(mean_val(cur_unitID_ind,:) - mean_val(n,:)) ./ ...
%                                   sqrt(0.5 * (std_val(cur_unitID_ind,:).^2 + std_val(n,:).^2));
%                 end
% 
%                 
%                 % Round
%                 precision = 0.05;
%                 average_dprime = sum(dprime(:,:),2)./num_sample_pnts;
%                 average_dprime = round(average_dprime ./ precision) .* precision;
% 
%                 
%                 passed_metric = average_dprime >= MIN_DPRIME_DISTANCE_THRESHOLD;
%                 passed_metric(cur_unitID_ind) = true;
%                 
%                 score = sum(passed_metric) == num_units;
%                 dprime_raw = find(~passed_metric)';
%                 
%                 
%                 % Sort dprime raw based on similiarity
%                 if ~isempty(dprime_raw)
%                     [~, sort_ind] = sort(average_dprime(dprime_raw));
%                     dprime_raw = dprime_raw(sort_ind);
%                 end
            
            end


            %% PLOT Functionality
            if property_struct.showplots && num_units > 1,
            
%                 figure
% fig_pos = get(gcf, 'Position');
% set(gcf, 'Position', [fig_pos(1), fig_pos(2), fig_pos(3)/2, fig_pos(4)/2]);
                hold all
                grid on
                unitY_ind = cur_unit_ind;
                curves = [];
                for unitX_ind = 1:num_units,

                    distributions = SSE(unitIndices(unitX_ind,:));

                    histogram = histc(distributions,X);

%                         curves(unitX_ind,:) = smooth(histogram,25);
                    handle(unitX_ind) = plot(X, histogram, '-', 'linewidth', 1);

                    norm_gamma_PDF = gamma_PDF(unitX_ind,:) * interval * sum(histogram);
                    plot(X, norm_gamma_PDF, '-', ...
                                            'color', get(handle(unitX_ind), 'color'), ...
                                            'linewidth', 2)
                end
                legend(handle, strread(num2str(unitID_list),'%s'))

                axis tight

%                 v = axis; 
%                 str1(1) = {sprintf('Score: %d', score)}; 
%                 str1(2) = {sprintf('Overlaps With: %s', num2str(dprime_raw))}; 
%                 str1(3) = {sprintf('Db: %s', num2str(db,2))}; 
%                 text(abs(v(1)-v(2))*.05 + v(1), v(4)*.9, str1, 'HorizontalAlignment', 'left');    
                

%                 if num_units > 1,
%                     hold all
%                     grid on
%                     
%                     for n = [2,6]   
%                         if passed_metric(n),
%                             plot(dprime(n,:))
%                         else
%                             plot(dprime(n,:), 'linewidth', 2)
%                         end
%                     end       
%                     
%                     ylim([0 5])
                    
%                     bar(average_dprime)
%                     plot(find(~passed_metric), average_dprime(~passed_metric), 'r*', 'linewidth', 2)
% 
%                     set(gca, 'XTick', [1:num_units]);
%                     xlabel('Units')
%                     ylabel('Averaged D'' Value Across All Samples')
% 
%                     plot([0, num_units+1], [MIN_DPRIME_DISTANCE_THRESHOLD, MIN_DPRIME_DISTANCE_THRESHOLD], 'r--', 'linewidth', 3)
%                     set(gca, 'XTickLabel', unitID_list);
%                 end
                
            end          
        end
        
    end
end