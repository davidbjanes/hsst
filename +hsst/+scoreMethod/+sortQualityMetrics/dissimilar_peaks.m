
% sortQualityMetric for statistically similar peaks for sortQualityClass
classdef dissimilar_peaks < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'Dissimilar Peaks';
        verbose_label = 'Dissimilar Peak Values'
    end
    
    methods (Static)

        function [score dprime_raw] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % Min/Max Peak Distribution Statistical Difference
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val

            %% VARIABLES FOR METRIC
            MIN_SAMPLE_DISTANCE_THRESHOLD = 3;
            MIN_DPRIME_DISTANCE_THRESHOLD = 1;
            
            if isfield(property_struct, 'wf') && ...
               isfield(property_struct, 'sort_IDs') && ...
               isfield(property_struct, 'cur_ID'),
            	
                all_wf   = property_struct.wf;
                sort_IDs = property_struct.sort_IDs;
                cur_ID   = property_struct.cur_ID;
                
                unitID_list   = unique(sort_IDs);
                
                num_waveforms = tabulate(sort_IDs);
                num_waveforms = num_waveforms(:,2);
            end
            
            %% METRIC ANALYSIS CODE

            if ~exist('unitID_list', 'var') || length(unitID_list) <= 1,
                score = 1;
                dprime_raw = [];
%                 warning('sortQualityMetric.dissimilar_peaks -> ONLY ONE UNIT !?!')
                return;
            end            
                        
            dprime_max = [];
            dprime_min = [];
            temporal_max_pk_dist = [];
            temporal_min_pk_dist = [];
            cur_unitID_ind = find(unitID_list == cur_ID);
            mean_smpl_pnt = [];
           
            for i = 1:2,  
                
                mean_val = [];
                std_val = [];

                if i == 1, [peak_val, peak_ind] = max(all_wf(:,:),[],2); end
                if i == 2, [peak_val, peak_ind] = min(all_wf(:,:),[],2); end 

                for n = 1:length(unitID_list),
                    unit_indices     = sort_IDs == unitID_list(n);
                    mean_val(n)      = mean(peak_val(unit_indices));
                    std_val(n)       = std(peak_val(unit_indices));
                    avg_wf           = mean(all_wf(unit_indices,:));
                    if i == 1,
                        [~,mean_smpl_pnt(i,n)] = max(avg_wf);
                    else
                        [~,mean_smpl_pnt(i,n)] = min(avg_wf);
                    end
%                     mean_smpl_pnt(i,n) = mode(peak_ind(unit_indices));
                    
%                     tabulate(peak_ind(unit_indices));
                    
                end

                oversorted_peaks = zeros([1 length(unitID_list)]);
                similar_peak_ind = zeros([1 length(unitID_list)]);

                % D'  (d-prime)
                dprime     = [];
                peak_dist  = [];
                precision = 0.1;
                for n = 1:length(unitID_list)
                    dprime(n) = abs(mean_val(cur_unitID_ind) - mean_val(n)) / ...
                                sqrt(0.5 * (std_val(cur_unitID_ind)^2 + std_val(n)^2));
                    dprime(n) = round(dprime(n) / precision) * precision;

                    peak_dist(n) = abs(mean_smpl_pnt(i,cur_unitID_ind) - mean_smpl_pnt(i,n));
                
                end
                               
                oversorted_peaks(dprime <= MIN_DPRIME_DISTANCE_THRESHOLD) = 1; 
                oversorted_peaks(cur_unitID_ind) = 0;

                similar_peak_ind(peak_dist <= MIN_SAMPLE_DISTANCE_THRESHOLD) = 1;
                similar_peak_ind(cur_unitID_ind) = 0;

                if i == 1, 
                    max_peak_good_units = ~(oversorted_peaks & similar_peak_ind); 
                    dprime_max           = dprime;
                    temporal_max_pk_dist = peak_dist;
                elseif i == 2, 
                    min_peak_good_units = ~(oversorted_peaks & similar_peak_ind); 
                    dprime_min           = dprime;
                    temporal_min_pk_dist = peak_dist;
                end
            end
            
            overlap_score = max_peak_good_units | min_peak_good_units;
            dprime_raw = find(overlap_score == 0);
            score = all(overlap_score);
                
%             if ~all(overlap_score)
%                 if any(num_waveforms(cur_unitID_ind) < num_waveforms(~overlap_score))
%                     score = 0;
%                 else
%                     score = 1;
%                 end
%             else
%                 score = 1;
%             end
                
%             score = all(overlap_score);
            

            %% PLOT Functionality
            if property_struct.showplots,
                
                hold on
                grid on
                
                num_units = length(unitID_list);
                
                bar([dprime_max; dprime_min]')
                plot([0.5, num_units+0.5], [1 1], 'r--', 'LineWidth', 3) 
                     
                ylabel('D prime Value')
                xlabel('Unit ID')
                xlim([0.5, num_units+0.5])
                ylim([0 2])
                set(gca, 'XTick', [1:num_units])
                set(gca, 'XTickLabel', strsplit(sprintf('(%d,%d) ', mean_smpl_pnt)))
                title(sprintf('D'' Distance Between Unit %d, and All Others', cur_ID))
                legend('Max Peak', 'Min Peak')
            end
            
            % OLD VERSION
            if false,                
                COLOR = [0 1 0; ... green
                         0 0 1; ... blue
                         1 0 0; ... red
                         0 1 1; ... cyan
                         1 0 1; ... magenta
                         1 1 0; ... yellow
                         rand([abs(length(unitID_list)-6) 3])]';  
                COLOR = COLOR(:,[1:length(unitID_list)]);
                     
                figure
                
                X = linspace(min(min(all_wf)), max(max(all_wf)), 500);
                mean_val = [];
                std_val = [];
                figures = [];
                for i = 1:4,
                	subplot(2,2,i)
                    if (i > 2)
                        set(gca, 'color', [0 0 0])
                        set(gca, 'XColor', [0.2 0.2 0.2]);
                        set(gca, 'YColor', [0.2 0.2 0.2]);
                    end
                    hold on
                    grid on
                    
                    for unit = unitID_list,
                        if mod(i,2) == 0, pk = + max(all_wf(sort_IDs == unit, :),[],2); end
                        if mod(i,2) == 1, pk = + min(all_wf(sort_IDs == unit, :),[],2); end

                        mean_val(length(mean_val)+1) = mean(pk);
                        std_val(length(std_val)+1) = std(pk);
                        
                        if i > 2,
                            plot(X, normpdf(X,mean_val(end),std_val(end)), ...
                                 'color', COLOR(:,unit), 'linewidth', 2);
                        else
                            hist(pk, X)
                        end
                        
                        if mod(i,2) == 1,
                            xlim([min(X) 0]);
                        else
                            xlim([0 max(X)]);
                        end
                    end
                    
                    if i < 3,
                        for j = 1:2,
                            subplot(2,2,j)
                            grid on
                            legend(strread(num2str(unitID_list),'%s'));

                            children = get(gca, 'child')';
                            for unit = 1:length(children),
                                set(children(unit), ...
                                    'FaceColor', COLOR(:,end-unit+1), ...
                                    'FaceAlpha', 0.5, ...
                                    'EdgeColor', COLOR(:,end-unit+1) );
                            end
                        end
                    end
                end
    
            end
            
            %% New Metric
            if false,
                
                num_wf = size(wf,1);
                num_sample_pnts = size(wf,2);
                num_units = length(unitID_list);
%                 sort_IDs
%                 cur_ID           
   
                mean_val = [];
                std_val = [];
                for n = 1:num_units
                    unit_indices    = sort_IDs == unitID_list(n);
                    mean_val(n,:)   = mean(all_wf(unit_indices,:));
                    std_val(n,:)    = std(all_wf(unit_indices,:));
                end
                
              
 
                % D'  (d-prime)
                dprime     = [];
                precision = 0.1;
                cur_unitID_ind = find(unitID_list == cur_ID);
                for n = 1:num_units   
                    dprime(n,:) = abs(mean_val(cur_unitID_ind,:) - mean_val(n,:)) ./ ...
                                  sqrt(0.5 * (std_val(cur_unitID_ind,:).^2 + std_val(n,:).^2));
                    dprime(n,:) = round(dprime(n,:) ./ precision) .* precision;
                end

                


                % PLOT
                figure
                hold all
                grid on
                for n = 1:num_units,
                    plot(sum(dprime(:,:),2)./num_sample_pnts, '*', 'linewidth', 2)
                end
%                 plot([1,num_sample_pnts], [MIN_DPRIME_DISTANCE_THRESHOLD,MIN_DPRIME_DISTANCE_THRESHOLD], 'r--', 'linewidth', 3)
%                 legend(strread(num2str(unitID_list'),'%s'));
                
            end
            
                
            
            %% Test New Algorithm
%             if false,
% 
%                 cur_unitID_ind = unitID_list == cur_ID;
%                 
%                 for n = 1:length(unitID_list),   
%                     mean_wf(n,:) = mean(all_wf(sort_IDs == unitID_list(n),:));
%                 end
% 
%                 number_waveforms = size(all_wf,1);
%                 xcorr_score = zeros([1 number_waveforms]);
%                 auto_corr_max = xcorr( mean_wf(cur_unitID_ind,:) );
%                 for n = 1:number_waveforms
%                     norm_cur_wf = all_wf(n,:);
%                     norm_cur_wf = (norm_cur_wf - min(norm_cur_wf))/(max(norm_cur_wf) - min(norm_cur_wf));
%                     
%                     norm_mean_wf = mean_wf(cur_unitID_ind,:);
%                     norm_mean_wf = (norm_mean_wf - min(norm_mean_wf))/(max(norm_mean_wf)-min(norm_mean_wf));
%                     
% %                     xcorr_score(n) = max( xcorr(norm_cur_wf, norm_mean_wf) );
% %                     xcorr_score(n) = max( xcorr(mean_wf(cur_unitID_ind,:), all_wf(n,:)) );
%                     
% %                     xcorr_score(n) = sum( [norm_cur_wf - norm_mean_wf].^2 );
%                     xcorr_score(n) = sum( [mean_wf(cur_unitID_ind,:) - all_wf(n,:)].^2 );
%                 end
%                       
%                 unit_mean_xcorr = zeros([1 length(unitID_list)]);
%                 unit_std_xcorr  = zeros([1 length(unitID_list)]);
%                 for n = 1:length(unitID_list),          
%                     ind = sort_IDs == unitID_list(n);
%                     
%                     unit_mean_xcorr(n) = mean( xcorr_score(ind) );
%                     unit_std_xcorr(n)  = std( xcorr_score(ind) );
%                     
%                 end
%                 
%                 for n = 1:length(unitID_list),        
%                     dprime(n) = abs(unit_mean_xcorr(cur_unitID_ind) - unit_mean_xcorr(n)) / ...
%                                 sqrt(0.5 * (unit_std_xcorr(cur_unitID_ind)^2 + unit_std_xcorr(n)^2));
%                 end
%                 
%                 
%                 % plot
%                 X = linspace(min(xcorr_score), max(xcorr_score), 100);
%                 COLOR = [0 1 0; ... green
%                          0 0 1; ... blue
%                          1 0 0; ... red
%                          0 1 1; ... cyan
%                          1 0 1; ... magenta
%                          1 1 0; ... yellow
%                          rand([abs(length(unitID_list)-6) 3])]';  
%                 COLOR = COLOR(:,[1:length(unitID_list)]);
%                 
%                 figure
%                 pos = get(gcf, 'Position');
%                 set(gcf, 'Position', pos.*([1 1 2 1]));
%                 subplot(1,2,1)
%                 set(gca, 'color', [0 0 0])
%                 set(gca, 'XColor', [0.2 0.2 0.2]);
%                 set(gca, 'YColor', [0.2 0.2 0.2]);
%                 hold on
%                 grid on
% 
%                 for n = 1:length(unitID_list),
%                     if cur_ID == unitID_list(n), lw = 3; else, lw = 2; end
%                     if dprime(n) > 1, format = ':'; else, format = '-'; end
%                                        
%                     hist(xcorr_score(sort_IDs == unitID_list(n)), X); hold on;
%                     histogram = hist(xcorr_score(sort_IDs == unitID_list(n)), X);
%                     
%                     children = get(gca, 'child')';
%                     set(children(1), ...
%                         'FaceColor', COLOR(:,n), ...
%                         'FaceAlpha', 0.5, ...
%                         'EdgeColor', COLOR(:,n) );
%                     
%                     plot(X, normpdf(X, unit_mean_xcorr(n), unit_std_xcorr(n))*sum(histogram)/10, ...
%                              'color', COLOR(:,n), 'linestyle', format, 'linewidth', lw);
%                 end
%                 
%                 subplot(1,2,2)
%                 set(gca, 'color', [0 0 0])
%                 set(gca, 'XColor', [0.2 0.2 0.2]);
%                 set(gca, 'YColor', [0.2 0.2 0.2]);
%                 hold on
%                 grid on
%                 for n = 1:length(unitID_list),
%                     if cur_ID == unitID_list(n), lw = 5; else, lw = 2; end
%                     if dprime(n) > 1, format = ':'; else, format = '-'; end
%                     plot(mean_wf(n,:), 'color', COLOR(:,n), ...
%                          'linestyle', format, 'linewidth', lw);
%                 end
%                 
%             end
%             
%             %% Test of a new metric
%             if false,
%             
%                 figure
%                 subplot(1,3,1)
%                     hold on
%                     to_plot = wf';
%                     plot(to_plot)
%                     plot(mean(to_plot,2), 'k:', 'linewidth',5)
%                     axis tight
%                 subplot(1,3,2)
%                     hold on
%                     to_plot = diff(wf');
%                     plot(to_plot)
%                     plot(mean(to_plot,2), 'k:', 'linewidth',5)
%                     axis tight
%                 subplot(1,3,3)
%                     hold on
%                     to_plot = diff(diff(wf'));
%                     plot(to_plot)
%                     plot(mean(to_plot,2), 'k:', 'linewidth',5)     
%                     axis tight
%                 
%             end
            
        end
        
    end
end