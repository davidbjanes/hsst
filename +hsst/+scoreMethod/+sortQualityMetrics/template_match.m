
% sortQualityMetric called template_match for sortQualityClass
classdef template_match < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'Template Match';
        verbose_label = 'Template Match Score'
    end
    
    methods (Static)

        function [score temp_match] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val
            
            import hsst.scoreMethod.sortQualityMetrics.*
            
            %% VARIABLES FOR METRIC
            NSIGMA = 1;     % threshold at which residuals are examined            
            MIN_PERCENT_OF_HISTOGRAM_PER_BIN = 0.01;    
            MIN_NUM_CONSECUTIVE_BINS_THRESHOLD = 3;
            DEPTH_OF_VALLEY_PERCENT = 0.40;
            DEBUG_PLOT_ERROR_FLAG = false;
            
            [~, ~, ALIGNMENT_POINT] = findThresholdCrossing(property_struct.wf);
            ALIGNMENT_POINT = ALIGNMENT_POINT + 1;
            
            
            %% METRIC ANALYSIS CODE
            
            p = normcdf([-NSIGMA NSIGMA]);
            template_threshold = p(2)-p(1);
            
            num_samples = size(wf,2); % Waveform Sample Length
            num_waveforms = size(wf,1); % Number of waveforms 

            mean_wf = mean(wf);
            detrended_wf = abs(bsxfun(@minus, wf, mean_wf));
            var_waveform = std( wf ,1 );
            
            normbyVar_rec_det_wf =  bsxfun(@rdivide, detrended_wf, var_waveform);
            percent_points_inside = mean(normbyVar_rec_det_wf' < NSIGMA, 2);
            
            points_below_threshold = find(percent_points_inside < template_threshold)';
            points_below_threshold(ismember(points_below_threshold, ALIGNMENT_POINT+[-1 0 1])) = [];
                          
            if ~isempty(points_below_threshold),
                
                % Sort low to high
                if ~property_struct.showplots,
                    [~, ind] = sort(percent_points_inside(points_below_threshold));
                    points_below_threshold = points_below_threshold(ind);
                end
               
                % For each point below the Gaussan Prediction
                bimodal_array = false([1 num_samples]);
                for pt_ind = [1:num_samples],
                                             
                    % Evaluate if waveform slice is bimodal
                    sample_points = wf(:,pt_ind);
                    bimodal_array(pt_ind) = bimodalityDetector(sample_points, DEBUG_PLOT_ERROR_FLAG);
                    
                    % Stop as soon as bimodality is found
                    if bimodal_array(pt_ind) & ~property_struct.showplots,
                        break;
                    end
                end
                              
                score = all(bimodal_array(points_below_threshold) == false);
                temp_match = bimodal_array;
                
            else
                bimodal_array = boolean(zeros([1 num_samples]));
                score = 1;
                temp_match = -1;            
            end
                                   
            
            %% PLOT Functionality
            if property_struct.showplots,

%                 resized_wf = wf';
%                 vec_wf = resized_wf(:)';
%                 wf_pair = [vec_wf; repmat([1:num_samples], 1, num_waveforms)];
% 
%                 xygrid = [round(sqrt(num_waveforms))*2 num_samples];
%                 binCount = hist3(wf_pair',xygrid);
% 
%                 max_y_val = max(vec_wf);
%                 min_y_val = min(vec_wf);
%                 Y = linspace( min_y_val, max_y_val, xygrid(1) );
%                 X = [1:xygrid(2)];
% 
%                 imagesc(X, (Y), (binCount));
%                 set(gca,'YDir','normal');
%                 colormap(gray)
%                 grid on
%                 set(gca, 'XTick', [0:num_samples]+0.5);
%                 set(gca, 'XTicklabel', []);
%                 
%                 pnts_handle = [];
%                 for n = 1:num_samples,
%                     pnts_handle(n) = line(n, mean_wf(n));
%                     userdata.wf = wf(:,n);
%                     userdata.sample_num = n;
%                     
%                     set(pnts_handle(n), 'userdata', userdata );
%                     
%                     if bimodal_array(n),
%                         set(pnts_handle(n), 'color', 'r', 'lineStyle', '*', 'linewidth', 2);
%                     else
%                         set(pnts_handle(n), 'color', 'k', 'lineStyle', '*');
%                     end                   
%                 end
%                 
%                 set(pnts_handle, 'ButtonDownFcn', @template_match.plotCallBack);
%                 
%                 v = axis; 
%                 str1(1) = {sprintf('Score: %d', score)}; 
%                 text(abs(v(1)-v(2))*.05 + v(1), v(4)*.9, str1, 'HorizontalAlignment', 'left');    
            
                
            %% Attempt 1
%                 figure
%                     subplot(2,1,1)
%                     set(gcf, 'renderer', 'zbuffer')
%                     fig_pos = get(gcf, 'Position');
%                     set(gcf, 'Position', [fig_pos(1), fig_pos(2), fig_pos(3)/2, fig_pos(4)]);
                    hold on
                    grid on
 
                % Plot Data
                h1 = plot(linspace(1,num_samples,num_samples), ...
                          ones([1 num_samples]) * template_threshold, 'r--', 'linewidth', 2);
                h2 = plot(percent_points_inside, 'b*-');
                                
                % Plot Axes
                axis tight
                ylim([0 1]);
                set(gca, 'YTick', [0, 0.3829, 0.6827, 0.8664, 0.9545]); %, 0.9973
                set(gca, 'YTickLabel',{'0', '38', '68', '87', '95'});
                ylabel('Percent [%]')
                
                set(gca, 'XTick', [1:10:31])
                set(gca, 'XTickLabel', strread(num2str([-10:10:21]/24, 2),'%s'))
%                 xlabel('Time')
                
                % Text
                v = axis; 
                str1(1) = {sprintf('Score: %d', score)}; 
                text(abs(v(1)-v(2))*.05 + v(1), v(4)*.9, str1, 'HorizontalAlignment', 'left');    

                
                % Plot Analysis
                pnts_handle = [];
                for n = 1:num_samples,
                    pnts_handle(n) = line(n, percent_points_inside(n));
                    userdata.wf = wf(:,n);
                    userdata.sample_num = n;
                    
                    set(pnts_handle(n), 'userdata', userdata );
                    
                    if bimodal_array(n),
                        set(pnts_handle(n), 'color', 'r', 'lineStyle', '*', 'linewidth', 2);
                    else
                        set(pnts_handle(n), 'color', 'g', 'lineStyle', '*');
                    end                   
                end
                
%                 title('Unit 1')
%                 legend([h1(1), pnts_handle(find(bimodal_array, 1, 'first')), pnts_handle(find(~bimodal_array, 1, 'first'))], ...
%                     'Ideal Sigma', 'Multimodal', 'Unimodal')
%                 set(gca,'FontSize',13)
%                 set(findall(gcf,'type','text'),'fontSize',13)
                set(pnts_handle, 'ButtonDownFcn', @template_match.plotCallBack);
                
            end
            
        end
        
        function plotCallBack( src, ~ )

            import hsst.scoreMethod.sortQualityMetrics.*

            userdata = get(src, 'userdata');
                    
            figure
            
            bimodalityDetector(userdata.wf, true);
                    
            title(sprintf('Sample Number %d', userdata.sample_num));
            
        end
               
    end
end
