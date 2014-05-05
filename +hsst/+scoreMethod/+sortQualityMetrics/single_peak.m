
% sortQualityMetric of SNR values for sortQualityClass
classdef single_peak < hsst.scoreMethod.sortQualityMetrics.metric
    
    properties (Constant)
        label = 'Single Peak';
        verbose_label = 'Single Peak'
    end
    
    methods (Static)

        function [score, num_critial_pts] = generateMetricScore(wf, ts, noiseEst, property_struct)
        % NxD       - wf                             
        % Nx1       - ts                              
        % ( * )     - property_struct.*
        % boolean   - score 
        % double    - snr_val         
            
            import hsst.scoreMethod.sortQualityMetrics.*

            %% VARIABLES FOR METRIC
            DEBUG_PLOT_ERROR_FLAG = false;
            
            
            %% METRIC ANALYSIS CODE

            warning ('off','all');
            
            mean_waveform = mean(wf);
            
            num_samples = size(wf,2);
            num_waveforms = size(wf,1);
            
            score = 0;
            num_critial_pts = -1;

            % Identify Shoulders in Data
            if num_samples > 5,
                
                shoulder_loc = single_peak_detectShoulder(mean_waveform, false);
                                
                score = sum(shoulder_loc) > 1;
                num_critial_pts = sum(shoulder_loc);
            end
            
            warning ('on','all');   
            
            
            %% PLOT Functionality
            if property_struct.showplots,
                
                set(gcf, 'Renderer', 'zbuffer');
                hold all
                grid on
                plot(wf', 'color', [0.7 0.7 0.7])                
                plot(mean_waveform, 'k--', 'linewidth', 5)
                
                plot(find(shoulder_loc), mean_waveform(shoulder_loc), 'c*', 'linewidth', 5) 

                xlabel('Sample Number')
                ylabel('Amplitude')
                
                axis tight            
                
            end
            
            if DEBUG_PLOT_ERROR_FLAG,
                single_peak_detectShoulder(mean_waveform, DEBUG_PLOT_ERROR_FLAG)
                
                subplot(2,1,1)
                    set(gcf, 'Renderer', 'zbuffer');
                    hold all
                    grid on
                    plot(wf', 'color', [0.7 0.7 0.7])                
                    plot(mean_waveform, 'k--', 'linewidth', 5)

                    plot(find(shoulder_loc), mean_waveform(shoulder_loc), 'c*', 'linewidth', 5) 
                    set(gca, 'XTick', [1:5:num_samples]);
                    axis tight
                    
                    xlabel('Sample Number')
                    ylabel('Amplitude')
                    
%                 fig_pos = get(gcf, 'Position');
%                 set(gcf, 'Position', [fig_pos(1) fig_pos(2), fig_pos(3)/2, fig_pos(4)]);
                    
            end            
                       
        end
       
        
    end
end

