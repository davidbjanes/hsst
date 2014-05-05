
% Abstract class to interface to sortQualityClass
classdef metric < handle
    % Interface for 
    
    properties (Abstract, Constant)
        label
        verbose_label
    end
    
    % Abstract Instance of the Metric Method
    methods (Abstract, Static)
        
        [bool_score raw_score] = generateMetricScore(wf, ts, noiseEst, property_struct);
                
    end
    
    methods 
        
        function plot_metric(this, wf_n, ts_n, noiseEst, property_struct, handle_axis)
            if isempty(handle_axis)
                handle_figure = figure;
                handle_axis   = gca;

            else
                obj_axis      = get(handle_axis);
                handle_figure = obj_axis.Parent;

            end

            figure(handle_figure);
            set(handle_figure, 'CurrentAxes', handle_axis); 
            cla(gca, 'reset');
                
            if property_struct.showplots,
                this.generateMetricScore(wf_n, ts_n, noiseEst, property_struct);
            else
                property_struct.showplots = true;
                this.generateMetricScore(wf_n, ts_n, noiseEst, property_struct);
                property_struct.showplots = false;
            end
            
        end
        
    end

    
    methods
        
        function [score raw_score] = scoreAllUnits(this, wf, ts, noiseEst, ...
                                                   sort_IDs, property_struct)
        % sortQualityMetric - this 
        % NxD               - wf                             
        % Nx1               - ts                             
        % boolean           - property_struct.verbose        
        % boolean           - property_struct.export_plot    
        % ( * )             - property_struct.*
        % Mx1               - score (M = number of units)
               
        unitID_list   = unique(sort_IDs);
                
        raw_output     = {};
        bool_output    = zeros([1 length(unitID_list)]);
        
        if property_struct.showplots, figure; end            
            
        for n = 1:length(unitID_list), 
            
            if property_struct.showplots, subplot(length(unitID_list), 1, n); end
            
            unit_index = (sort_IDs == unitID_list(n));
            
            wf_n = wf(unit_index, :);
            ts_n = ts(unit_index);
            
            property_struct.cur_ID = unitID_list(n);
            
            [bool_output(n), raw_output{n}] = this.generateMetricScore(wf_n, ts_n, noiseEst, property_struct);
            
        end
        
        score     = bool_output;
        raw_score = raw_output;

        if property_struct.verbose, fprintf([this.verbose_label ': %s\n'], num2str(snr_val)); end
        
        end
        
    end
end

