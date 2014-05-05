function plotWaveforms(snips, sortCode, V_thresh, u_k)
    
    if nargin < 4,
        V_thresh = [];
        u_k = [];
    elseif nargin < 5,
        u_k = [];
    end 
    
    K = length(sortCode);
    
    figure
        hold all
        grid on
        
        format = {'m', 'r', 'b', 'g', 'c', 'y', };
        format_dots = {'g', 'y', 'm', 'r', 'b', 'c'};
        
        if ~isempty(V_thresh)
            plot([1:size(snips,2)], ones([size(snips,2) 1])*V_thresh, 'r:') 
        end
        
        for k = unique(sortCode),
            if ~isempty(snips(sortCode == k, :)),
                line_color = format{mod(k, length(format)) + 1};
                if k == 0, 
                    line_color = 'w'; 
                end
                plot([1:size(snips,2)], snips(sortCode == k, :)', 'Color', line_color);
            end
        end

        for k = 1:K,
            if ~isempty(u_k),
                plot([1:size(snips,2)], u_k(:,k), '*', ...
                        'Color', format_dots{k}, 'LineWidth', 3)
            end
        end        
        
        xlabel('Sample Number')
        ylabel('Amplitude [uV]')
        title('Voltage vs Time of Snips')
        axis tight
        set(gca, 'color', [0 0 0])
        set(gca, 'xcolor', [0.4 0.4 0.4])
        set(gca, 'ycolor', [0.4 0.4 0.4])
    
end