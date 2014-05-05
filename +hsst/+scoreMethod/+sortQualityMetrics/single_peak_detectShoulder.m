
% Determine Detect Shoulder 
function shoulder_loc = single_peak_detectShoulder(data, PLOT_FIGURE)

    N = length(data); 
    shoulder_loc = [];
    LEFT = 1;
    RIGHT = 2;

    if PLOT_FIGURE, figure; end

    for direction = [LEFT],
        if direction == LEFT, data = data; end
        if direction == RIGHT, data = fliplr(data); end


        all_slopes = diff(data);
        slope_adjacnt = all_slopes(1:end-1);
        slope_endpnt = bsxfun(@minus, data(1), data(2:end-1)) ./ -[1:(length(data)-2)];

        diff_btw_slopes = abs(slope_adjacnt - slope_endpnt);
        signfcnt_diff_btw_slopes = diff_btw_slopes > std(diff_btw_slopes);

        [~, slope_max_ind] = max([slope_adjacnt(1:end); slope_endpnt(1:end)]);
        slope_crossings_temp = slope_max_ind(1:end-1) + slope_max_ind(2:end) == 3;
        slope_crossings = [0 slope_crossings_temp];

        signfcnt_diff_btw_slopes(logical(slope_crossings)) = 0;
        signfcnt_groups = bwlabel(signfcnt_diff_btw_slopes);
        unsignfcnt_groups = bwlabel(~signfcnt_diff_btw_slopes);

        critial_points = [];

        num_sig_groups = length(unique(signfcnt_groups)) - 1;
        num_unsig_groups = length(unique(unsignfcnt_groups)) - 1;
        for group = 1:num_sig_groups,

            if find(unsignfcnt_groups == group, 1) < find(signfcnt_groups == group, 1),
                a_ind = find(unsignfcnt_groups == group); 
                a_ind(end+1) = a_ind(end)+1;
                slope_crossing_before = a_ind(find(slope_crossings(a_ind), 1, 'last'));

                if ~isempty(slope_crossing_before),
                    critial_points(end+1) = find(signfcnt_groups == group, 1, 'first')-1;
                end
            end

            if num_unsig_groups >= group+1,
                b_ind = find(unsignfcnt_groups == group+1); b_ind(end+1) = b_ind(1)-1;
                slope_crossing_after = b_ind(find(slope_crossings(b_ind), 1, 'first'));

                if ~isempty(slope_crossing_after),
                    critial_points(end+1) = find(signfcnt_groups == group, 1, 'last')+1;
                end         
            end
        end


        if ~isempty(critial_points), 
            critial_points(abs(diff(critial_points)) == 0) = [];

            if direction == LEFT, 
                loc = critial_points+1;
            elseif direction == RIGHT, 
                loc = (N) - critial_points;
            end

            shoulder_loc(end+1:end+length(loc)) = loc;
        end


        if PLOT_FIGURE, 
            subplot(2,1,direction+1)
            hold on
            grid on
            if direction == LEFT, X = [2:N-1]; end
            if direction == RIGHT, X = [N-1:-1:2]; end   

            plot(X, slope_adjacnt, 'r*-')
            plot(X, slope_endpnt, 'b*-')

            slope_diff = mean([slope_endpnt; slope_adjacnt]);
            plot(X(signfcnt_diff_btw_slopes), slope_diff(signfcnt_diff_btw_slopes), 'k*', 'linewidth',2)

            plot(X(critial_points), slope_diff(critial_points), 'g*', 'linewidth', 1)

            legend('Adjacent', 'EndPoint')
            if direction == LEFT, 
%                 title('Left');
            elseif direction == RIGHT, 
%                 title('Right'); 
%                         set(gca, 'XDir', 'reverse')
            end
            
            xlabel('Sample Number')
            set(gca, 'XTick', [1:5:length(data)]);
            axis tight
            xlim([1 N]);
        end            

        if direction == RIGHT, data = fliplr(data); end
    end

    % Remove concurrent/duplicate values
    shoulder_loc = sort(shoulder_loc);
    shoulder_loc(find(abs(diff(shoulder_loc)) <= 1)+1) = [];

%     if PLOT_FIGURE, 
%         subplot(3,1,2)
%         hold on
%         grid on
%         
%         plot(data, 'k-*')
%         plot(shoulder_loc, data(shoulder_loc), 'g*', 'linewidth', 2)
%         xlim([1 N])
%     end        

    % Build mask
    shoulder_mask = zeros([1 N]);
    shoulder_mask(shoulder_loc) = 1;
    shoulder_loc = logical(shoulder_mask);

end
