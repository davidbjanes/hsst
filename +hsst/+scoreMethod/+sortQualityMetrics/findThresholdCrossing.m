
function [pos_threshold, th_value, th_pnt] = findThresholdCrossing(wf)

    if size(wf,1) > size(wf,2), wf = wf'; end

    diff_wf = diff(wf, [], 1);
    same_direction = sign(diff_wf);
    same_direction(same_direction == -1) = 0;
    same_direction = logical(same_direction);
    
    possible_pos_threshold_pts = find(all(same_direction,2));
    possible_neg_threshold_pts = find(all(~same_direction,2));
    
    if ~isempty(possible_pos_threshold_pts)
        pos_threshold = true;   
        th_pnt = possible_pos_threshold_pts(1);
        th_value = max(wf(th_pnt,:));
        
    elseif ~isempty(possible_neg_threshold_pts)
        pos_threshold = false;
        th_pnt = possible_neg_threshold_pts(1);
        th_value = min(wf(th_pnt,:));
        
    else
        pos_threshold = [];
        th_pnt = [];
        th_value = [];
    
    end
    
end
