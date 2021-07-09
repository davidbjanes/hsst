
function [ snips, time_stamps, property ] = getSnippets( rawWaveform, sample_freq, threshold )
    
    %% Variables
    window_size = 32;
    window_align = round(window_size/3);
    
    % Estimate Threshold from raw waveform if none is set
    noise_estimate = [];
    if nargin < 3 || isempty(threshold),
        noise_estimate = hsst.extractorMethod.thresholding.getNoiseEstimate( rawWaveform, sample_freq );
        threshold = noise_estimate;
    end
    
    %% Code
    
    % Get Snips
    snips = [];
    time_stamps = [];
    n = length(rawWaveform);

    ind = 1:n;
    time_vector = 0:1/sample_freq:(n/sample_freq-1/sample_freq); 

    % threshold is positive
    if sign(threshold) > 0,
        rawWaveform = -rawWaveform;
        threshold = -threshold;
        invert_flag = true;
    else
        invert_flag = false;
    end
    
    % Grab Snips
    window_ind = [];
    for k = ind(rawWaveform < threshold),

        lower_ind = k-window_align+1;
        upper_ind = lower_ind+window_size-1;
        
        ind_within_upper_bound = upper_ind <= length(rawWaveform);
        ind_within_lower_bound = lower_ind > 0;
              
        if ind_within_upper_bound && ind_within_lower_bound, 
            positive_slope_incpt = rawWaveform(k-1) >= threshold;
            
            if positive_slope_incpt,
                snips(size(snips,1)+1, :) = rawWaveform(lower_ind:upper_ind);
                time_stamps(length(time_stamps)+1) = time_vector(k);
                
                window_ind(:, size(window_ind,2)+1) = [lower_ind, upper_ind];
            end
        end   
    end
    
    % Handle Overlapping Snippets
    window_ind = window_ind';
    x = window_ind(2:end,1)-window_ind(1:end-1, 2);
    x = sign(x) == -1;
    overlapping_snips = [0 x' 0 0] | [0 0 x' 0];
    
    if sum(overlapping_snips) > 1,
        overlapping_ind(:,1) = find(diff(overlapping_snips) == 1);
        overlapping_ind(:,2) = find(diff(overlapping_snips) == -1)-1;
        snips_to_remove = boolean(zeros([length(window_ind), 1]));
        for k = 1:size(overlapping_ind,1),

            overlapping_ind_k = overlapping_ind(k,1):overlapping_ind(k,2);

            [~, ind] = max(snips(overlapping_ind_k,:), [], 2);
            [~, ind] = min(abs(window_align - ind));
            overlapping_ind_k(ind) = [];
            snips_to_remove(overlapping_ind_k) = 1;  

        end

        snips(snips_to_remove, :)    = [];
        time_stamps(snips_to_remove) = [];
    end

    % if V_thresh had been positive, re-invert wf/V_thresh    
    if invert_flag,
        rawWaveform = -rawWaveform;
        snips = -snips;
        threshold = -threshold;
    end
    
    % Output Info
    property.threshold = threshold;
    property.noise_estimate = noise_estimate;
    property.negative_threshold = invert_flag;
    property.window_size = window_size;
    property.window_align = window_align;
        
end

