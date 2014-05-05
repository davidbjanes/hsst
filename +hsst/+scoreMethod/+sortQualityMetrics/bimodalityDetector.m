
% Helper Function for Metric "Template Match"
function bimodal = bimodalityDetector(data, PLOT_FLAG)
    global MIN_PERCENT_OF_HIST_PER_BIN   
    global MIN_NUM_CONSECUTIVE_BINS_THRESHOLD
    global DEPTH_OF_VALLEY_PERCENT
    global INSIG_PK_VALLEY_PAIR_DIST
    
    %% CONSTANTS
    MIN_PERCENT_OF_HIST_PER_BIN = 0.01;
    MIN_NUM_CONSECUTIVE_BINS_THRESHOLD = 3;
    DEPTH_OF_VALLEY_PERCENT = 0.40;
    INSIG_PK_VALLEY_PAIR_DIST = 0.25;
    
    
    %% CODE
    
    data = data(:);
    
    % Check that there is data to run tests on
    if length(data) < 50,
        bimodal = 0;
        return;
    end
    
    % Histogram Data
    [histogram, X, fitHistogram, fitX, fitData] = binData(data);

    % Check that there is significant data to run tests on
    if length(fitData) < 50 | length(fitX) < 3,
        bimodal = 0;
        return;
    end
    
    % Smooth Histogram
    smoothedFitHistogram = smoothHistogram(fitHistogram);

    % Check if any points dip below MIN_PERCENT_OF_HISTOGRAM_PER_BIN
    bimodalThresholdTest = checkThresholdCrossings(smoothedFitHistogram);
    
    % Check for sigificant Peak/Valley Pairs
    [min_locs, max_locs] = findLocalPeaks(smoothedFitHistogram);
    [bimodalPeaksTest, significantLocs, nvd] = findSignificantPeaks(min_locs, max_locs, ...
                                                               smoothedFitHistogram);
    
    % Check (log_)GMM model
    [min_locs, max_locs] = removeInsigPkValleyPairs(smoothedFitHistogram, ...
                                                    min_locs, max_locs);
    [bimodalFitTest, model_pdf] = fitDataWithModel(smoothedFitHistogram, fitX, ...
                                                   min_locs, max_locs, fitData);
    
    % Score
    bimodal = any(~[bimodalThresholdTest bimodalPeaksTest bimodalFitTest]);
    
    
    % Plot
    if PLOT_FLAG, %& bimodal,
        debuggingPlotter(histogram, X, smoothedFitHistogram, fitX, model_pdf, ...
                         min_locs, max_locs, significantLocs, nvd, ...
                         bimodal);
    end
end


%% Plotting Functions

function debuggingPlotter(histogram, X, ...
                          smoothedFitHistogram, fitX, ...
                          model_pdf, ...
                          min_locs, max_locs, significantLocs, norm_valley_depth, ... 
                          bimodal)
    
    set(gcf, 'Renderer', 'zbuffer');
    hold all
    grid on

    % Plot Original Data
    bar(X, histogram);
    plot(fitX, smoothedFitHistogram, 'r-', 'linewidth', 2)

    % Plot Peaks
    plot(fitX(max_locs), smoothedFitHistogram(max_locs), 'y*', 'linewidth', 2)
    plot(fitX(min_locs), smoothedFitHistogram(min_locs), 'g*', 'linewidth', 2)
    plot(fitX(significantLocs), smoothedFitHistogram(significantLocs), 'r*', 'linewidth', 2) 

    % Model PDF
    plot(fitX, model_pdf, 'c', 'linewidth', 2 )
           
    
    % Axis and Labels
    v = axis;
    ylim([0 v(4).*1.25])
    xlim([min(X), max(X)])

    v = axis;
    str1(1) = {sprintf('Depth: %s', num2str(norm_valley_depth, '%2.2f '))}; 
    str1(2) = {sprintf('Score: %d', ~bimodal)}; 
    text(abs(v(1)-v(2))*.05 + v(1), v(4)*.9, str1, 'HorizontalAlignment', 'left');                
end


%% Bimodal Test Functions

% Fits data to model 
function [bimodalFitTest, model_pdf] = fitDataWithModel(histogram, X, ...
                                                        min_locs, max_locs, data)  
        
    % Identify Std of Peaks
    sigma = [];
    for k = 1:length(max_locs),
        min_pk_L_ind = find(min_locs < max_locs(k), 1, 'last');
        min_pk_R_ind = find(max_locs(k) < min_locs, 1, 'first');

        all_peaks = X([min_pk_L_ind max_locs(k) min_pk_R_ind]);
        dist2Valley = abs(diff(all_peaks));
        
        sigma(k) = mean(dist2Valley);
    end
    
    mu = X(max_locs);
    dx = abs(diff(X([2 1])));

%     if any(data <= 0),
        [ model_pdf, ~ ] = fitGMMtoData(data, X, mu, sigma);
%     else
%         [ model_pdf, ~ ] = fitLogGMMtoData(data, X, peaks);
%     end

    model_pdf = model_pdf .* dx .* sum(histogram);

    [~, max_pdf_locs] = findLocalPeaks(model_pdf);

    bimodalFitTest = length(max_pdf_locs) < 2;
   
end

% Checks for Significant Peak/Valley Pairs, Passes if Finds None
function [bimodalPeaksTest, locs, nvd] = findSignificantPeaks(min_locs, max_locs, data)
    global DEPTH_OF_VALLEY_PERCENT
    global MIN_PERCENT_OF_HIST_PER_BIN
    
    shift_data = (data - MIN_PERCENT_OF_HIST_PER_BIN);
    norm_data = shift_data ./ max(shift_data);
    
    significant_pk_height = false(0);
    norm_valley_depth = [];
    for k = 1:length(min_locs),
        max_pk_L_ind = find(max_locs < min_locs(k), 1, 'last');
        max_pk_R_ind = find(min_locs(k) < max_locs, 1, 'first');

        if isempty(max_pk_L_ind) | isempty(max_pk_R_ind),
            norm_valley_depth(k) = 0;
        else
            max_pk_L = norm_data(max_locs(max_pk_L_ind));
            max_pk_R = norm_data(max_locs(max_pk_R_ind));
            min_pk   = norm_data(min_locs(k));

            norm_valley_depth(k) = min([max_pk_L max_pk_R]) - min_pk;
        end

        significant_pk_height(k) = norm_valley_depth(k) >= DEPTH_OF_VALLEY_PERCENT; 
        
    end
    
    locs = min_locs(significant_pk_height);
    nvd = norm_valley_depth;
    bimodalPeaksTest = all(~significant_pk_height);
    

    %% Attempt 1  
%     global SIGNIFICANT_AREA
%     global SIGNIFICANT_NORM_HEIGHT
% 
%     gmm_area_norm = sum(model_components, 2);
% 
%     significant_pk_area = 0;
%     significant_pk_height = 0;
%     max_pks = smoothedFitHistogram(max_locs);
%     min_pks = smoothedFitHistogram(min_locs);
%     
%     for n = 2:length(max_pks),
%         valley_found = max_locs(n-1) < min_locs & min_locs < max_locs(n);
% 
%         if any(valley_found),
%             peak_height = max_pks(n+[-1 0]) + min_pks(valley_found);
%             norm_peak_height = (peak_height - 0)/(max(norm_gmm_pdf_fit) - 0);
% 
%             significant_pk_area(n-1) = all(gmm_area_norm > SIGNIFICANT_AREA);
%             significant_pk_height(n-1) = all(norm_peak_height > SIGNIFICANT_NORM_HEIGHT);
%         end
%     end
%         
%     bimodalPeaks = any(significant_pk_height) | any(significant_pk_area);

end

% Checks for threshold crossings
function bimodalThreshold = checkThresholdCrossings(histogram)
    global MIN_PERCENT_OF_HIST_PER_BIN
    global MIN_NUM_CONSECUTIVE_BINS_THRESHOLD
    
    % Normalize Histogram
    if sum(histogram) > 1,
        histogram = histogram ./ sum(histogram);
    end
    
    % Min_percent crossing
    threshold_crossing = histogram > MIN_PERCENT_OF_HIST_PER_BIN;
    groups = bwlabel(~threshold_crossing);
    table = tabulate(groups);
    group_crossing = table(2:end,2) <= MIN_NUM_CONSECUTIVE_BINS_THRESHOLD;

    bimodalThreshold = all(group_crossing);
end


%% Helper Functions
function [min_locs, max_locs] = findLocalPeaks(data)

    if size(data,1) > size(data,2),
        data = data';
    end
    
    % Find all local maximum
    N = length(data);
    [~, max_locs] = findpeaks(data);
    [~, min_locs] = findpeaks(-data);    

    if ~isempty(max_locs) | ~isempty(min_locs)
        
        all_locs = sort([max_locs, min_locs]);
        
        if data(1) < data(all_locs(1)),
            min_locs = [1 min_locs];
        else
            max_locs = [1 max_locs];
        end

        if data(N) < data(all_locs(end)),
            min_locs = [min_locs N];
        else
            max_locs = [max_locs N];
        end        
    
    else
        if data(1) < data(N),
            min_locs = 1;
            max_locs = N;
        else
            min_locs = N;
            max_locs = 1;
        end
        
    end
    
end

function [min_locs, max_locs] = removeInsigPkValleyPairs(data, min_locs, max_locs)
    global INSIG_PK_VALLEY_PAIR_DIST
    global MIN_PERCENT_OF_HIST_PER_BIN
           
    % Normalize Data
    shift_data = (data - MIN_PERCENT_OF_HIST_PER_BIN);
    norm_data = shift_data ./ max(shift_data);
    N = length(norm_data);

    all_locs = sort([max_locs, min_locs]);
    is_max_loc = ismember(all_locs,max_locs);
            
    locs_to_remove = [nan];
    while ~isempty(locs_to_remove),

        % Check if any PK/Valley Pairs are insignificant
        norm_pk_val = norm_data(all_locs);
        diff_norm_pk_val = abs(diff(norm_pk_val));

        % Remove first Insignificant Pair
        pk_ind_to_remove = diff_norm_pk_val < INSIG_PK_VALLEY_PAIR_DIST;            

        all_locs_ind = find(pk_ind_to_remove, 1, 'first');
        locs_to_remove = all_locs([all_locs_ind+0 all_locs_ind+1]);
       
        all_locs(ismember( all_locs, locs_to_remove)) = [];
   
    end
    
    sig_locs = [ all_locs ];
    all_locs = sort([max_locs, min_locs]);
    is_max_loc = [is_max_loc(2) is_max_loc(ismember(all_locs,sig_locs)) is_max_loc(end-1)];
    sig_locs = [1 sig_locs N];
    
    remove_dups = ([ diff(sig_locs([1:end-1])) == 0, false, diff(sig_locs([end-1:end])) == 0 ]);
    sig_locs(remove_dups) = [];
    is_max_loc(remove_dups) = [];
    
    max_locs = sig_locs(is_max_loc);
    min_locs = sig_locs(~is_max_loc);
end

function [histogram, X, fitHistogram, fitX, fitData] = binData(data)
    global MIN_PERCENT_OF_HIST_PER_BIN

    % Bin Histogram
    length_X = round(sqrt(length(data)));  % 100;
    X = linspace(min(data), max(data), length_X);     
    histogram = histc(data, X)';
    
    % Normalize Histogram
    norm_histogram = histogram./sum(histogram);

    % Find Outliers 
    % Discard edge values of bins below min. percent (stop when two
    % consecutive bins are above min. percent)
    invaluable_points = norm_histogram > MIN_PERCENT_OF_HIST_PER_BIN;
    start_ind = find(invaluable_points == 1);
    start_ind = start_ind(find(diff(start_ind) == 1, 1, 'first'));
    end_ind = find(invaluable_points == 1);
    end_ind = end_ind(find(diff(end_ind) == 1, 1, 'last')+1);
    
    invaluable_points(:) = false;
    invaluable_points(start_ind:end_ind) = true;
    
    % If Binning did a bad job (less than half the points where
    % significant
    if mean(invaluable_points) < 0.5 & length(invaluable_points) > 3,
        
        % Rebin Histogram
        valuable_length = abs(start_ind-end_ind);
        if isempty(valuable_length),
            valuable_length = 2;
        end  
        
        new_interval = (max(data)-min(data)) / (length_X * (length_X/valuable_length));
        min_interval_of_data = abs(min(diff(unique(data))));
        new_X = [min(data):max([min_interval_of_data new_interval]):max(data)];
        new_histogram = histc(data, new_X)';

        % Normalize Histogram
        new_norm_histogram = new_histogram./sum(new_histogram);

        % Discard edge values of bins below min. percent
        new_invaluable_points = new_norm_histogram > MIN_PERCENT_OF_HIST_PER_BIN;
        new_start_ind = find(new_invaluable_points == 1);
        new_start_ind = new_start_ind(find(diff(new_start_ind) == 1, 1, 'first'));
        new_end_ind = find(new_invaluable_points == 1);
        new_end_ind = new_end_ind(find(diff(new_end_ind) == 1, 1, 'last')+1);

        new_invaluable_points(:) = false;
        new_invaluable_points(new_start_ind:new_end_ind) = true;
        
        if sum(new_invaluable_points) > sum(invaluable_points)
            invaluable_points = new_invaluable_points;
            histogram = new_histogram;
            X = new_X;
            start_ind = new_start_ind;
            end_ind = new_end_ind;
        end
    end   
    
    % If some outlier is preventing proper binning
    if all(invaluable_points == 0),
        med_sub_data = data - median(data);
        rec_med_sub_data = abs(med_sub_data);
        data( rec_med_sub_data > 6*std(data)) = [];
        
        [histogram, X, fitHistogram, fitX, fitData] = binData(data);
        
    else
        % Discard Outliers
        fitHistogram = histogram(invaluable_points);
        fitX = X(invaluable_points);
        fitData = data( X(start_ind) <= data & data <= X(end_ind) );
    end
    
end
    
function smoothedFitHistogram = smoothHistogram(data)

    N = length(data);

    % 3-Sample-Median Smoothing
    for i = 2:N-1,
        smoothedData(i) = median(data(i-1:i+1));
    end
    smoothedData([1 N]) = data([1 N]);
    data = smoothedData;            


    % Sliding Window Smoothing
    smoothedData = smooth(data,5);
    data = smoothedData';
    
    smoothedFitHistogram = data;
end

function [gmm_pdf, gmm_components] = fitGMMtoData(data, X, mu, sigma)

    warning('off', 'all')        
    
    dx = abs(diff(X([2 1])));
    
    if nargin < 3,
        S.mu = [X(1), X(end)]';
        K = 2;
        S.Sigma = ones(1,1,K) .* dx;
    else
        S.mu = mu(:)'; 
        K = length(S.mu);
        S.Sigma = ones(1,1,K) .* 2 * dx; %reshape(sigma,1,1,K);
    end

    if K > 1, 
        % GMM model
        S.PComponents = ones([K 1]) .* 1/K;
        options = statset('TolFun', 1e-2); 

        try
            gmm = gmdistribution.fit(data, K, ...
                                     'regularize', abs(min(diff(data)))/1000, ...
                                     'covtype', 'diagonal', ...
                                     'Start', S, ...
                                     'Options', options);
        catch
            gmm = gmdistribution.fit(data, K, ...
                                     'regularize', abs(min(diff(data)))/1000, ...
                                     'covtype', 'diagonal', ...
                                     'Replicates', 10, ...
                                     'Options', options);
        end
                                     
        model_components = [];
        for k = 1:K,
            gmm_piecewise_fit = normpdf( X', gmm.mu(k), sqrt(gmm.Sigma(:,:,k)));
            model_components(k,:) = gmm_piecewise_fit .* gmm.PComponents(k);
        end
        
        [~, sort_ind] = sort(gmm.mu);
        model_components = model_components(sort_ind, :);

    else
        model_components = normpdf(X, mean(data), std(data));

    end                   
                         
    gmm_components = model_components;
    gmm_pdf = sum(gmm_components, 1);
    
    warning('on', 'all') 
end


%% LOG GMM
% function [log_gmm_pdf, log_gmm_components] = fitLogGMMtoData(data, X, mu)
% 
%     warning('off', 'all')
% 
%     log_X = log(X);
%     log_data = log(data);
%     dx = abs(diff(X([2 1])));
%     
%     if nargin < 3,
%         S.mu = log_X([1 end])';
%         K = 2;
%     else
%         if size(mu,1) < size(mu,2), S.mu = log(mu)'; end
%         if size(mu,1) >= size(mu,2), S.mu = log(mu); end
%         K = length(S.mu);
%     end
%     
%     S.Sigma = ones(1,1,K) .* 2*dx;
%     S.PComponents = ones([K 1]).*1/K;
%     options = statset('TolFun', 1e-2); % 'MaxIter',10, 
%     
%     if K > 1 & length(log_data) > 10, 
%         % LOG space
%         gmm = gmdistribution.fit(log_data, K, ...
%                                  'regularize', abs(min(diff(log_data)))/1000, ...
%                                  'covtype', 'diagonal', ...
%                                  'Start', S, ...
%                                  'Options', options);
%         model_components = [];
%         for k = 1:K,
%             gmmLOG = makedist('Lognormal', 'mu', gmm.mu(k), ...
%                                            'sigma', sqrt(gmm.Sigma(:,:,k)));
%             gmm_piecewise_fit = gmm.PComponents(k) * pdf(gmmLOG, X');
%             model_components(k,:) = gmm_piecewise_fit .* dx;
%         end
% 
%         [~, sort_ind] = sort(gmm.mu);
%         model_components = model_components(sort_ind, :);
% 
%     else
%         model_components = normpdf(X, mean(data), std(data));
% 
%     end                   
% 
%     log_gmm_components = model_components;
%     log_gmm_pdf = sum(log_gmm_components, 1);
% 
%     warning('on', 'all') 
% end
