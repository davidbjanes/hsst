
function sortCode = MSorter_DavidMod(wf, thresh, sortparameter)

    %% set parameters
    % S scales: 1st S for template generation, 2nd S for detection
    Svalue = [7,3]; 
    % template number
    tempclass = sortparameter; 
    % tau
    Thr = [thresh thresh]; % detection threshold (V)
    mode = 'neg'; % threshold type, 'neg' for negative spikes, 'pos' for positive spikes, 'both' for both kinds

    % supplementary parameters
    Fr = 40.000; % sampling frequency (kHz)
    Thr_min = 0; % Minimum peak-to-peak value (V)
    Thr_max = inf; % Maximum peak value (V)
    corr_sig = .5; % Correlation coefficient significance level
    dist_sig = 10; % Distance to the template (V square)

    %% M-Sorter
    
    % Normalize Waveforms
    max_value = max(max(abs(wf)));
    wf = wf ./ max_value*2;
    
    if size(wf,1) < size(wf,2),
        spikes = wf';
    else
        spikes = wf;
    end
            
    num_samples = size(spikes,2);
    num_wf = size(spikes,1);
    
    index = [1 num_wf];
    duration = num_samples * ones([1 num_wf]);
    sortCode = ones([1, num_wf]);
    t_marker = 0;

    for S1=1:2
%         tic
        if S1==1
%             disp('Template Generating...');
            if size(spikes,1) > 10
                fname = '';
                template = MSorter_template_generation(fname,Fr,spikes,tempclass,Thr_max);
                t_marker = 1;
            else
%                 disp('poor recording quality: no qualified spikes')
                t_marker = 0;
            end
        elseif t_marker == 1
%             disp('Detecting...');
%             disp('Clustering...');
            fname = '';
            sortCode = MSorter_template_matching(fname,template,Fr,spikes,index,duration,Thr_min,dist_sig,corr_sig);
%             disp('Done!');
        end
%         toc

    end
    
    [sortCode_val, sortCode] = max(sortCode, [], 1);
    sortCode(sortCode_val == 0) = 0;
    
end