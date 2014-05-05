function block = UltraMegaSort2000_wrapper(catName, trial, channels, varargin)
% DUKESPIKESORTING_WRAPPER Wrapper for duke spike sorting code that applies metrics to determine best sort
%
% INPUTS
% =========================================================================
%   catName  - (char) Cat to sort
%   trial    - (numeric) block to sort
%   channels - (numeric) channels to sort
%
% VARARGIN
% =========================================================================
%   sortparameters         - (numeric) Default: [10^6 10^5 10^4 10^3 10^2]. Additive Noise Parameters
%   sortparameterslabels   - (cell)    Default: {'10^6', '10^5', '10^4', '10^3', '10^2'};
%   spikeDataObj           - (numeric) Default: 1. SpikeData Object Number (TankSort)
%   force_replace          - (logical) Default: false.  Overwrites Exisiting Data in either DB or generated data
%   force_score            - (logical) Default: false.
%   save_to_generated_data - (logical) Default: true. Save All SortCodes and SortScores to GENERATED_DATA *.mat file
%   save_to_qnap           - (logical) Default: false. Save sortCode to database
%   verbose                - (logical) Default: true. Display Print Statements to MATLAB console
%   showplots              - (logical) Default: false.
%   use_operators          - (logical) Default: false. Whether spike opera
%   unsorted_code          - (numeric) Default: 0. Value of unsorted data
%   sort_name              - (char)    Default: 'DukeSort'. Name of sort in database

% Input Variables

% TODO: varargin pass through to dukeSpikeSorting(wf, varargin)

varargin = sanitizeVarargin(varargin);
DEFINE_CONSTANTS
% Default Additive Noise Parameters
sortparameters          = 5e-6*[1:10 50 100]*1e-4;
sortparameterslabels    = [];
spikeDataObj = 1;
force_replace           = false;
force_score             = false;
save_to_generated_data  = true;
save_to_qnap            = false;
verbose                 = false;
showplots               = false;
use_operators           = false;
unsorted_code           = 0;
sort_name               = 'UltraMegaSort2000';
END_DEFINE_CONSTANTS
if isempty(sortparameterslabels)
    sortparameterslabels = arrayfun(@num2str,sortparameters,'un',0);
end
name = 2;
code = 1;

% Load Data from QNAP

if verbose,
    fprintf('Loading Cat: ''%s''.  Trial #%d.  Channel(s): %s \r', ...
        catName, trial, num2str(channels));
end

% Initialize Varibles

% Get main trial object
tr = getDBObj('trial', catName, trial);

% Number of Channels in that trial
spike_data      = tr.spikeData(spikeDataObj);
db_channels     = spike_data.st.channel;
db_channels     = [db_channels{:}];
number_channels = length(db_channels);


% Throw Error if Any Channels Are Invalid
mask = ismember(channels,db_channels);

assert(all(mask),'Invalid Channel Number(s): [%s]',array2strMatlabStyle(channels(~mask)));

% Establish File Name and GENERATE DATA

[C ~] = setupConvPathForCat(catName);
createFolderIfNoExist(C.ANALYSIS_PATHS.DUKESORT_ROOT);

blockName = sprintf('block-%d', trial);
filePath = fullfile(fileparts(C.ANALYSIS_PATHS.DUKESORT_ROOT),'UMS',blockName);
createFolderIfNoExist(filePath);

% Sorting and Saving
for this_chan = channels,
    fprintf('%02d/%02d\n',this_chan,length(channels));
    save_flag = false;	% boolean flag indicating save is necessary
    fullFilePath = fullfile(filePath, [blockName, sprintf('-%d.mat', this_chan)]);
    % File exists.  Load it.
    if exist(fullFilePath, 'file') == 2,
        load(fullFilePath);
        if verbose, fprintf('Loaded ''%s'' file. \r', fullFilePath); end
        block.updated{spikeDataObj} = false(1,length(block.sortScores{spikeDataObj}));
    else
        block = struct;
        block.sortCodes  = {};
        block.sortScores = {};
        block.sortParameterInDataBase = [];
        block.spike_mask = [];
        save_flag = true;
        %formattedWarning('MATLAB:FileNotFound', 'Creating: %s', fullFilePath);
    end
    if lengthprop(spike_data.st(this_chan),'ts') == 0
        continue;
    end
    
    % Sort Data
    
    % retrieve existingSortCodes
    try
        existingSortCodes = [block.sortCodes{spikeDataObj}{name, :}];
    catch MSerror
        if strcmp(MSerror.identifier,'MATLAB:nonExistentCellElement') || ...
                strcmp(MSerror.identifier,'MATLAB:badsubscript'),
            block.sortCodes{spikeDataObj} = {};
        end
        
        existingSortCodes = [block.sortCodes{spikeDataObj}{name, :}];
    end
    
    % Find previous sorts with matching names
    param_idx = [];   % array of indices where each sort is found, 0-if not found
    nmbr_exst_SortCodes = length(existingSortCodes);
    index_of_exst_SortCodes = [1:nmbr_exst_SortCodes];
    for iiParam = 1:length(sortparameterslabels),
        sortLabel = sortparameterslabels{iiParam};
        index_match = index_of_exst_SortCodes(strcmp(existingSortCodes, sortLabel));
        
        if isempty(index_match)
            index_match = 0;
        end
        
        param_idx(iiParam) = index_match;
    end
    
    % If any sorts are not found, or force_replace is true
    if any(param_idx == 0) || force_replace,
        
        % Saving to GENERATED_DATA is necessary due to changes to block
        save_flag = true;
        
        if isprop(tr,RawData.listAs)
            noise_est = tr.rawData.wf(this_chan).noiseEstimate;
        elseif isprop(tr,PDecData.listAs)
            noise_est = estimateNoise(tr.pdecData.wf(this_chan));
        else
            noise_est = tr.spikeData.st(this_chan).threshold;
        end
        Fs = spike_data.fs;
        
        spikes = ss_default_params(Fs);
        % UMS2000 Interface ===============================================
        if ~isprop(tr,RawData.listAs)
            rethresholded   = false;
            if use_operators
                [ts,spike_mask] = getSpikesUsingOperators(spike_data.st(this_chan));
                aligned_wf      = spike_data.st(this_chan).snips.wf;
                aligned_wf      = aligned_wf(:,spike_mask);
            else
                aligned_wf  = spike_data.st(this_chan).snips.wf;
                ts          = spike_data.st(this_chan).ts;
                spike_mask  = true(length(ts),1);
            end
            aligned_wf       = alignment_function(aligned_wf,11);
            if isempty(aligned_wf)
                % this can be empty if the operators filtered all spikes
                continue;
            end
            nSnips   = size(aligned_wf,2);
            nSamples = size(aligned_wf,1);
            
            spikes.waveforms  = aligned_wf';
            spikes.spiketimes = ts;
            spikes.trials     = ones(1,nSnips);
            % CAA TODO These aren't really used
            spikes.params.max_jitter   = 3/Fs*1e3;
            spikes.params.window_size  = (nSamples)/Fs*1e3;   % ms, width of a spike
            spikes.params.shadow       = 2/Fs*1e3;          % ms, enforced dead region after each spike
            spikes.params.cross_time   = 0/Fs*1e3;         % ms, alignment point for peak of waveform
            
            % spikes.info.detect.stds          = 28.8141
            spikes.info.detect.thresh        = spike_data.st(this_chan).threshold;
            spikes.info.detect.event_chann = ones(nSnips,1);
            spikes.info.detect.dur           = tr.duration;
            spikes.info.detect.align_sample  = spike_data.st(this_chan).snips.nSamplesPreThreshold+1;
            spikes.unwrapped_times           = ts;
            % these values are ordinarily generated by the
            % align function
            [pca.u,pca.s,pca.v] = svd(detrend(spikes.waveforms(:,:),'constant'), 0);             % SVD the data matrix
            spikes.info.pca = pca;
        else
            % turn operators off because the spike times are extracted from
            % the raw data
            use_operators = false;
            rethresholded = true;
            spikes      = ss_detect({tr.rawData.wf(this_chan).data},spikes);
            spikes      = ss_align(spikes);
            spike_mask  = true(length(spikes.spiketimes),1);
            
        end
        
        % disable progress bars to reduce runtime
        options.progress = false;
        % ss_kmeans stalls when there are too few spikes, in this case
        % assign all spikes to unsorted
        if size(spikes.waveforms,1) > 10
            spikes = ss_kmeans(spikes,options);
            spikes = ss_energy(spikes,options);
            
            % For each sort in the sortparameters list
            for iiParam = 1:length(sortparameterslabels),
                spikes_copy = spikes;
                % if the sort is already in GEN_DATA or force replace = true
                if param_idx(iiParam) == 0 || force_replace,
                    
                    % if the sort isn't in GEN_DATA, index it
                    if param_idx(iiParam) == 0,
                        nmbr_exst_SortCodes  = nmbr_exst_SortCodes +1;
                        param_idx(iiParam) = nmbr_exst_SortCodes;
                    end
                    
                    % Additive Noise Parameter
                    agg_param = sortparameters(iiParam);
                    label = sortparameterslabels(iiParam);
                    try
                        spikes_copy.params.agg_cutoff = agg_param;
                        spikes_copy = ss_aggregate(spikes_copy);
                        sort_ids    = spikes_copy.assigns;
                    catch ME
                        simpleExceptionDisplay(ME);
                        sort_ids = ones([1 length(spikes.spiketimes)]);
                        formattedWarning('!!!! SORTCODE with Param = %s has triggered this fault: ''%s'' !!!!', ...
                            label{1,1}, ME.message);
                    end
                    
                    % Sort Score generated
                    score_object = sortQualityClass(sort_ids, spikes.waveforms', spikes.spiketimes, ...
                        'max_num_units',    6,      ...
                        'template_sigma',   2,      ...
                        'string_label',     label,  ...
                        'verbose',          false,  ...
                        'mean_threshold',   noise_est);
                    
                    sort_ids = int8(sort_ids);
                    block.sortCodes{spikeDataObj}{code, param_idx(iiParam)} = sort_ids;
                    block.sortCodes{spikeDataObj}{name, param_idx(iiParam)} = label;
                    block.sortScores{spikeDataObj}(param_idx(iiParam))      = score_object;
                    block.spike_mask{spikeDataObj}                          = spike_mask;
                    block.updated{spikeDataObj}(param_idx(iiParam))         = true;
                    block.rethresholded                                    = rethresholded;
                    block.spikes{spikeDataObj}                              = spikes;
                    % if sort is found in GEN_DATA, and force replace = false
                elseif param_idx(iiParam) > 0 && ~force_replace,
                    % Do Nothing!  Sort is already finished
                else
                    % Error State
                    error('Problem');
                end
            end
        else
            % unclear whether this is necessary
            sort_ids = unsorted_code*ones(size(spikes.spiketimes));
            for iiParam = 1:length(sortparameterslabels)
                % if the sort isn't in GEN_DATA, index it
                if param_idx(iiParam) == 0,
                    nmbr_exst_SortCodes  = nmbr_exst_SortCodes +1;
                    param_idx(iiParam) = nmbr_exst_SortCodes;
                end
                label = sortparameterslabels(iiParam);
                score_object = sortQualityClass(sort_ids, spikes.waveforms', spikes.spiketimes, ...
                    'max_num_units',    6,      ...
                    'template_sigma',   2,      ...
                    'string_label',     label,  ...
                    'verbose',          false,  ...
                    'mean_threshold',   noise_est);
                
                sort_ids = int8(sort_ids);
                block.sortCodes{spikeDataObj}{code, param_idx(iiParam)} = sort_ids;
                block.sortCodes{spikeDataObj}{name, param_idx(iiParam)} = label;
                block.sortScores{spikeDataObj}(param_idx(iiParam))      = score_object;
                block.spike_mask{spikeDataObj}                          = spike_mask;
                block.updated{spikeDataObj}(param_idx(iiParam))         = true;
                block.spikes{spikeDataObj}                              = spikes;
                block.rethresholded                                     = rethresholded;
            end
        end
    end
    
    
    % Retreive Existing SortScore Versions
    try
        existingSortVers = [block.sortScores{spikeDataObj}.VERSION_NUMBER];
    catch MSerror
        if strcmp(MSerror.identifier,'MATLAB:nonExistentCellElement')
            block.sortScores{spikeDataObj} = {1};
        end
        existingSortVers = [block.sortScores{spikeDataObj}.VERSION_NUMBER];
    end
    
    % Boolean Array - out of data scores
    if ~force_score
        outOfDateScore = existingSortVers < sortQualityClass.MOST_CURRENT_VERSION_NUMBER;
    else
        outOfDateScore = ~block.updated{spikeDataObj};
    end
    % if any scores are out of date
    if any(outOfDateScore),
        
        if ~(any(param_idx == 0) || force_replace)
            spike_mask  = block.spike_mask{spikeDataObj};
            aligned_wf  = spike_data.st(this_chan).snips.wf(:,spike_mask);
            ts          = spike_data.st(this_chan).ts(spike_mask);
            aligned_wf  = alignment_function(aligned_wf,11);
            if isprop(tr,RawData.listAs)
                noise_est = tr.rawData.wf(this_chan).noiseEstimate;
            elseif isprop(tr,PDecData.listAs)
                noise_est = estimateNoise(tr.pdecData.wf(this_chan));
            else
                noise_est = tr.spikeData.st(this_chan).threshold;
            end
        end
        
        % Saving to GENERATED_DATA is necessary due to changes to block
        save_flag = true;
        
        % update each out of date score
        for iiParam = 1:length(outOfDateScore),
            if outOfDateScore(iiParam) || force_replace,
                
                sort_ids = block.sortCodes{spikeDataObj}{code, iiParam};
                label    = block.sortCodes{spikeDataObj}{name, iiParam};
                
                % Sort Score generated
                [score_object] = sortQualityClass(sort_ids, aligned_wf, ts, ...
                    'max_num_units',    6,      ...
                    'template_sigma',   2,      ...
                    'string_label',     label,  ...
                    'verbose',          false,  ...
                    'mean_threshold',   noise_est);
                
                block.sortScores{spikeDataObj}(iiParam) = score_object;
                
            end
        end
    end
    
    % Index of Highest Scoring sortcode
    existingSortVers = [block.sortScores{spikeDataObj}.score];
    [~, block.sortParameterInDataBase] = max(existingSortVers);
    
    % Save *.mat file to GENERATED_DATA folder
    if save_flag && save_to_generated_data,
        save(fullFilePath, 'block');
    end
    
end

% Save to QNAP
if save_to_qnap,
    
    sortCodes     = {};
    chann_numbers = [];
    num_channels  = length(channels);
    
    index = 1;
    if ~block.rethresholded
        for this_chan = channels,
            
            fullFilePath = fullfile(filePath, [blockName, sprintf('-%d.mat', this_chan)]);
            
            if exist(fullFilePath, 'file') == 2,
                load(fullFilePath);
                
                sortIndex = block.sortParameterInDataBase;
                % map duke sort codes to linearly increasing from 1
                [uCode,~,sort_ids]  = unique_2011b(block.sortCodes{spikeDataObj}{code, sortIndex});
                
                if use_operators
                    [ts,spike_mask]  = getSpikesUsingOperators(spike_data.st(this_chan));
                    sortCodes{index} = zeros(size(spike_mask));
                    sortCodes{index}(spike_mask) = sort_ids;
                else
                    sortCodes{index} = sort_ids;
                end
                chann_numbers(index) = this_chan;
                index = index + 1;
            end
        end
        
        if verbose, fprintf('Saving sortIDs to QNAP. \r'); end
        
        saveSortCodes(spike_data, ...
            sort_name, ...
            sortCodes, ...
            chann_numbers, ...
            'force_replace', force_replace);
    else
        mutex      = lockExperimentDatabase('Erol')
        % CAA TODO This should be a static member of SpikeDAta
        if length(tr.spikeData) == 1
            % verify that this sort does not exist
            mask           = strcmp(tr.spikeData.sortCodes,sort_name);
            assert(~any(mask) == 1,'Sort ''%s'' already exists in the database, cannot rethreshold under this name',sort_name)
            addobj(tr,'SpikeData');
            spike_data = tr.spikeData(end);
            setSortCodes(spike_data,{sort_name});
        else
            % object exists, check if old objects are being overwritten
            spike_data  = tr.spikeData;
            sort_mask   = arrayfun(@(x)any(strcmp(sort_name,x.sortCodes)),spike_data);
            assert(sum(sort_mask) == 1,'Sort ''%s'' appears in more than one object',sort_name)
            if any(sort_mask)
                % clear out the data from the existing sort
                %                 if ~force_replace
                %                     val = questdlg(sprintf('Overwrite existing sort for ''%s''',sort_name),'Overwrite?','yes','no','no');
                %                     if strcmp(val,'no')
                %                         return;
                %                     end
                %                 end
                spike_data = tr.spikeData(sort_mask);
                if isprop(spike_data,'st')
                    nChannels   = lengthprop(spike_data,'st');
                    for iiChannel = nChannels:-1:1
                        remobj(spike_data,'SpikeTimes',iiChannel);
                    end
                end
            else
                % this is a new sort, add it
                addobj(tr,'SpikeData');
                spike_data = tr.spikeData(end);
                setSortCodes(spike_data,{sortName});
            end
        end
        spike_data.fs = tr.rawData.fs;
        save(spike_data);
        %         nChannels = size(
        nChannels = length(channels);
        for iiChannel = 1:nChannels
            fprintf('%02d/%02d\n',iiChannel,nChannels);
            this_chan   = channels(iiChannel);
            mask        = ismember(tr.rawData.channels,this_chan);
            raw_channel = tr.rawData.wf(mask);
            fullFilePath = fullfile(filePath, [blockName, sprintf('-%d.mat', this_chan)]);
            
            if exist(fullFilePath, 'file') == 2,
                load(fullFilePath);
                
                % CREATE DB OBJECTS =======================================
                % SpikeTimes ==============================================
                addobj(spike_data,'SpikeTimes');
                spike_channel = spike_data.st(iiChannel);
                spikes = block.spikes{spikeDataObj};
                spike_channel.ts        = spikes.spiketimes;
                spike_channel.channel   = this_chan;
                spike_channel.signature = raw_channel.signature;
                spike_channel.threshold = spikes.info.detect.thresh;
                
                keepFlag = true;
                idx      = block.sortParameterInDataBase;
                initSorts(spike_channel,{sort_name},block.sortCodes{spikeDataObj}{idx}(:),true);
                
                spike_channel = addlink(spike_channel,'chanMap',getRecChannelMap(spike_channel));
                
                % Snippet =================================================
                addobj(spike_channel,'Snippet');
                snip    = spike_channel.snips;
                snip.wf = spikes.waveforms';
                
                %                 snip.nSamplesPreThreshold = spikes,;
                %                 snip.nSamplesTotal        = nan;
                %                 snip.nSamplesOmitted      = 0
                close(raw_channel);
            end
        end
        save(tr)
        delete(mutex)
    end
end
end % dukeSpikeSorting function