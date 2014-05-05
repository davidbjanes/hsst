classdef interfacer < handle
    
    properties
        
    end
    
    methods
        
        function output = sort( this, cat_name, block_num, channels, varargin)
        % VARARGIN HANDLER --------------------------------------------------------

        varargin = sanitizeVarargin(varargin);

        DEFINE_CONSTANTS
        sortparameters          = 5e-6*[1:10 50 100]*1e-4;  % Default Additive Noise Parameters
        sortparameterslabels    = [];
        spikeDataObj            = 1;        % Default SpikeData Object Number (TankSort)
        force_replace           = false;    % Automatically Overwrites Exisiting Data in either QNAP or GENERATED_DATA file
        force_score             = false;
        save_to_generated_data  = true;     % Save All SortCodes and SortScores to GENERATED_DATA *.mat file
        save_to_qnap            = false;    % Save sortCode to database
        verbose                 = false;    % Display Print Statements to MATLAB console
        showplots               = false;
        use_operators           = false;
        unsorted_code           = 0;
        END_DEFINE_CONSTANTS


        % Clean Up Varargins ------------------------------------------------------
            if isempty(sortparameterslabels)
                sortparameterslabels = arrayfun(@num2str,sortparameters,'un',0);
            end
            name = 2;
            code = 1;

        % Initialize Varibles -----------------------------------------------------
            % Get main trial object
            tr         = getDBObj('trial', cat_name, block_num);
            cat_name   = cat_name;
            block_num  = block_num;
            channels   = channels;

            % Number of Channels in that trial
            spike_data      = tr.spikeData(spikeDataObj);
            db_channels     = spike_data.st.channel;
            db_channels     = [db_channels{:}];
            number_channels = length(db_channels);


        % Throw Error if Any Channels Are Invalid ---------------------------------
            mask = ismember(channels,db_channels);

            assert(all(mask),'Invalid Channel Number(s): [%s]', array2strMatlabStyle(channels(~mask)));


        % Establish File Name and GENERATE DATA -----------------------------------
            [C ~] = setupConvPathForCat(cat_name);
            blockName = sprintf('block-%d', block_num);
            GENERATED_DATA_DIR = fullfile(C.ANALYSIS_PATHS.ROOT, ...
                                          this.sortMethodLabel, ...
                                          blockName);
            createFolderIfNoExist(GENERATED_DATA_DIR);      


        % Set up for channel sort (channel by channel) ----------------------------

            for this_chan = channels,

                % Verbose Output Indicating Processing of Channel Num
                fprintf('%02d/%02d\n',this_chan,length(channels));

                save_flag = false;	% boolean flag indicating save is necessary
                fullFilePath = fullfile(GENERATED_DATA_DIR, ...
                                        [blockName, sprintf('-%d.mat', this_chan)]);

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

        % Sort Data --------------------------------------------------------------- 

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

                    if use_operators
                        [ts,spike_mask] = getSpikesUsingOperators(spike_data.st(this_chan));
                        wf              = spike_data.st(this_chan).snips.wf;
                        wf              = wf(:,spike_mask);
                    else
                        wf              = spike_data.st(this_chan).snips.wf;
                        ts              = spike_data.st(this_chan).ts;
                        spike_mask      = true(length(ts),1);
                    end            

                    thresh        = spike_data.st(this_chan).threshold;
                    align_sample  = spike_data.st(this_chan).snips.nSamplesPreThreshold+1;
                    dur           = tr.duration;

                    if ~isprop(tr,RawData.listAs)
                        raw_data        = [];
                        rethresholded   = false;

                    else
                        raw_data        = {tr.rawData.wf(this_chan).data};
                        use_operators   = false;
                        rethresholded   = true;

                        spike_mask  = true(length(ts),1);
                    end

                    % For each sort in the sortparameters list
                    for iiParam = 1:length(sortparameterslabels),

                        % if the sort is already in GEN_DATA or force replace = true
                        if param_idx(iiParam) == 0 || force_replace,

                            % if the sort isn't in GEN_DATA, index it
                            if param_idx(iiParam) == 0,
                                nmbr_exst_SortCodes  = nmbr_exst_SortCodes +1;
                                param_idx(iiParam) = nmbr_exst_SortCodes;
                            end

                            agg_param = sortparameters(iiParam);
                            label = sortparameterslabels(iiParam);
                            [sort_ids, wf, ts, property] = this.sortMethod(wf, ...
                                                             ts, Fs, thresh, ...
                                                             align_sample, dur, ...
                                                             raw_data, agg_param);

                            if length(sort_ids) < 10,
                                sort_ids = unsorted_code * ones(size(sort_ids));
                            end

                            % Sort Score generated
                            score_object = sortQualityClass(sort_ids, wf, ts, ...
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
                            block.rethresholded                                     = rethresholded;
                            block.property{spikeDataObj}                            = property;


                        % if sort is found in GEN_DATA, and force replace = false
                        elseif param_idx(iiParam) > 0 && ~force_replace,
                            % Do Nothing!  Sort is already finished
                        else
                            % Error State
                            error('Problem');
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

        end
    end
    
end

