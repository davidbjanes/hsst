%% Duke Sorting Algorithm Wrapper for QNAP

function [sortCode] = dukeSpikeSorting_wrapper(catName, trial, channels, varargin)
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

    %% Input Variables
    
    % TODO: varargin pass through to dukeSpikeSorting(wf, varargin)
    
    varargin = sanitizeVarargin(varargin);
    DEFINE_CONSTANTS
        % Default Additive Noise Parameters
        sortparameters          = [10^6 10^5 10^4 10^3 10^2];
        sortparameterslabels    = {'10^6', '10^5', '10^4', '10^3', '10^2'};

        % Default SpikeData Object Number (TankSort)
        spikeDataObj = 1;              

        % Automatically Overwrites Exisiting Data in either QNAP or GENERATED_DATA file
        force_replace           = false;        

        force_score             = false;
        
        % Save All SortCodes and SortScores to GENERATED_DATA *.mat file
        save_to_generated_data  = true;

        % Save sortCode to database
        save_to_qnap            = false;

        % Display Print Statements to MATLAB console
        verbose                 = true;
        showplots               = false;
        
        use_operators           = false;
        
        unsorted_code           = 0;
        sort_name               = 'DukeSort';
    END_DEFINE_CONSTANTS
    
    name = 2;
    code = 1;

    %% Load Data from QNAP

    if verbose,
        fprintf('Loading Cat: ''%s''.  Trial #%d.  Channel(s): %s \r', ...
                catName, trial, num2str(channels));  
    end

    %% Initialize Varibles

    % Get main trial object
    tr = getDBObj('trial', catName, trial);

    % Number of Channels in that trial
    number_channels = size(tr.spikeData(spikeDataObj).st, 2);


    %% Throw Error if Any Channels Are Invalid

    if any(channels > number_channels) | any(channels < 1),
        error('Error: Invalid Channel Number')
    end

    %% Establish File Name and GENERATE DATA

    [C ~] = setupConvPathForCat(catName);
    createFolderIfNoExist(C.ANALYSIS_PATHS.DUKESORT_ROOT);

    blockName = sprintf('block-%d', trial);
    filePath = C.ANALYSIS_PATHS.DUKESORT_ROOT;
    filePath = [filePath, '\', blockName];

    if exist(filePath, 'dir') ~= 7,
        mkdir(filePath);
    end
    

    %% Sorting and Saving 
    for chann = channels,
        
        save_flag = false;	% boolean flag indicating save is necessary

        fullFilePath = [filePath, '\', blockName, sprintf('-%d.mat', chann)];
        
        % File exists.  Load it.
        if exist(fullFilePath, 'file') == 2,
            load(fullFilePath);
            if verbose, fprintf('Loaded ''%s'' file. \r', fullFilePath); end

        % File doesn't exist.  Create it.
        else    
            block = struct;
            block.sortCodes = {};
            block.sortScores = {};
            block.sortParameterInDataBase = [];
            
            save_flag = true;
            warning('MATLAB:FileNotFound', 'Creating: %s', fullFilePath);
        end
        
    
        %% Sort Data
        
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
        match_found = [];   % array of indices where each sort is found, 0-if not found
        nmbr_exst_SortCodes = length(existingSortCodes);
        index_of_exst_SortCodes = [1:nmbr_exst_SortCodes];
        for i = 1:length(sortparameterslabels),
            sortLabel = sortparameterslabels{i};
            index_match = index_of_exst_SortCodes(strcmp(existingSortCodes, sortLabel));
            
            if isempty(index_match)
                index_match = 0;
            end
            
            match_found(i) = index_match;
        end 
        
        % Check if user wants to override existing sorts
%         if any(matches) > 0 && ~force_replace,
%             force_replace = strcmp(questdlg('Replace Existing Sort Codes in GENERATED_DATA?', ...
%                         'Replace Sorts Warning','yes','no','no'), 'yes')
%         end
        
        % If any sorts are not found, or force_replace is true
        if any(~match_found) == true || force_replace,
            
            % Saving to GENERATED_DATA is necessary due to changes to block
            save_flag = true;
            
            if use_operators
                [ts,spike_mask] = getSpikesUsingOperators(tr.spikeData(spikeDataObj).st(chann));
                wf        = tr.spikeData(spikeDataObj).st(chann).snips.wf;
                wf        = wf(:,spike_mask);
            else
                wf = tr.spikeData(spikeDataObj).st(chann).snips.wf;
                ts = tr.spikeData(spikeDataObj).st(chann).ts;
            end
            
            if isprop(tr,RawData.listAs)
                noise_est = tr.rawData.wf(chann).noiseEstimate;
            elseif isprop(tr,PDecData.listAs)
                noise_est = estimateNoise(tr.pdecData.wf(chann));
            else
                noise_est = tr.spikeData.st(chann).threshold;
            end
            
            % For each sort in the sortparameters list
            for i = 1:length(sortparameterslabels),

                % if the sort is already in GEN_DATA or force replace = true
                if match_found(i) == 0 || force_replace,

                    % if the sort isn't in GEN_DATA, index it 
                    if match_found(i) == 0,
                        nmbr_exst_SortCodes = nmbr_exst_SortCodes +1;
                        match_found(i) = nmbr_exst_SortCodes;                    
                    end

                     % Additive Noise Parameter
                    ANP = sortparameters(i);
                    label = sortparameterslabels(i);
                    
                    try 
                        % Sort Code Created 
                        [sortCode,aligned_wf] = dukeSpikeSorting(wf, ...
                                                    'ANP', ANP, ...
                                                    'verbose', verbose, ...
                                                    'showplots', showplots);
                    catch error_message
                        simpleExceptionDisplay(error_message);
                        sortCode = ones([1 length(wf)]);
                        warning('!!!! SORTCODE with Param = %s has triggered this fault: ''%s'' !!!!', ...
                                        label{1,1}, error_message.message);
                    end

                    % Sort Score generated
                    [scoreObject] = sortQualityClass(sortCode, aligned_wf, ts, ...
                                                    'string_label', label, ...
                                                    'verbose',      false, ...
                                                    'mean_threshold', noise_est);                        


                    sortCode = int8(sortCode);                      
                    block.sortCodes{spikeDataObj}{code, match_found(i)} = sortCode;
                    block.sortCodes{spikeDataObj}{name, match_found(i)} = label;
                    block.sortScores{spikeDataObj}(match_found(i))      = scoreObject;

                % if sort is found in GEN_DATA, and force replace = false
                elseif match_found(i) > 0 && ~force_replace,

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
            outOfDateScore = true(1,length(existingSortVers));
        end
        % if any scores are out of date        
        if any(outOfDateScore), 
            
            if ~(any(~match_found) == true || force_replace),
                    
                if use_operators
                    [ts,spike_mask] = getSpikesUsingOperators(tr.spikeData(spikeDataObj).st(chann));
                    wf        = tr.spikeData(spikeDataObj).st(chann).snips.wf;
                    wf        = wf(:,spike_mask);
                else
                    wf = tr.spikeData(spikeDataObj).st(chann).snips.wf;
                    ts = tr.spikeData(spikeDataObj).st(chann).ts;
                end
                wf = alignment_function(wf,11);
                if isprop(tr,RawData.listAs)
                    noise_est = tr.rawData.wf(chann).noiseEstimate;
                elseif isprop(tr,PDecData.listAs)
                    noise_est = estimateNoise(tr.pdecData.wf(chann));
                else
                    noise_est = tr.spikeData.st(chann).threshold;
                end
            end

            % Saving to GENERATED_DATA is necessary due to changes to block           
            save_flag = true;
            
            % update each out of date score
            for i = 1:length(outOfDateScore),
                if outOfDateScore(i) == 1 || force_replace,

                    sortCode = block.sortCodes{spikeDataObj}{code, i};
                    label    = block.sortCodes{spikeDataObj}{name, i};
                    
                    % Sort Score generated
                    [scoreObject] = sortQualityClass(sortCode, wf, ts, ...
                                                    'string_label',     label, ...
                                                    'verbose',          false, ...
                                                    'mean_threshold',   noise_est);                        
                   
                    block.sortScores{spikeDataObj}(i) = scoreObject;

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
        
        sortCodes = {};
        channel_numbers = [];
        num_channels = length(channels);
        
        index = 1;
        for chann = channels,
            
            fullFilePath = [filePath, '\', blockName, sprintf('-%d.mat', chann)];
            
            if exist(fullFilePath, 'file') == 2,
                load(fullFilePath);
        
                sortIndex = block.sortParameterInDataBase;
                % map duke sort codes to linearly increasing from 1
                [uCode,~,sortCode]  = unique_2011b(block.sortCodes{spikeDataObj}{code, sortIndex});
                
                if use_operators
                    [ts,spike_mask]  = getSpikesUsingOperators(tr.spikeData(spikeDataObj).st(chann));
                    sortCodes{index} = zeros(size(spike_mask));
                    sortCodes{index}(spike_mask) = sortCode;
                else
                    sortCodes{index} = sortCode;
                end
                channel_numbers(index) = chann;
                index = index + 1;
            end
        end
        
        if verbose, fprintf('Saving sortIDs to QNAP. \r'); end
        
        saveSortCodes(tr.spikeData(spikeDataObj), ...
                      sort_name, ...
                      sortCodes, ...
                      channel_numbers, ...
                      'force_replace', force_replace);
    end
    

    sortCode = block;    
    
end % dukeSpikeSorting function


