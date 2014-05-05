
function this = generateIterativeOutput(this)
           
    import hsst.scoreMethod.*
    import stdlib.*


    %% Setup Variables
    DEBUG_PLOTS = false;

    CONTINUE_SORT_FLAG = true;
    all_sortIDs = zeros(size(this.ts));

    UNSORTED_UNITS = 0;
    SORTED_UNITS = [];
    NOISE_UNITS = [];

    K = 1;
    diff_btw_prev_and_updated = [];

    %% Iterative
    while CONTINUE_SORT_FLAG,

        % Previous SortCode 
        prev_all_sortIDs = all_sortIDs;

        % Sort Any Unsorted Waveforms
        unsortedWf_ind = ismember(all_sortIDs, UNSORTED_UNITS);
        [recent_sortIDs] = this.SoMObj.sort(this.wf(:,unsortedWf_ind), ...
                                            this.ts(unsortedWf_ind), ...
                                            this.fs, ...
                                            this.noiseEstimate, ...
                                            [], [], [], ...
                                            K);                                        

        % Reintegrate recently sorted Waveforms into final_sortCode
        nextID = max(unique(all_sortIDs))+1 - min(unique(recent_sortIDs));
        all_sortIDs(unsortedWf_ind) = nextID + recent_sortIDs;


        % Score ReIntegrated Waveforms/sortIDs
        sortClassName = class(this.SoMObj);
        sortMethodLabel = eval(sprintf('%s.sortMethodLabel',sortClassName) );

        scoreClassName = class(this.ScMObj);
        sortLabel = sprintf('%s: %s', sortMethodLabel, 'Iterative');
%                 sortScoreObj = feval(scoreClassName,all_sortIDs, ...
%                                                     this.wf, ...
%                                                     this.ts, ...
%                                                     this.noiseEstimate, ...
%                                                     'string_label', sortLabel, ...
%                                                     'normalize_num_wf', true);
        sortScoreObj = sortQualityClass(all_sortIDs, ...
                                        this.wf, ...
                                        this.ts, ...
                                        this.noiseEstimate, ...
                                        'string_label', sortLabel, ...
                                        'normalize_num_wf', true);

        % Get Info from sortScoreObj
        unitIDs         = sortScoreObj.unitIDs;               
        is_resortedUnit = unitIDs > nextID;

        is_noise        = sortScoreObj.noise_units;
        is_undersorted  = sortScoreObj.undersorted_units;
        is_oversorted   = sortScoreObj.oversorted_units;                
        is_goodUnit     = ~is_noise & ~is_undersorted & ~is_oversorted;

        oversorted_ind  = [sortScoreObj.raw_values{5}; sortScoreObj.raw_values{6}];

        if DEBUG_PLOTS, plotPCA_WF(all_sortIDs, this.wf); end
%                 if DEBUG_PLOTS, plotPCA_WF(recent_sortIDs, this.wf(:,unsortedWf_ind)); end                


        %% Sort "Logic" for resorted Units

%                 % Set Noise
%                 NOISE_UNITS = [NOISE_UNITS, unitIDs(is_resortedUnit & is_noise)];
% 
%                 % Set "Good" Units
%                 SORTED_UNITS = [SORTED_UNITS unitIDs(is_resortedUnit & is_goodUnit)];

        % UnderSorted Units
        UNDERSORTED_UNITS = unitIDs(is_resortedUnit & is_undersorted);

        % OverSorted Units
        oversorted_unit_indices = find(is_resortedUnit & is_oversorted & ~is_undersorted);
        OVERSORTED_UNITS = unitIDs(oversorted_unit_indices);

        % Resort These Units on the next iteration
        UNSORTED_UNITS = [OVERSORTED_UNITS, UNDERSORTED_UNITS];


        while ~isempty(OVERSORTED_UNITS),
            OS_unitID           = OVERSORTED_UNITS(1);
            OS_with_unitIDs     = unitIDs( [oversorted_ind{:, oversorted_unit_indices}] );

            % For all units overlaping with this "oversorted" unit
            for OS_with_unitID = OS_with_unitIDs,                        
                sortID_OS_with_unitIndex = ismember(all_sortIDs, OS_with_unitID);
                all_sortIDs(sortID_OS_with_unitIndex) = OS_unitID;
                OVERSORTED_UNITS(OVERSORTED_UNITS == OS_with_unitID) = [];
            end
            OVERSORTED_UNITS(OVERSORTED_UNITS == OS_unitID) = [];

            UNSORTED_UNITS = [UNSORTED_UNITS OS_unitID];
        end

        UNSORTED_UNITS = unique(UNSORTED_UNITS);

        if DEBUG_PLOTS, plotPCA_WF(all_sortIDs, this.wf); end


        %% Continue Sorting?  
        diff_btw_prev_and_updated(end+1) = calculateACC(prev_all_sortIDs, all_sortIDs)

        END_REACHED = false;
        if length(diff_btw_prev_and_updated) > 2,
            END_REACHED = all(diff_btw_prev_and_updated(end-1:end) == 0);
        end
        CONTINUE_SORT_FLAG = ~(isempty(UNSORTED_UNITS) | END_REACHED);

        if CONTINUE_SORT_FLAG
            if K > length(UNSORTED_UNITS),
                if length(UNSORTED_UNITS) == 1,
                    K = 2;
                else
                    K = 1;
                end
            else
                K = K + 1;
            end
        end

    end


    %% Record Info
    all_sortIDs(ismember(all_sortIDs, NOISE_UNITS)) = 0;
    sortCode = all_sortIDs;

    sortLabel = sprintf('%s: %s', sortMethodLabel, 'Iterative');
    scoreObj = feval(scoreClassName,sortCode, ...
                                    this.wf, ...
                                    this.ts, ...
                                    this.noiseEstimate, ...
                                    'string_label', sortLabel, ...
                                    'normalize_num_wf', true);

    this.resultObj.parameter  = [1];                        
    this.resultObj.sortCode   = sortCode;
    this.resultObj.scoreObj   = scoreObj;
    this.resultObj.score      = scoreObj.getScore();

end