
function this = generateIterativeOutput_v2(this)
           
    import stdlib.*

    %% Setup Variables
    DEBUG_PLOTS = true;

    CONTINUE_SORT_FLAG = true;
    all_sortIDs = zeros(size(this.ts));

    UNSORTED_UNITS = 0;
    SORTED_UNITS = [];
    NOISE_UNITS = [];

    K = 1;
    diff_btw_prev_and_updated = [];

    %% Iterative
    while CONTINUE_SORT_FLAG,

        % Save Previous SortCode 
        prev_all_sortIDs = all_sortIDs;
        
        
        %% Sort Unsorted DATA
        [resorted_sortCode, nextID] = sortData(this, all_sortIDs, UNSORTED_UNITS, K);


        %% Score Unsorted DATA
        [unitIDs, ...
         UNDERSORTED_UNITS, ...
         OVERSORTED_UNITS, ...
         OS_pairs] = scoreUpdatedSort(this, resorted_sortCode, nextID);
        
        
        %% PLot?
%         if DEBUG_PLOTS, plotPCA_WF(prev_all_sortIDs, this.wf); end
        if DEBUG_PLOTS, plotPCA_WF(resorted_sortCode, this.wf); end
%         if DEBUG_PLOTS, plotPCA_WF(recent_sortIDs, this.wf(:,unsortedWf_ind)); end  


        %% Merge/Mark Units for resorting based on Score
        rearranged_sortCode = resorted_sortCode;
        
        UNSORTED_UNITS = [OVERSORTED_UNITS, UNDERSORTED_UNITS];
        
        while ~isempty(OVERSORTED_UNITS),
            OS_unitID           = OVERSORTED_UNITS(1);
            OS_unitIndex        = unitIDs == OS_unitID;
            OS_with_unitIDs     = OS_pairs{OS_unitIndex}(:)';
            OS_with_unitIDs(OS_with_unitIDs == 0) = [];

            % For all units overlaping with this "oversorted" unit
            for OS_with_unitID = OS_with_unitIDs,                        
                sortID_OS_with_unitIndex = ismember(rearranged_sortCode, OS_with_unitID);
                rearranged_sortCode(sortID_OS_with_unitIndex) = OS_unitID;
                OVERSORTED_UNITS(OVERSORTED_UNITS == OS_with_unitID) = [];
            end
            OVERSORTED_UNITS(OVERSORTED_UNITS == OS_unitID) = [];

            UNSORTED_UNITS = [UNSORTED_UNITS OS_unitID];
        end

        UNSORTED_UNITS = unique(UNSORTED_UNITS);

        
        %% PLot?
        if DEBUG_PLOTS, plotPCA_WF(rearranged_sortCode, this.wf); end


        %% Continue Sorting?  
        all_sortIDs = rearranged_sortCode;
        diff_btw_prev_and_updated(end+1) = calculateACC(prev_all_sortIDs, all_sortIDs)

        if length(diff_btw_prev_and_updated) > 2,
            END_REACHED = all(diff_btw_prev_and_updated(end-1:end) == 0);
        else
            END_REACHED = false;
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
    
    
    sortClassName = class(this.SoMObj);
    sortMethodLabel = eval(sprintf('%s.sortMethodLabel',sortClassName) );
    sortLabel = sprintf('%s: %s', sortMethodLabel, 'Iterative');
    
    scoreClassName = class(this.ScMObj);
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


%% Helper Functions
function [updated_sortCode, nextID] = sortData(this, current_sortCode, UNITS_TO_SORT, K)

    % Sort Any Unsorted Waveforms
    unsortedWf_ind = ismember(current_sortCode, UNITS_TO_SORT);
    [recent_sortIDs] = this.SoMObj.sort(this.wf(:,unsortedWf_ind), ...
                                        this.ts(unsortedWf_ind), ...
                                        this.fs, ...
                                        this.noiseEstimate, ...
                                        [], [], [], ...
                                        K);                                        

    % Reintegrate recently sorted Waveforms into final_sortCode
    nextID = max(unique(current_sortCode))+1 - min(unique(recent_sortIDs));
    current_sortCode(unsortedWf_ind) = nextID + recent_sortIDs;
    
    updated_sortCode = current_sortCode;

end

function [unitIDs, UNDERSORTED, OVERSORTED, OS_pairs] = scoreUpdatedSort(this, sortCode, nextID)

    import hsst.scoreMethod.*
        
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
    sortScoreObj = sortQualityClass(sortCode, ...
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

    OS_pairs = {};
    for unit_ind = 1:length(unitIDs),
        oversorted_Unit_Pairs = [cell2mat(sortScoreObj.raw_values{6}(unit_ind)), ...
                                 cell2mat(sortScoreObj.raw_values{5}(unit_ind))]';
        OS_pairs{unit_ind} = unique(oversorted_Unit_Pairs, 'stable');
    end           

    
    %% Sort "Logic" for resorted Units

%     % Set Noise
%     NOISE_UNITS = [unitIDs(is_noise)];
% 
%     % Set "Good" Units
%     SORTED_UNITS = [unitIDs(is_goodUnit)];

    % UnderSorted Units
    UNDERSORTED = unitIDs(is_resortedUnit & is_undersorted);

    % OverSorted Units
    OVERSORTED = unitIDs(is_resortedUnit & is_oversorted & ~is_undersorted);   
    
end



















        