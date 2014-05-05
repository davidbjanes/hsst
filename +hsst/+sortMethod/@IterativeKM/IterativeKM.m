
classdef IterativeKM < hsst.sorter
    
    properties (Constant)  
        sortMethodLabel = 'Iterative KM';
        defaultParameters = [1];
    end
    
    methods (Static)
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, ~, thresh, ...
                                                     ~, ~, ~, ...
                                                     ~)
            import hsst.scoreMethod.*
            import hsst.sortMethod.*
            import stdlib.*
                                             
            % Handle Inputs
            if size(wf, 1) ~= length(ts),
                wf = wf';
                if size(wf, 1) ~= length(ts),
                    error('Time Stamp and Number of Waveform Mismatch');
                end
            end
                                    
            num_wf = size(wf, 1);
            num_samples = size(wf, 2);
            
            DEBUG_PLOTS = true;
                
            
            %% Iterative KM
            CONTINUE_SORT_FLAG = true;
            all_sortIDs = zeros([1 num_wf]);
            
            UNSORTED_UNITS = 0;
            SORTED_UNITS = [];
            NOISE_UNITS = [];
            
            K = 1;
            diff_btw_prev_and_updated = [];
            
            while CONTINUE_SORT_FLAG,
                                
                % Previous SortCode 
                prev_all_sortIDs = all_sortIDs;
                
                % Sort Any Unsorted Waveforms
                unsortedWf_ind = ismember(all_sortIDs, UNSORTED_UNITS);
                [recent_sortIDs, SUCCESS_FLAG] = IterativeKM.sortByKM(wf(unsortedWf_ind,:), K);
                
                % Reintegrate recently sorted Waveforms into final_sortCode
                nextID = max(unique(all_sortIDs))+1 - min(unique(recent_sortIDs));
                all_sortIDs(unsortedWf_ind) = nextID + recent_sortIDs;
                    
                % Score Integrated Sort
                sortScoreObj = sortQualityClass(all_sortIDs, ...
                                                wf, ...
                                                ts, ...
                                                thresh, ...
                                                'string_label', [], ...
                                                'showplots', false);
                                 
                unitIDs         = sortScoreObj.unitIDs;               
                is_resortedUnit = unitIDs > nextID;
                
                is_noise        = sortScoreObj.noise_units;
                is_undersorted  = sortScoreObj.undersorted_units;
                is_oversorted   = sortScoreObj.oversorted_units;                
                is_goodUnit     = ~is_noise & ~is_undersorted & ~is_oversorted;
                
                oversorted_ind  = [sortScoreObj.raw_values{5}; sortScoreObj.raw_values{6}];

                if DEBUG_PLOTS, plotPCA_WF(all_sortIDs, wf); end
%                 if DEBUG_PLOTS, plotPCA_WF(recent_sortIDs, wf(unsortedWf_ind,:)); end                
                
                
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
                    
                    UNSORTED_UNITS = [UNSORTED_UNITS OS_unitID OS_with_unitIDs];
                end
                
                UNSORTED_UNITS = unique(UNSORTED_UNITS);
                
                if DEBUG_PLOTS, plotPCA_WF(all_sortIDs, wf); end
                
                
                %% Continue Sorting?  
                diff_btw_prev_and_updated(end+1) = calculateACC(prev_all_sortIDs, all_sortIDs);
                
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
            
            all_sortIDs(ismember(all_sortIDs, NOISE_UNITS)) = 0;
                        
            sort_ids = all_sortIDs;
                        
        end
        
        function [sortCode, SUCCESS_FLAG] = sortByKM(wf, K)
            
            if K > 1,
                warning('off', 'all');
                KM_SUCCESS_FLAG = false([1]);

                N = 5;
                while all(~KM_SUCCESS_FLAG) & length(KM_SUCCESS_FLAG) < 5,
                    try 
                        [sortCode, ~, ~] = kmeans(wf, K, ...
                                                  'start', 'cluster', ...
                                                  'Replicates', N);
                        KM_SUCCESS_FLAG(end+1) = true;

                    catch
                        fprintf('KMeans Failure, will try again');
                        KM_SUCCESS_FLAG(end+1) = false;
                    end
                end
                warning('on', 'all');

                if all(~KM_SUCCESS_FLAG),
                    sortCode = ones([1 num_snips]);
                end


                SUCCESS_FLAG = any(KM_SUCCESS_FLAG);
                
            else
                sortCode = ones([1 length(wf)]);
                SUCCESS_FLAG = true;
            
            end
            
        end           
        
    end
    
end
