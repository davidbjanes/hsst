function generateParameterScores_IterativeNoiseRemoval(this)

    sortCode = ones([length(sortParams{sort_index}) length(perfectSortCode)]);

    %TODO: ITERATIVE Noise Removal Deal
    FLAG_NOISE_TO_REMOVE = true;
    while FLAG_NOISE_TO_REMOVE,

        % Remove Known Noise Waveforms
        noise_wf_ind = all(sortCode) == 0;
        noise_removed_wf = wf(:, ~noise_wf_ind);
        noise_removed_ts = ts(~noise_wf_ind);
        sortCode(:,noise_wf_ind) = [];

        if ~isempty(sortCode),

            % Sorting
            tic
            fprintf('Sorting via %s... ', SORT_OBJ{sort_index}.sortMethodLabel);
            for i = 1:length(sortParams{sort_index}),
                [sortCode(i, :)] = SORT_OBJ{sort_index}.sortMethod(noise_removed_wf, ...
                                                                   noise_removed_ts, ...
                                                                   sample_Hz, ...
                                                                   noiseEstimate, ...
                                                                   10, [], [], ...
                                                                   sortParams{sort_index}(i));

            end
            fprintf(' Finished in: %2.2f sec, ', toc)

            % Scoring
            fprintf('Scoring... ');
            temp_error_rate = [];
            sortScoreObj = sortQualityClass([1], [1], [1]);
            for i = 1:length(sortParams{sort_index}),
                sortScoreObj(i) = sortQualityClass(sortCode(i, :), ...
                                                   noise_removed_wf, ...
                                                   noise_removed_ts, ...
                                                   noiseEstimate, ...
                                                   'string_label', sprintf('%s: %d', output_text{sort_index}, sortParams{sort_index}(i)), ...
                                                   'normalize_num_wf', true);
            end

            % Check for Noise Units            
            fprintf('Check for Noise... ');
            FLAG_NOISE_TO_REMOVE = false;
            for i = 1:length(sortParams{sort_index}),
                noiseUnitIDs = sortScoreObj(i).unitIDs(sortScoreObj(i).noise_units);
                if ~isempty(noiseUnitIDs),
                    for unitID = noiseUnitIDs,
                        ind = sortCode(i,:) == unitID;
                        sortCode(i, ind) = 0;
                    end
                    FLAG_NOISE_TO_REMOVE = true;                    
                end
            end      

        else
            FLAG_NOISE_TO_REMOVE = false;
        end
    end

    % Re-add noise waveforms back into sort, and rescore
    signal_wf_ind = ismember(ts, noise_removed_ts);
    remerged_sortCode = zeros([length(sortParams{sort_index}) length(perfectSortCode)]);
    remerged_sortCode(:,signal_wf_ind) = sortCode;
    sortCode = remerged_sortCode;

    for i = 1:length(sortParams{sort_index}),
            sortScoreObj(i) = sortQualityClass(sortCode(i, :), ...
                                               wf, ...
                                               ts, ...
                                               noiseEstimate, ...
                                               'string_label', sprintf('%s: %d', output_text{sort_index}, sortParams{sort_index}(i)), ...
                                               'normalize_num_wf', true);
        % Difference between Spike Trains / KMEANs
        temp_error_rate(i) = calculateACC(perfectSortCode, sortCode(i, :));
    end
    
end
