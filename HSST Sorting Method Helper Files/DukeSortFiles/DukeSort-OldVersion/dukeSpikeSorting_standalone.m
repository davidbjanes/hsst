%% Duke Sorting Algorithm Wrapper for QNAP

function [output] = dukeSpikeSorting_standalone(wf, ts, noise_est, varargin)

    %% Input Variables
    
    % TODO: varargin pass through to dukeSpikeSorting(wf, varargin)
    
    varargin = sanitizeVarargin(varargin);
    DEFINE_CONSTANTS
        % Default Additive Noise Parameters
        sortParameters          = [10^6 10^5 10^4 10^3 10^2];
        sortParametersLabels    = {'10^6', '10^5', '10^4', '10^3', '10^2'};

        % Display Print Statements to MATLAB console
        verbose                 = true;
        showplots               = false;
    END_DEFINE_CONSTANTS
   
    %% SORT
    
    data = struct;
    data.wf = wf;
    data.ts = ts;
    data.noiseEstimate = noise_est;
    
    % For each sort in the sortparameters list
    for i = 1:length(sortParametersLabels),

         % Additive Noise Parameter
        ANP = sortParameters(i);
        label = sortParametersLabels(i);

        % Sort Code Created 
        sortCode = dukeSpikeSorting(wf, ...
                                    'ANP', ANP, ...
                                    'K', min(size(wf)), ...
                                    'NegativePeak', round(min(size(wf))/3), ...
                                    'verbose', verbose, ...
                                    'showplots', showplots);

        % Sort Score generated
        [scoreObject] = sortQualityClass(sortCode, wf, ts, ...
                                        'string_label', label, ...
                                        'verbose',      false, ...
                                        'mean_threshold', noise_est);                        

        sortCode = int8(sortCode);  
                            
        data.sortCode(i, :)     = sortCode;
        data.label(i)           = label;
        data.score(i)           = scoreObject;

    end     
        
    % Index of Highest Scoring sortcode
    sortScores = [data.score.score];
    [~, data.maxScoreIndex] = max(sortScores);

    output = data;    
    
end % dukeSpikeSorting function


