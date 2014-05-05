
function this = generateParameterScores(this, sortParameters)
            
    if this.verbose,
        tic
        fprintf('Sorting/Scoring via %s... ', this.SMObj.sortMethodLabel);
    end

    % Sorting at every parameter value
    for n = 1:length(sortParameters),
        fprintf('Sorting w/Param: %d \n', sortParameters(n))
        sortCode = this.SoMObj.sortViaDataObject(this.dataObj, sortParameters(n));
        sortCode = this.cleanUpSortCode(sortCode);
        
        sortClassName = class(this.SoMObj);
        sortMethodLabel = eval(sprintf('%s.sortMethodLabel',sortClassName) );

        scoreClassName = class(this.ScMObj);
        sortLabel = sprintf('%s: %d', sortMethodLabel, sortParameters(n));

        scoreObj = feval(scoreClassName,sortCode, ...
                                        this.dataObj.wf_snippets, ...
                                        this.dataObj.time_stamps, ...
                                        this.dataObj.noise_estimate, ...
                                        'string_label', sortLabel);
        
        % Update dataObj
        this.dataObj.addSortCode( sortCode, sortMethodLabel, sortParameters(n), scoreObj );
        
    end

    if this.verbose,
        fprintf('Finished in: %2.2f sec, ', toc)
        fprintf('Scoring... ');
    end


end