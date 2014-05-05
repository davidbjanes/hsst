
% Abstract sorting class to interface with the QNAP database
classdef sorter < handle

    % Public Properties
    properties (Abstract, Constant)         
        sortMethodLabel         % I.e. 'Dukesort', 'UMS', etc
        defaultParameters
    end

    % Public Settable Properties
    properties (Access = public)
        verbose = false;
    end       

    % Abstract Methods
    methods (Abstract, Static)

        % OUTPUTS ---------------------------------------------------------
        % sort_ids      - array of ID numbers labeling each waveform
        % wf            - waveforms 
        % ts            - time stamps
        % property      - [struct] method dependant
        % INPUTS ----------------------------------------------------------
        % wf            - waveforms
        % ts            - timestamps
        % fs            - sampling frequency
        % th            - threshold/noise value
        % align_sample  - sample point waveforms are aligned at
        % dur           - duration (only needed if raw_data is given)
        % raw_data      - unfiltered recorded waveform
        % sortparameter - sort method dependant
        [sort_ids, wf, ts, property] = sort(wf, ts, fs, th, ...
                                            align_sample, dur, raw_data, ...
                                            sortParameter)
    end
    % Public Methods
    methods

        % Additional way to sort 
        function [sort_ids, wf, ts, property] = sortViaDataObject(this, dataObj, sortParameter)

            wf              = dataObj.wf_snippets;
            ts              = dataObj.time_stamps;
            fs              = dataObj.sampling_freq;
            th              = dataObj.threshold_value;
            align_sample    = dataObj.alignment_sample_num;
            dur = [];
            raw_data        = dataObj.wf_continuous;
            
            try
                [sort_ids, wf, ts, property] = this.sort(wf, ts, fs, th, ...
                                                         align_sample, dur, raw_data, ...
                                                         sortParameter);
            catch error_obj
                fprintf(error_obj.message)
                sort_ids = ones([1 length(ts)]);
            end
        end 

        function obj = sorter(verbose)

            if nargin > 1,
                obj.verbose = verbose;
            end

            if obj.verbose,
                fprintf('Sort Object Created: %s \n', obj.sortMethodLabel)
            end

        end

    end

end