classdef sortQualityClass < hsst.scorer
     
    % Public Properties
    properties (Constant)         
        scoreMethodLabel = 'HSST';
    end
    
    % Visible (but nonEditable) Properties
    properties (SetAccess = private)

    	VERSION_NUMBER
        
        % Inherited
        score           % Double
    
        num_units     	%#ok<*SAGROW> %Int  
        num_wf        	% Int Array
        normByWf      	% Boolean
        label         	% Char Array
        unitIDs                 % unit ID numbers
        
        % Sort Metrics
        final_scores    % Boolean Array (NxD, N-number of units, D-number of metrics)
        raw_values      % Raw Outputs from metrics
        METRICS        = {};
        metric_labels  = {};
        
        num_of_metrics  % Int 
        
        % General Labels
        noise_units         % Boolean Array
        undersorted_units   % Boolean Array
        oversorted_units    % Boolean Array
        
    end

    % Hidden Constants
    properties (Hidden, Constant)
        
        MOST_CURRENT_VERSION_NUMBER = 4.0;

    end
    
    % Hidden global vars
    properties (Access = private)
            
        property_struct     % contains other varibles to be passed to metrics
                
        % Variables
        % property_struct.wf    % waveform (NxD)
        % property_struct.ts    % timestamps (Nx1)
        sort_IDs                % sort code     
        noiseEstimate           % noise estimate
        total_num_wf            % total number of waveforms    
        unit_indices            % Cell Array - indices of each unit
        num_smpl_pts            % size of each waveform
        
        max_num_units    = inf; %(numeric) maximum number of units before score
        % is set to 0
    end
    
	% Public Methods
    methods
    
        % Constructor
        function obj = sortQualityClass(sortCode, waveforms, timeStampData, ...
                                        noiseEstimate, varargin)
            
            obj@hsst.scorer(sortCode, waveforms, timeStampData, noiseEstimate)
            
            % Variable Initialization
            obj.property_struct = buildParameterStruct(obj, varargin);
            
            % General Variables
            obj.property_struct.wf                      = waveforms;
            obj.property_struct.ts                      = timeStampData(:)';
            obj.property_struct.sort_IDs                = sortCode(:)';
            
            % Other Varibles
            obj.label           = obj.property_struct.label;
            obj.normByWf        = obj.property_struct.normByWf;
            obj.max_num_units	= obj.property_struct.max_num_units;
            obj.sort_IDs        = sortCode; 
            obj.noiseEstimate   = noiseEstimate;
            
            % Metrics 
            import hsst.scoreMethod.sortQualityMetrics.*
            obj.METRICS    = {};
            % Noise Metrics
            obj.METRICS{1}  = snr();
            obj.METRICS{2}  = single_peak();
            obj.METRICS{3}  = isi_violations();
            obj.METRICS{4}  = isi_exp_fit();
            % Over Sorting
            obj.METRICS{5}  = dissimilar_peaks();
            obj.METRICS{6}  = mean_wf_similarity();
            % Under Sorting
            obj.METRICS{7}  = cross_correlation();
            obj.METRICS{8}  = template_match();
            obj.METRICS{9}  = bimodal_pk();
            obj.METRICS{10} = threshold_slope();
            
           
            obj.num_of_metrics = length(obj.METRICS);
            for i = 1:obj.num_of_metrics,
                obj.metric_labels{i} = obj.METRICS{i}.label;
            end
            
            % Update Version Number
            obj.VERSION_NUMBER = obj.MOST_CURRENT_VERSION_NUMBER;

            calcUnitInfo(obj);
            evaluateSQMetrics(obj);
            calculateFinalScore(obj);
            
        end

        function property_struct = buildParameterStruct(this, varargin)
            
            % Initializes Constants Inputs
            varargin = sanitizeVarargin(varargin{1});
            DEFINE_CONSTANTS                                
                string_label            = 'SortCode Name: [Here]';

                normalize_num_wf        = true;	% Else normalize by number of Units

                verbose                 = false;    % Outputs to cmd line
                showplots               = false;    % Outputs to plots

                max_num_units           = inf;
                
            END_DEFINE_CONSTANTS
            
            property_struct = struct();
            property_struct.verbose                 = verbose;
            property_struct.showplots               = showplots;
            
            % Other Vars
            property_struct.label                   = string_label;
            property_struct.normByWf                = normalize_num_wf;
            property_struct.max_num_units           = max_num_units;
            
        end
        
        function plotMetric(this, unitID, metric_index, axis_handle)
            
            if nargin < 4,
                axis_handle = [];
            end
            
            unit_index = find(this.unitIDs == unitID);
            unit_wf = this.property_struct.wf(this.unit_indices{unit_index}, :);
            unit_ts = this.property_struct.ts(this.unit_indices{unit_index});
            
            this.property_struct.cur_ID = unitID;
            
            this.METRICS{metric_index}.plot_metric(unit_wf, ...
                                                   unit_ts, ...
                                                   this.noiseEstimate, ...
                                                   this.property_struct, ...
                                                   axis_handle); 
        end
                
        function score = getMetricScore(this, metric_index)
            
            if metric_index <= this.num_of_metrics,
                score = mean(this.final_scores(:,metric_index));
            else
                score = [];
            end
            
        end
        
        function score = getUnitScoreByMtrc(this, unitID)
            
            unit_index = find(this.unitIDs == unitID);
            
            % if unitID given is invalid or doesn't exist
            if isempty(unit_index),
               score = [];
            
            else
                score = this.final_scores(unit_index,:)';
                
            end            
        end
        
        function score = getUnitScore(this, unitID)
            
            score = mean(getUnitScoreByMtrc(this, unitID));
            
        end
        
        function score = getScore(this)
           
            score = this.score;
            
        end
        
    end
    
    % Helper Methods to Constructor
    methods (Hidden, Access = private)
        
        % Calculate Basic Stats on SortCode
        function this = calcUnitInfo(this)
            
            % Santitize Inputs to Sort Quality Class
            
            % Change SortCode to Row Vector
            if size(this.sort_IDs,1) > size(this.sort_IDs,2)
                this.sort_IDs = this.sort_IDs';
            end

            % If length(sortCode) != #(waveforms) || #(timestamps) == ERROR
            if size(this.property_struct.wf,1) ~= length(this.property_struct.ts)
                this.property_struct.wf = this.property_struct.wf';
            end
            
            % If negatively thresholded, invert waveforms
            if this.noiseEstimate < 0,
                this.property_struct.wf = -this.property_struct.wf;
                this.noiseEstimate = -this.noiseEstimate;
            end

            % If #(waveforms) != #(timestamps) then ERROR
            if size(this.sort_IDs, 2) ~= size(this.property_struct.wf, 1) || ...
               size(this.sort_IDs, 2) ~= length(this.property_struct.ts),
                error('Sort Code length, number of Waveforms and number of timeStamps must match!')
            end
                
            % Unit IDs
            [this.unitIDs] = unique(this.sort_IDs);
            this.num_units = max(size(this.unitIDs));

            % Number of waveforms per ID and indices of each unit
            this.num_wf = zeros(1, this.num_units);
            this.unit_indices = cell(1, this.num_units);
            for unit = 1:this.num_units,
                this.num_wf(unit) = max(size(find(this.sort_IDs == this.unitIDs(unit))));
                this.unit_indices{unit} = find(this.sort_IDs == this.unitIDs(unit));

            end 

            % Total Number of Waveforms
            this.total_num_wf = length(this.sort_IDs);

            % Waveform Sample Length
            this.num_smpl_pts = size(this.property_struct.wf,2);
            
        end
              
        % Metric Evaluation
        [this] = evaluateSQMetrics(this)
        
        % Calculate Score
        this = calculateFinalScore(this)
                               
    end
    
end