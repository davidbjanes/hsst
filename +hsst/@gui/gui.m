classdef gui < handle
    
    % Public Properies
    properties
        fig_handle
        h
        
    end
     
    % Global Constants
    properties (GetAccess = private, Hidden, Constant)
        DEFAULT_SORT_PARAM = 1;
        DEFAULT_UNIT_SEL = 1;
        
        GENERATED_SORT_PARAMETERS = [10^6]; %, 10^5, 10^4, 10^3, 10^2];
        
        GREY = [0.5 0.5 0.5];
        DEFAULT_BKGND_COLOR = 'k';              % Default background color
        DEFAULT_FORGND_COLOR = [0.5 0.5 0.5];   % Default foreground color
        
        %STATIC VARIABLES
        COLORS = ['w','g','b','r','c','m','y'];
        SYMBOLS = ['.','o','x','+','*','s','d','v','^','<','>','p','h'];
    end
    
    % Private Global Vars
    properties (GetAccess = private, SetAccess = private)
        % QNAP Reference Objects
        tr;
                
        status_text = 'Running...';
        
        % Default Sort Object (holds, wf, ts, sortCode, etc)
        sort = struct(...
                        'catName', '',...
                        'block', [],...
                        'channel', [],...
                        'spikeDataObj', 1,...
                        'wf', [],...
                        'ts', [],... 
                        'code', [],...
                        'label', '',...
                        'score', [],...
                        'valid', true);
        sortList
        cur_index
        
        % List Box - Sort Select
        selected_SDO      % Spike Data Obj
        selected_sort     % Sort index
        
        % Variables Currently Displayed on the GUI 
        curSortCode
        curSortScore
        validSortCodes          % TOBE discontinued
        
        % Handles to Waveform Line Objects, and PCA scatterplots
        % Syntax: handle_object{sortParam}(unit number)
        h_wf;
        h_pca;
        h_text_wf;
        h_text_pca;
        h_metricLegend = 0.00001;
        h_metricButtons = [];
        
        % Metric Display
        selected_metric_index = 5;
        colormapping = [0.5 0 0; ...
                        1 0   0; ...
                        1 0.5 0; ...
                        1 1   0; ...
                        0 1   0; ...
                        0 1   1; ...
                        0 0   1; ...
                        0.5 0 1; ...
                        1 0   1; ...
                        0.2 0.2 0.2];
        
        % Unit Information on Display on the GUI 
        unitFormat = {};
        unitLabels = [];
        unit
        unitIDs
        unitVisibility 
        unitHGGroup_WF
        unitHGGroup_PCA
        selected_units_unitListBox                
        
        INITIALIZED_FLAG = false;
        
        % DATABASE CONNECTION
        QNAP_online = true;
        filepathToLOADED_DATA
        loaded_data_struct
        
        % Test Purposes
        test = 'Refresh Done: ''hello world''';
        verbose = true;
        debug = false; 
    end
    
    
    % INTERNAL METHODS
    % =====================================================
    methods
        
        % Constructor -----------------------------------------------------
        function obj = gui(dataObj)
            
            %quick cheap dirty code
            class_path = which('hsst.gui');
            gui_path   = fullfile(fileparts(class_path), 'hsstGUI.fig'); 
            
            uiopen(gui_path, 1);
            obj.fig_handle = gcf;
            obj.h = guihandles(obj.fig_handle);    
            setappdata(obj.fig_handle,'obj',obj);
            
            % UnImplemented Features
%             set(obj.h.import_workButton, 'Enable', 'off');
            set(obj.h.import_QNAPButton, 'Enable', 'off');
            set(obj.h.export_OfflineSorterButton, 'Enable', 'off');
                        
            % Set CallBack functions
            % Input Panel
%             set(obj.h.import_QNAPButton,'Callback',@obj.import_QNAP_CB);
            set(obj.h.import_workButton,'Callback',@obj.import_work_CB);
            set(obj.h.clearAllPlotsButton,'Callback',@obj.clearAllPlotsButton_CB);
            
            % Select Sort Panel
            set(obj.h.spikeDataObjListBox,'Callback',@obj.spikeDataObjListBox_CB);
            set(obj.h.sortListBox,'Callback',@obj.sortListBox_CB);
            set(obj.h.generateDSButton,'Callback',@obj.generateDSButton_CB);
            set(obj.h.compareButton,'Callback',@obj.compareButton_CB);
            
            % Export Panel
            set(obj.h.export_workButton,'Callback',@obj.export_workButton_CB);
            set(obj.h.export_QNAPButton,'Callback',@obj.export_QNAPButton_CB);
            set(obj.h.export_OfflineSorterButton,'Callback',@obj.export_OfflineSorterButton_CB);

            % Units Panel
            set(obj.h.unitListbox,'Callback',@obj.unitListbox_CB);
            set(obj.h.combineButton,'Callback',@obj.combineButton_CB);
            set(obj.h.uncombineAllButton,'Callback',@obj.uncombineAllButton_CB);
            
            % Build Metric Buttons
            import hsst.scoreMethod.*
            temp_sortQualityClass_instance = sortQualityClass(1,1,1,1);
            for i = 1:temp_sortQualityClass_instance.num_of_metrics,
                obj.h_metricButtons(i) = uicontrol(... obj.h, 'metric_1', ...
                        'String', sprintf('%d', i), ...
                        'Style', 'pushbutton', ....
                        'Units', 'pixels', ...
                        'Position', [325+(i-1)*40 50 35 40], ...
                        'BackgroundColor', obj.colormapping(i,:), ...
                        'Callback',@obj.selectMetric_CB);
            end
            
            % Exit Button
            set(obj.h.exitButton,'Callback',@obj.exitButton_CB);
            
            % Disable All Buttons Before Data is Imported
            useGUI(obj, 'restart');
            updateStatus(obj, 'Waiting for User Input...');
            
            if nargin > 0,
                reinitialize_base_workspace(obj, dataObj)
                updateWF_PCA_PLOTS(obj);
            end
        end
        
        % Helper Function to Constructor ----------------------------------
        function status = reinitialize_base_workspace(this, dataObj)

            this.filepathToLOADED_DATA = '';
            this.sortList = [];
            
            % LOADING FROM WORKSPACE --------------------------------------
                
            this.QNAP_online = false;
            
            try     
                if nargin < 2,
                    [FileName,PathName,~] = uigetfile;
                    fullFilePath = [PathName FileName];
                    data = load(fullFilePath);
                    field_strings = fieldnames(data);
                    dataObj = data.(field_strings{1});
                else
                    fullFilePath = 'Loaded From Memory';
                end
            
                % Check if object is a class "dataObject"
                classStrName = class(dataObj);
                if ~strcmp(classStrName,'hsst.dataObject')
                    error('Wrong Class, variable passed in must be an instance of the class "hsst.dataObject()"')
                end
                
                N = length(dataObj);
                if N > 1, 
                    error('Only 1 object is allowed to be passed, no arrays.')
                end
                
                wf              = dataObj.wf_snippets;
                ts              = dataObj.time_stamps;
                noiseEstimate   = dataObj.noise_estimate;

                label = dataObj.dataLabel;

                this.filepathToLOADED_DATA = fullFilePath;
                this.loaded_data_struct = dataObj;

                num_of_sorts = dataObj.num_sortCodes;

            catch error_object
                this.INITIALIZED_FLAG = false;

                msgbox(error_object.message, 'LOAD ERROR!');

                status = false;
                return;
            end

            if num_of_sorts == 0,
                addSort(this, wf, ts, noiseEstimate, [], [], {label}); 
            else             
                for n = 1:num_of_sorts,
                    newsort_label = {dataObj.sortScoreObjList{n}.label};
                    newsort_sortCode = dataObj.sortCodeList{n};
                    newsort_score = dataObj.sortScoreObjList{n};
                    
                    addSort(this, wf, ts, noiseEstimate, newsort_sortCode(:)', ...
                            newsort_score, newsort_label);     
                end
            end                
                
            %% LOADING FROM QNAP DATABASE ---------------------------------------
%             else
%                 
%                 if isempty(this.chann) || isempty(this.block),
%                     error('Block Number or Channel Number Invalid.')
%                 end
% 
%                 % Check if QNAP DATABASE is connected
%                 try
%                     % Get Waveform and Timestamp information from QNAP
%                     this.tr = getDBObj('trial', this.catName, this.block);
% 
%                     % Variables for saving GENERATED DATA
%                     [C ~] = setupConvPathForCat(this.catName);
%                     createFolderIfNoExist(C.ANALYSIS_PATHS.DUKESORT_ROOT);
% 
%                     blockName = sprintf('block-%d', this.block);
%                     filePath = C.ANALYSIS_PATHS.DUKESORT_ROOT;
%                     filePath = [filePath, '\', blockName];
% 
%                     fullFilePath = [filePath, '\', blockName, sprintf('-%d.mat', this.chann)];
% 
%                     this.QNAP_online = true;
% 
%                 % QNAP Offline
%                 catch EM
%                     this.QNAP_online = false;
% 
%                     user_resp = questdlg(['Connection Error : QNAP database not found!                              ' ...
%                                           'For SORTED data, please select appropiate *.mat file.                    ' ...
%                                           'For UN-SORTED data, please cancel and select load from workspace.'], ...
%                                          'WARNING!', 'Continue', 'Cancel', 'Cancel');
% 
%                     fullFilePath = 0;
%                     if strcmp('Continue', user_resp)
%                         fullFilePath = uigetfile;
%                     end
% 
%                     % If User Cancels, cancel whole operation.
%                     if fullFilePath == 0 | strcmp('Cancel', user_resp),
%                         this.INITIALIZED_FLAG = false;
% 
%                         this.catName_changed = false;
%                         this.block_changed = false;
%                         this.chann_changed = false;
% 
%                         status = false;
%                         return;
%                     end
% 
%                 end
% 
% 
%                 % Attempt to LOAD GENERATED DATA 'block' exists.  Load it.
%                 if exist(fullFilePath, 'file') == 2,
%                     load(fullFilePath);
%                     if this.verbose, fprintf('Loaded ''%s'' file. \r', fullFilePath); end
% 
%                     GEN_DATA = block;  %#ok<*CPROP,*PROP>
% 
%                 % LOAD FAILED: 'block' doesn't exist.  Create a temporary variable 
%                 else           
%                     GEN_DATA = [];  
%                 end
% 
% 
%                 % Initialize sortList
%                 this.sortList = [];
% 
% 
%                 % Load Data from GENERATED_DATA (if it exists)
%                 if ~isempty(GEN_DATA),
% 
%                     existingSpikeDataObj = [1:size(GEN_DATA.sortCodes, 2)]; %#ok<*NBRAK>
%                     for spikeDataObj = existingSpikeDataObj,
% 
%                         if ~isempty(GEN_DATA.sortCodes{spikeDataObj}),
% 
%                             existingSortCodes = [1:length(GEN_DATA.sortCodes{spikeDataObj}(1,:))];
%                             for sort_index = existingSortCodes,
%                                 
%                                 sortCode = GEN_DATA.sortCodes{spikeDataObj}{1, sort_index};
%                                 
%                                 wf = this.tr.spikeData(spikeDataObj).st(this.chann).snips.wf;
%                                 ts = this.tr.spikeData(spikeDataObj).st(this.chann).ts;
%                                 noiseEstimate = this.tr.rawData.wf(this.chann).noiseEstimate;
%                                 
%                                 score = GEN_DATA.sortScores{spikeDataObj}(sort_index);
%                                 
%                                 label = strcat('DukeSort: ', GEN_DATA.sortCodes{spikeDataObj}{2, sort_index});
%                                 
%                                 addSort(this, this.catName, this.block, this.chann, spikeDataObj, ...
%                                     wf, ts, noiseEstimate, sortCode, score, ...
%                                     label);
%                             end
%                         end
%                     end
% 
%                 end
% 
% 
%                 % if DATABASE is connected, find "All Sorts" in Database
%                 if this.QNAP_online,
% 
%                     validSpikeDataObj = [1:length(this.tr.spikeData)];
%                     for spikeDataObj = validSpikeDataObj,
% 
%                         validSortCodeNames = this.tr.spikeData(spikeDataObj).getAvailableSortCodes;
%                         for sortName = validSortCodeNames,
% 
%                             sortCode = this.tr.spikeData(spikeDataObj).st(this.chann).getSortIDs(sortName);
%                             
%                             wf = this.tr.spikeData(spikeDataObj).st(this.chann).snips.wf;
%                             ts = this.tr.spikeData(spikeDataObj).st(this.chann).ts;
%                             noiseEstimate = this.tr.rawData.wf(this.chann).noiseEstimate;
%                             
%                             score = [];
%                             
%                             label = sortName;
%                             
%                             addSort(this, this.catName, this.block, this.chann, spikeDataObj, ...
%                                     wf, ts, noiseEstimate, sortCode, score, ...
%                                     label);        
%                         end
%                     end
%                 end

%             end % ---------------------------------------------------------

            
            % ListBox Selection Defaults 
            this.selected_SDO = this.DEFAULT_SORT_PARAM;     
            this.selected_sort = this.DEFAULT_SORT_PARAM;


            % Set current sort code and current sort score
            setSortParam(this, this.selected_SDO, this.selected_sort);
            this.curSortCode = this.sortList(this.cur_index).code;
            this.curSortScore = this.sortList(this.cur_index).score;

            
            % Set ListBox Selection to Default Value
            this.selected_units_unitListBox = this.DEFAULT_UNIT_SEL;
            

            % Initialize GUI
            initializeGUI(this);
            
            status = true;
            return;
        end
        
        % Add Sort to SORT LIST
        function new_sort = addSort(this, wf, ts, noiseEstimate, sortCode, score, label)
        
            new_sort = this.sort;
            new_sort.valid = true;

            new_sort.wf = wf;
            new_sort.ts = ts;
            new_sort.noise_est = noiseEstimate;
            
            new_sort.label = label;
            new_sort.spikeDataObj = 1;
            
            % Lump all waveforms in same unit if unsorted
            if isempty(sortCode)
                sortCode = ones([length(ts), 1]);
            end
            
            % If Waveform length or Ts or SortCode lengths aren't
            % equal
            if length(sortCode) == length(new_sort.ts) && ...
                        length(sortCode) == length(new_sort.wf),

                % Get from Database and generate sortQuality scores
                new_sort.code = sortCode; 
                
                import hsst.scoreMethod.*
                if isempty(score) || score.VERSION_NUMBER ~= sortQualityClass.MOST_CURRENT_VERSION_NUMBER,
                    
                    [scoreObject] = sortQualityClass(sortCode, ...
                                                     new_sort.wf, ...
                                                     new_sort.ts, ...
                                                     new_sort.noise_est, ...
                                                     'string_label', label, ...
                                                     'normalize_num_wf', true);
                    score = scoreObject;
                end
                
                new_sort.score = score;
                                
                % Index for next insertion of sort into List
                save_sort_index = length(this.sortList) + 1;

                % Initialize Array if array is empty, or append if
                % array is already initialized
                if isempty(this.sortList),
                    this.sortList = new_sort;
                else
                    this.sortList(save_sort_index) = new_sort;
                end
            
            % Invalid Import
            else
                new_sort.valid = false;
            end
        
        end       
                
        % Initalize All GUI Objects ---------------------------------------
        function initializeGUI(this)
            
            % Local Temp Variables
            waveformAxis = this.h.waveformAxis;
            pcaAxis      = this.h.pcaAxis;
            wf           = this.sortList(this.cur_index).wf;

            
            % WAVEFORM LINES
            set(this.fig_handle, 'CurrentAxes', waveformAxis); cla
            set(waveformAxis, 'Color', this.DEFAULT_BKGND_COLOR);
            grid on
            set(gca, 'XColor', this.DEFAULT_FORGND_COLOR);
            set(gca, 'YColor', this.DEFAULT_FORGND_COLOR);
            hold all
            
            timeStepVector = repmat([1:size(wf,1)]', 1, size(wf, 2));
            wf_handles = line(timeStepVector, wf(:,:));
            xlim([1 size(wf,1)]);

            
            % PCA POINTS
            set(this.fig_handle, 'CurrentAxes', pcaAxis); cla
            set(pcaAxis, 'Color', this.DEFAULT_BKGND_COLOR);
            grid on
            set(gca, 'XColor', this.DEFAULT_FORGND_COLOR);
            set(gca, 'YColor', this.DEFAULT_FORGND_COLOR);
            
            [~, score{1}] = princomp(wf');
            PCA_2D_point_X = repmat(score{1}(:,1)', 2, 1);    
            PCA_2D_point_Y = repmat(score{1}(:,2)', 2, 1); 
            pca_handles = line(PCA_2D_point_X, PCA_2D_point_Y, 'LineStyle', '*');
                       
            
            % Update Existing SpikeDataObject (from database) ListBox            
            [validSpikeDataObj, indices] = unique([this.sortList.spikeDataObj]);
            set(this.h.spikeDataObjListBox, 'String',  [this.sortList(indices).label]);
            set(this.h.spikeDataObjListBox, 'Value', this.selected_SDO);
            
            % Update Existing SpikeTimes (from the default selected 
            % SpikeDataObject above) ListBox   
            indices = [this.sortList.spikeDataObj] == validSpikeDataObj(this.selected_SDO);
            set(this.h.sortListBox, 'String',  [this.sortList(indices).label]);    
            set(this.h.sortListBox, 'Value', this.selected_sort);
            
            
            % Save handles to global varibles
            this.h_wf = wf_handles;
            this.h_pca = pca_handles;
               
            
            % Renable All Button functions
            if ~this.INITIALIZED_FLAG, 
                useGUI(this, 'on');
                
                this.INITIALIZED_FLAG = true;
            end

        end
                
        % Setter Function for Additive Noise Parameter
        function setSortParam(this, spikeDataObj_index, sort_index)
            
            % Update global varibles
            this.selected_SDO = spikeDataObj_index;    % Spike Data Obj
            this.selected_sort = sort_index;           % Sort index
            
            % Find all indices of selected spike data object
            [validSpikeDataObj] = unique([this.sortList.spikeDataObj]);
            indices = find([this.sortList.spikeDataObj] == validSpikeDataObj(spikeDataObj_index));
            
            % update global variables of current sort and score
            this.cur_index = indices(sort_index);
            this.curSortCode = this.sortList(this.cur_index).code;
            this.curSortScore = this.sortList(this.cur_index).score;
            
            % change sort code to column vector if needed
            if size(this.curSortCode, 1) > size(this.curSortCode, 2),
                this.curSortCode = this.curSortCode';
            end   
            
        end
        
        % Setter Function for Additive Noise Parameter
        function setListBoxSelection(this, k)
            
            % update global varible
            this.selected_units_unitListBox = k; 

        end
        
        % Updates Waveform and PCA plots
        function updateWF_PCA_PLOTS(this)
                        
            % Initialize Base Workspace if flags have been tripped
            if false,
%                 this.block_changed || this.chann_changed || this.catName_changed,
                status = reinitialize_base_workspace(this, false);
            
                if ~status,
                    return;
                end
            else
            
                % Local Temp Variables
                curSortCode = this.curSortCode;
                curSortScore = this.curSortScore;
%                 figure(this.fig_handle);
%                 gcf = this.fig_handle;
                
                % Generate Colors for Waveform and PCA plots
                % Organize Waveforms Lines & PCA points, Update GUI on
                % colors and graphing order
                
                [this.unitIDs, ~] = unique_2011b(curSortCode);
                this.unitLabels = [];   
                this.unitFormat = {};

                number_wf_perUnit = this.curSortScore.num_wf;
                sort_result = sortrows([number_wf_perUnit; double(this.unitIDs)]', 1);
                defaultRankOrder = sort_result(:, 2)';
                for unit_Id = defaultRankOrder,
                    i = find(this.unitIDs == unit_Id);
                    
                    this.unitFormat{i} = strcat(this.COLORS(1+rem(i, numel(this.COLORS))), ...
                                             this.SYMBOLS(1+rem(3*i, numel(this.SYMBOLS))));
                    
                    set(this.h_wf(curSortCode == unit_Id), 'Color', this.unitFormat{i}(1), ...
                                                       'Visible', 'on');
                    set(this.h_pca(curSortCode == unit_Id), 'Color', this.unitFormat{i}(1), ...
                                                        'LineStyle', this.unitFormat{i}(2:end));

                    uistack(this.h_wf(this.curSortCode == unit_Id), 'bottom');
                    uistack(this.h_pca(this.curSortCode == unit_Id), 'bottom');
                    
                end  
                
                % Set all Units Visiblity == true
                this.unitVisibility = repmat({'on'}, length(this.unitIDs), 1);
                
                
                % Update WAVEFORM UNIT Text
                text_handle = zeros(1, length(this.unitIDs));
                textspace = 0.065;
                for unit_Id = this.unitIDs, 
                    i = find(this.unitIDs == unit_Id);
                    
                    this.unitLabels = [this.unitLabels; {strcat('Unit: ', num2str(unit_Id))}];
                    
                    text_handle(i) = text(0.95-((i-1)*textspace), 0.99, ...
                                          sprintf('%s\n%d', this.unitLabels{i}, curSortScore.num_wf(i)), ...
                                          'Color',           this.unitFormat{i}(1), ...
                                          'HorizontalAlign', 'center', ...
                                          'VerticalAlign',   'top',    ...
                                          'parent',          this.h.waveformAxis, ...
                                          'Units',           'normalized' ...
                                          );
                                      
                    set(text_handle(i), 'ButtonDownFcn', @sortQualityEvalGUI.unitVisiblity_CB);
                    set(text_handle(i), 'UserData', i);
                end
                if ishandle(this.h_text_wf)
                    delete(this.h_text_wf);
                end
                this.h_text_wf = text_handle;
                
                
                % Update PCA Legend
                text_handle = zeros(1, length(this.unitIDs));
                textspace = 0.03;
                for unit_Id = this.unitIDs, 
                    i = find(this.unitIDs == unit_Id);
                    
                    text_handle(i) = text(0.975, 0.99-((i-1)*textspace), ...
                                          sprintf('%s %s', this.unitFormat{i}(2), this.unitLabels{i}), ...
                                          'Color',           this.unitFormat{i}(1), ...
                                          'HorizontalAlign', 'right', ...
                                          'VerticalAlign',   'top',    ...
                                          'parent',          this.h.pcaAxis, ...
                                          'Units',           'normalized' ...
                                          );
                end
                if ishandle(this.h_text_pca)
                    delete(this.h_text_pca);
                end
                this.h_text_pca = text_handle;
  
                
                % Update List Unit 
                if any(this.selected_units_unitListBox > length(this.unitLabels)),
                    set(this.h.unitListbox, 'Value', 1);
                    setListBoxSelection(this,1);
                end
                set(this.h.unitListbox, 'String', this.unitLabels);
                
                % Would be nice to add feature Coloring units in the list
                % box:
                % this.unitLabels = [this.unitLabels; {sprintf('<HTML><BODY bgcolor="%s">%s', ...
                %                      this.unitFormat{i}(1), ...
                %                      strcat('Unit: ', num2str(unit_Id)))}];
           
                
                % Update SpikeTimes ListBox (from selected SpikeDataObj) 
                indices = [this.sortList.spikeDataObj] == this.sortList(this.cur_index).spikeDataObj;
                set(this.h.sortListBox, 'String',  [this.sortList(indices).label]);    
                set(this.h.sortListBox, 'Value', this.selected_sort);
                
                
                % Sort Quality Plot
                scoreObject = this.curSortScore;
                metric_score = zeros(scoreObject.num_of_metrics, scoreObject.num_units);
                for i = 1:length(this.unitIDs), 
                    metric_score(:,i) = scoreObject.getUnitScoreByMtrc(this.unitIDs(i));
                end
                
                figure(this.fig_handle);
                set(this.fig_handle, 'CurrentAxes', this.h.varAxis);
                cla
                
                updateSortScoreMatrix(this, metric_score, this.unitIDs, scoreObject.metric_labels);
           
                set(this.h.totalScoreText, ...
                    'String', sprintf('Current Sort Score: %0.3f', scoreObject.score));
                               
                % Update Sort Quality Plots
                updateSortQuality_PLOTS(this);
                
            end
            
        end
        
        % TODO: Document
        function updateSortScoreMatrix(this, score_matrix, y_axis_labels, legend_str)
        % MxN, score_matrix, N - number of units, M - number of metrics
                               
            for x = 1:size(score_matrix,1)
                patch([0 0], [0 0], this.colormapping(x, :));
            end

            h_legend = legend(legend_str, 'Location','NorthEastOutside');
            set(h_legend,'FontSize',8);

            gap = 0.1;
            score_matrix = fliplr(score_matrix);
            for x = 1:size(score_matrix,1)
                for y = 1:size(score_matrix,2)

                    if score_matrix(x, y), 
                        color = this.colormapping(x, :);
                    else
                        color = 'w';
                    end

                    h = patch([.5+gap+(x-1) 0.5+gap+(x-1) 1.5-gap+(x-1) 1.5-gap+(x-1)], ...
                              [(y-1)    y         y         (y-1)], ...
                              color);

                end
            end

            axis([0.5 size(score_matrix,1)+0.5 0 size(score_matrix,2)])
            set(gca, 'YTick', [0.5:size(score_matrix,2)])
            set(gca, 'YTickLabel', [fliplr(y_axis_labels)])
            set(gca, 'XTickLabel', {''})
            
        end
        
        % Updates ISI and Current Sort Quality of Selected Unit plots
        function updateSortQuality_PLOTS(this)
            
            % Local Temp Variables
            isiAxis = this.h.isiAxis;
            unit_index = this.selected_units_unitListBox;
            set(this.h.unitListbox, 'Value', this.selected_units_unitListBox);

            % Sort Quality Metric
            selected_unitID = this.unitIDs(unit_index(1));
            
            try
                this.curSortScore.plotMetric(selected_unitID, ...
                          	this.selected_metric_index, isiAxis);
            
            % Handle Error for Single waveform in Unit: results in invalid
            % ISI Analysis
            catch ME
%                 if verbose, sprintf(ME.message); end
                fprintf(sprintf('%s, %s \n', ME.identifier, ME.message));
                set(this.fig_handle, 'CurrentAxes', isiAxis); cla
                str1(1) = {'Invalid Metric Analysis'};
                v = axis;
                text(abs(v(1)-v(2))*.05, v(4)*.9, str1, 'HorizontalAlignment', 'left');
                
                rethrow(ME);
            end
            
            % Sort Quality Score                
            set(this.h.curUnitScoreText, ...
                    'String', sprintf('Current Unit Score: %0.3f', ...
                    this.curSortScore.getUnitScore(selected_unitID)));

        end
        
        % Groups Units together to create new Sort Code
        function groupUnitIDs(this)
            
            curSortCode = this.curSortCode;
            
            selected_units_ind = this.selected_units_unitListBox;
            
            base_unit = this.unitIDs(selected_units_ind(1));
            
            for unit_index = this.selected_units_unitListBox(2:end),
                unit = this.unitIDs(unit_index);
                indices = curSortCode == unit;
                curSortCode(indices) = base_unit;
            end
            this.unitLabels(selected_units_ind(2:end)) = [];
            this.unitFormat(selected_units_ind(2:end)) = [];
            this.unitIDs(selected_units_ind(2:end))    = [];
            this.setListBoxSelection(selected_units_ind(1));
           
            import sqm.scoreMethod.*
            [scoreObject] = sortQualityClass(curSortCode, ...
                                 this.sortList(this.cur_index).wf, ...
                                 this.sortList(this.cur_index).ts, ...
                                 this.sortList(this.cur_index).noise_est, ...
                                 'normalize_num_wf', true);
                                                           
            this.curSortScore = scoreObject;
            this.curSortCode = curSortCode;
        end
        
        % Updates Unit Waveform Visiblity
        function updateUnitVisibility(this, k)
            
            unit_Id = this.unitIDs(k);
            curSortCode = this.curSortCode;
            
            if strcmpi(this.unitVisibility{k}, 'on')
                this.unitVisibility{k} = 'off';
                set(this.h_text_wf(k), 'Color', this.GREY);
            else
                this.unitVisibility{k} = 'on';
                set(this.h_text_wf(k), 'Color', this.unitFormat{k}(1));
            end
            
            set(this.h_wf(curSortCode == unit_Id), 'Visible', this.unitVisibility{k});
%             set(this.unitHGGroup_WF(k), 'Visible', this.unitVisibility{k});

        end
        
         % Exports sort struct to workspace and saves as *.mat file
        function exportSortCodeToWorkspace(this, name)
               
            export_sort = this.sortList(this.cur_index);
            
            export_sort.sortCode = this.curSortCode;
            export_sort = rmfield(export_sort, 'code');
            
            export_sort.score = this.curSortScore;
            export_sort.label = name;
            
            export_sort.noiseEstimate = export_sort.noise_est;
            export_sort = rmfield(export_sort, 'noise_est');
            
            save(name, 'export_sort');            
        end
        
        % TODO: Document
        function exportSortCodeToQNAP(this, name)
               
            if this.QNAP_online,
                export_sort = this.sortList(this.cur_index);

                sortCode = {this.curSortCode};
                spikeDataObj = export_sort.spikeDataObj;
                channel = export_sort.channel;

                saveSortCodes(this.tr.spikeData(spikeDataObj), ...
                          name, ...
                          sortCode, ...
                          [channel]);      
            else
                msgbox('Connection Error : QNAP database not found!', 'WARNING');
            end
                
        end
        
        % Update GUI Status Bar
        function updateStatus(this, arg)
            
           set(this.h.statusText, 'String', arg); 
           
           drawnow update
           
        end
        
        % Enable/Disable Buttons during GUI 'busy'
        function useGUI(this, arg)
            
            restart_flag = 0;
            if strcmp(arg, 'restart')
                arg = 'off';
                restart_flag = 1;
            end
            
            if strcmp(arg, 'on'),
                updateStatus(this, 'Idle...');
            else
                updateStatus(this, 'Busy!');
            end
            
            % All User Input functions
            handles = [...
                this.h.spikeDataObjListBox, ...
                this.h.sortListBox , ...
                this.h.unitListbox  , ...
                this.h.combineButton , ...
                this.h.uncombineAllButton , ...
                this.h.export_workButton , ...
                this.h.export_QNAPButton , ...
                this.h.clearAllPlotsButton , ...
                this.h.generateDSButton, ...
                this.h.compareButton, ...
                this.h.import_workButton, ...
                this.h_metricButtons]; 
%                 this.h.import_QNAPButton, ...
%                 this.h.export_OfflineSorterButton, ...

            set(handles,'Enable', arg); 
            
            if restart_flag
%                 set(this.h.import_QNAPButton,'Enable', 'on'); 
                set(this.h.import_workButton,'Enable', 'on'); 
            end
            
        end
        
        % TODO: Document
        function generateDukeSort(this)
            
            % Find all indices of selected spike data object
            [validSpikeDataObj] = unique([this.sortList.spikeDataObj]);
            spikeDataObj = validSpikeDataObj( this.selected_SDO );
            
            if strcmp(questdlg('Sort Current Waveforms with DukeSort?  WARNING: May Take Upwards of 5-20 Minutes...','Generate Sort Codes Warning','yes','no','no'),'yes'),
                if this.QNAP_online,
                    
                    dukeSpikeSorting_wrapper('Vinnie', ...
                                         this.block, ...
                                         this.chann, ...
                                         'spikeDataObj', spikeDataObj, ...
                                         'verbose', true, ...
                                         'save_to_GENERATED_DATA', true, ...
                                         'force_replace', true);
                                     
                    % Reboot Workspace
                    status = reinitialize_base_workspace(this, true);
                    
                else
                    
                    wf = this.sortList(this.cur_index).wf;
                    ts = this.sortList(this.cur_index).ts;
                    noiseEstimate = this.sortList(this.cur_index).noise_est;
                    
                    [min_size, ind] = min(size(wf));
                    num_params = length(this.GENERATED_SORT_PARAMETERS);
                    
                    DUKEsortCode = ones([max(size(wf)), 1]);
                    
                    for k = 1:num_params,
                        DUKEsortCode = dukeSpikeSorting(wf, ...
                                         'ANP', this.GENERATED_SORT_PARAMETERS(k), ...
                                         'K', min_size, ...
                                         'NegativePeak', round(min_size/3), ...
                                         'verbose', true);
                       
                                     
                        addSort(this, [], [], [], [], ...
                                    wf, ts, noiseEstimate, DUKEsortCode, [], ...
                                    {sprintf('DS: Param %d', k)});
                    end
                    
                    % TODO: save generated sort back to workspace???

                end
                
                % Update GUI
                updateWF_PCA_PLOTS(this);

            end
            
        end
        
        % TODO: Implement and Document
        function compareSortCodeScores(this)
                        
            % Get Indices of SpikeTimes within selected SpikeDataObj
            [validSpikeDataObj] = unique([this.sortList.spikeDataObj]);
            indices = find([this.sortList.spikeDataObj] == validSpikeDataObj(this.selected_SDO));

            scoreObject = [];
            sort_types = {};
            all_sort_types = {};
            
            for i = [1:length(indices)],
                if (i == 1),
                    scoreObject = this.sortList(indices(i)).score;
                else
                    scoreObject(i) = this.sortList(indices(i)).score;
                end
                
                [new_sort_type,~] = strsplit(this.sortList(indices(i)).label{1}, ':');
                all_sort_types{i} = new_sort_type{1};
                if ~ismember(all_sort_types{i},sort_types),
                    sort_types{length(sort_types)+1} = all_sort_types{i};
                end
            end
            
            [selected_sort_type] = strsplit(this.sortList(this.selected_sort).label{1}, ':');
            sort_type_ind = strcmp(all_sort_types,selected_sort_type{1});
            scoreObject = scoreObject(sort_type_ind);
            
            number_sorts = length(scoreObject);

            score = zeros(1, number_sorts);
            percent_good_units = [];

            for k = 1:number_sorts,   
                score(k) = scoreObject(k).score;
                metric_score = zeros(1, scoreObject(k).num_of_metrics);
                for i = 1:scoreObject(k).num_of_metrics, 
                    metric_score(i) = scoreObject(k).getMetricScore(i);
                end
                percent_good_units(k, :) = metric_score;       
            end
            
            if strcmp('char', class(scoreObject(1).label) ) 
                SortLabels = {scoreObject.label};
            elseif strcmp('char', class(scoreObject(1).label{1}) ) 
                SortLabels = [scoreObject.label];
            else
                SortLabels = [];
            end

            figure
            subplot(2,1,1)
                hold on
                grid on
                plot(score, 'o', 'Color', 'b')
                plot(find(score == max(score)), ...
                                    max(score), '*', 'Color', 'r', 'LineWidth', 10)
                set(gca, 'XTick', [1:number_sorts])
                axis([0.5 (number_sorts+0.5) 0 1])
                xlabel('Tuning Parameter Values')
                set(gca, 'XTickLabel', SortLabels)
                ylabel('Score [Range: 0-1]')


            subplot(2,1,2)
                grid on
                bar(percent_good_units(:,:))
                legend(scoreObject.metric_labels, ...
                       'Location','NorthEastOutside')
                xlabel('Tuning Parameter Values')
                set(gca, 'XTickLabel', SortLabels)
                ylabel('% of "Clean" Units')
       
        end
        
    end
    


    % STATIC METHODS
    % =====================================================
    methods (Static)
        
        function selectMetric_CB(button_handle,~)
            
            obj = getappdata(gcbf, 'obj');
            
            obj.selected_metric_index = str2int(get(button_handle, 'String'));
            
clc
updateSortQuality_PLOTS(obj);

        end
        
        % TODO: Document
        function generateDSButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            % Disable All Buttons Until Loading is Complete
            useGUI(obj, 'off');
            
            generateDukeSort(obj);
                        
            % Renable All Buttons after Loading is Complete
            useGUI(obj, 'on');
            
        end
        
        % TODO: Implement
        function export_OfflineSorterButton_CB(~,~)
            
        end
        
        % TODO: Document
        function compareButton_CB(~,~)
          
            obj = getappdata(gcbf, 'obj');
            
            compareSortCodeScores(obj);
            
        end
        
        % TODO: Document
        function unitVisiblity_CB(text_handle, ~)
            
            obj = getappdata(gcbf, 'obj');
            
            text_handle_index = get(text_handle, 'UserData');
            
            updateUnitVisibility(obj, text_handle_index);
            
        end 
            
        % Combine Button Push
        function combineButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');            
            groupUnitIDs(obj);
            
            % Update GUI
            updateWF_PCA_PLOTS(obj);
            
        end
        
        % Export current sort code to workspace as varible
        function export_workButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            name = get(obj.h.exportTextbox, 'String');
            
            exportSortCodeToWorkspace(obj, name);
            
        end
        
        % Export current sort code to QNAP (database)
        function export_QNAPButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            name = get(obj.h.exportTextbox, 'String');
            
            exportSortCodeToQNAP(obj, name);
            
        end
        
        % Uncombine All Button Push
        function uncombineAllButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            setSortParam(obj, obj.selected_SDO, obj.selected_sort);
            
            % Update GUI
            updateWF_PCA_PLOTS(obj);
            
        end
        
        % CallBack Function
        function unitListbox_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            index = get(obj.h.unitListbox, 'Value');
            
            setListBoxSelection(obj, index);
            
            updateSortQuality_PLOTS(obj);
            
        end
        
        % Clear All Button CallBack Function
        function clearAllPlotsButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            set(obj.fig_handle,'CurrentAxes',obj.h.waveformAxis);   cla
            set(obj.fig_handle,'CurrentAxes',obj.h.pcaAxis);        cla
            set(obj.fig_handle,'CurrentAxes',obj.h.isiAxis);        cla
            set(obj.fig_handle,'CurrentAxes',obj.h.sortQualityAxis);cla
            set(obj.fig_handle,'CurrentAxes',obj.h.varAxis);        cla
        
        end
        
        % Change in Sort Parameter Radio Buttons Function
        function spikeDataObjListBox_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            % Disable GUI
            useGUI(obj, 'off')
            
            index = get(obj.h.spikeDataObjListBox, 'Value');
            
            setSortParam(obj, index, obj.DEFAULT_SORT_PARAM); 
            
            % Initialize GUI
            initializeGUI(obj);   
            
            % Update GUI
            updateWF_PCA_PLOTS(obj);  
            
            % Enable GUI
            useGUI(obj, 'on')

        end
        
         % Change SpikeTimes Obj (Listbox titled 'Sort')
        function sortListBox_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            % Disable GUI
            useGUI(obj, 'off')
            
            index = get(obj.h.sortListBox, 'Value');
            
            setSortParam(obj, obj.selected_SDO, index); 
            
            % Update GUI
            updateWF_PCA_PLOTS(obj);   
            
            % Enable GUI
            useGUI(obj, 'on')

        end
              
        % TODO: Implement
        function import_work_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            % Disable All Buttons Until Loading is Complete
            useGUI(obj, 'off');
            
            % Initialization
            status = reinitialize_base_workspace(obj); 
            
            if status
                % Update GUI
                updateWF_PCA_PLOTS(obj);
                
                % Renable All Buttons after Loading is Complete
                useGUI(obj, 'on');
            else
                % Renable All Buttons after Loading is Complete
                useGUI(obj, 'restart');    
            end
        end
        
        % Import From Database Button Callback Function
        function import_QNAP_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            
            % TODO: Implement
        end
        
        % Channel Number Textbox Callback Function        
        function catTextBox_CB(textboxObj, ~)
            
            obj = getappdata(gcbf, 'obj');
            
            textboxObj = get(textboxObj);
            obj.catName = textboxObj.String;
            obj.catName_changed = true;

        end
        
        % Channel Number Textbox Callback Function        
        function channelTextBox_CB(textboxObj, ~)
            
            obj = getappdata(gcbf, 'obj');
            
            textboxObj = get(textboxObj);
            obj.chann = str2double(textboxObj.String);
            obj.chann_changed = true;

        end
        
        % Block Number Textbox Callback Function
        function blockTextBox_CB(textboxObj, ~)
            
            obj = getappdata(gcbf, 'obj');
            
            textboxObj = get(textboxObj);
            obj.block = str2double(textboxObj.String);
            obj.block_changed = true;

        end
        
        % Exit Button Callback Function
        function exitButton_CB(~,~)
            
            obj = getappdata(gcbf, 'obj');
            delete(obj);
            
            close gcbf
            
            % remove deleted handle from base workspace      
            for var = evalin('base','who')',
                value = evalin('base', var{1});
                try 
                    if eq(value, obj),
                        evalin('base', sprintf('clear %s', var{1}));
                    end
                end
            end
        end
        
    end
    
end
