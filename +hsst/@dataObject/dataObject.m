
classdef dataObject < handle
% dataObject - stores all information about data
% Written by David Bjanes (dab194@pitt.edu)
%  
% INPUTS
%   wf_c	(1xN) : raw waveform 
%   wf_s	(DxN) : raw waveform 
%   ts      (1xN) : vector of time stamps 
%   fs      (1x1) : sampling frequency (Hz)
%   th      (1x1) : value at which wf_continuous was threshold
%   ne      (1x1) : estimate of noise floor
%
% PROPERTIES (Viewable only)
%   wf_continuous       (1xN) : raw waveform 
%   wf_snippets         (DxN) : raw waveform 
%   time_stamps         (1xN) : vector of time stamps 
%   sample_freq         (1x1) : sampling frequency (Hz)
%   threshold_value     (1x1) : value at which wf_continuous was threshold
%   noise_estimate      (1x1) : estimate of noise floor
%
%   number_waveforms    (1x1) : equals N
%   number_samples      (1x1) : equals D
%   alignment_sample_num(1xM) : value between [1,D] which waveforms are 
%                               aligned about
%   sortMethodObjList   {1xP} : cell array of sortMethodObjects
%   sortParameterList   {1xP} : cell array of parameters
%   sortScoreObjList    {1xP} : cell array of scoreMethodObjects
%   sortCodeList        {1xP} : cell array of sortCodes, each (1xN)
%   num_sortCodes       (1x1} : equal to P
%

    %% Public Viewable Properties
    properties (SetAccess = public)       
        dataLabel
    end
    
    %% Public Viewable Properties
    properties (GetAccess = public, SetAccess = private)
        wf_continuous
        wf_snippets
        time_stamps
        sampling_freq
        threshold_value
        noise_estimate
        
        number_waveforms
        number_samples
        alignment_sample_num
        
        sortMethodNameList
        sortParameterList
        sortScoreObjList
        sortCodeList
        
        num_sortCodes
    end
    
    %% Public Methods
    methods
        
        function obj = dataObject(wf_c, wf_s, ts, fs, th, nE)
            obj.wf_continuous   = wf_c;
            obj.wf_snippets     = wf_s;
            obj.time_stamps     = ts;
            obj.sampling_freq   = fs;
            obj.threshold_value = th;
            obj.noise_estimate  = nE;
            
            obj.sortMethodNameList = {};
            obj.sortParameterList = {};
            obj.sortScoreObjList  = {};
            obj.sortCodeList      = {};
            obj.num_sortCodes     = 0;

            obj = formatInputs(obj);
        end
                
        function success = addSortCode(this, sortCode, sortMethodName, sortParameter, sortScoreObj)
            
            success = true;
            sortLabel = sprintf('%s: %d', sortMethodName, sortParameter);
            
            if length(sortCode) == length(this.time_stamps),               
                if ~hsst.isObjectScoreMethod(sortScoreObj),
                    import hsst.scoreMethod.*
                    sortScoreObj = sortQualityClass(sortCode(:)', ...
                                                    this.wf_snippets, ...
                                                    this.time_stamps, ...
                                                    this.noise_estimate, ...
                                                    'string_label', sortLabel);
                end
                
                this.removeSortCode(sortMethodName, sortParameter);
                
                this.num_sortCodes = this.num_sortCodes + 1;
                this.sortCodeList{this.num_sortCodes}       = sortCode(:)';
                this.sortMethodNameList{this.num_sortCodes} = sortMethodName;
                this.sortParameterList{this.num_sortCodes}  = sortParameter;

                this.sortScoreObjList{this.num_sortCodes}   = sortScoreObj;
            else
                success = false;
            end            
        end
        
        function success = removeSortCode(this, sortMethodName, sortParameter)
            
            name_ind = find(strcmp(sortMethodName, this.sortMethodNameList));
            
            if ~isempty(sortParameter),
                param_ind = find(cellfun(@(x) isequal(x,sortParameter), this.sortParameterList));
                ind = find(ismember(name_ind,param_ind), 1, 'first');
            else
                param_ind = [1:this.num_sortCodes];
                ind = find(ismember(name_ind,param_ind));
            end          
            
            sort_ind = name_ind(ind);
            
            success = false;
            if ~isempty(sort_ind)
                this.sortMethodNameList(sort_ind) = [];
                this.sortParameterList(sort_ind) = [];
                this.sortScoreObjList(sort_ind) = [];
                this.sortCodeList(sort_ind) = [];
                this.num_sortCodes = length(this.sortCodeList);
                success = true;
            end
            
        end
        
        function [optimalSortCodeList, method_names, indices] = optimalSortCodeBYMETHOD(this)
            
            method_names = unique(this.sortMethodNameList, 'stable');
            
            indices = [];
            for n = 1:length(method_names)
                method_indices = find(strcmp(this.sortMethodNameList, method_names{n})); 
                max_index = hsst.optimizer.chooseOptimalParameter(this.sortScoreObjList(method_indices));
            
                indices(end+1) = method_indices(max_index);
            end
            
            optimalSortCodeList = this.sortCodeList(indices);
            
        end
        
        function [snr] = getMeanUnitSNR(this, sortCodeIndices)
            if nargin < 2,
                sortCodeIndices = 1:this.num_sortCodes;
            end
            
            snr = {};
            for n = sortCodeIndices,
                sortCode = this.sortCodeList{n};
                unique_unitIDs = unique(sortCode);
                
                sortCode_n_snr = [];
                for unit = unique_unitIDs
                    temp_snr = mean(max(abs(this.wf_snippets(:,sortCode == unit))));
                    temp_snr = temp_snr / abs(this.noise_estimate);
                    
                    sortCode_n_snr(end+1) = temp_snr;
                end
                
                snr{n==sortCodeIndices} = sortCode_n_snr;
            end            
        end
        
        function error_struct = calculateErrorFromKnownSortCode(this, knownSortCode)
            
            N = this.num_sortCodes;
            M = length(unique(knownSortCode));
            total_error = zeros([1 N]);
            error = ones([M N]);
            acc = zeros([M N]);
            tpr = zeros([M N]);
            spc = zeros([M N]);
            
            for n = 1:N,
                [temp_total_error, ...
                 temp_error, ...
                 mapping, ...
                 temp_acc, ...
                 temp_tpr, ...
                 temp_spc] = calculateACCMatrix(knownSortCode, this.sortCodeList{n});
                
                total_error(n) = temp_total_error;
                
%                 unit_found_ind = find(temp_error ~= 1);
%                 mat_ind = sub2ind(size(temp_acc),unit_found_ind,mapping(temp_error ~= 1));
%                 
%                 acc(unit_found_ind,n) = temp_acc(mat_ind);
%                 tpr(unit_found_ind,n) = temp_tpr(mat_ind);
%                 spc(unit_found_ind,n) = temp_spc(mat_ind);
                
                [~,mapping2] = min(temp_error,[],2);
                unit_found_ind2 = 1:size(error,1);
                mat_ind2 = sub2ind( size(temp_acc), unit_found_ind2, mapping2' );

                error(unit_found_ind2,n) = temp_error(mat_ind2);
                acc(unit_found_ind2,n) = temp_acc(mat_ind2);
                tpr(unit_found_ind2,n) = temp_tpr(mat_ind2);
                spc(unit_found_ind2,n) = temp_spc(mat_ind2);
            end
            
            meanSNR = this.getMeanUnitSNR(1);
            
            error_struct.total_error = total_error;
            error_struct.error = error;
            error_struct.acc = acc;
            error_struct.tpr = tpr;
            error_struct.spc = spc;
            error_struct.meanUnitSNR = meanSNR{1};
            
        end
        
        function error_struct = calculateErrorFromKnownSortCodeBYMETHOD(this, knownSortCode)
            error_struct = calculateErrorFromKnownSortCode(this, knownSortCode);
                        
            [~, method_names, indices_to_save] = optimalSortCodeBYMETHOD(this);
            
            error_struct.total_error = error_struct.total_error(indices_to_save);
            error_struct.error       = error_struct.error(:,indices_to_save);
            error_struct.acc         = error_struct.acc(:,indices_to_save);
            error_struct.tpr         = error_struct.tpr(:,indices_to_save);
            error_struct.spc         = error_struct.spc(:,indices_to_save);
            error_struct.methodNames = method_names;
        end
            
        function success = updateDataObj(this, str_var_name, var_value)
        
            if any(ismember(properties(this), str_var_name)),
                try
                    this.(deblank(str_var_name)) = var_value;
                    success = true;
                catch error_message
                    fprintf('%s \n', error_message.message)
                    success = false;
                end  
            end
        end       
        
        function plotPCA_WF(this, sortCode_ind)
        
            %% Constants
            XAxisLabel = 'Data Samples';
            YAxisLabel = 'Amplitude';

            PCA_Title = 'PCA';
            Waveform_Title = 'Waveforms';

            wfAxis = [];
            pcaAxis = [];
            h = [];
            
            wf = this.wf_snippets';
            sortCode = this.sortCodeList{sortCode_ind};
            sortCode = sortCode(:)';

            %% Code
            unitIDs = unique(sortCode);
            number_units = max(size(unitIDs));

            unit_num_wf = [];
            for unit_ind = 1:number_units,
                unit_num_wf(unit_ind) = sum(sortCode(:) == unitIDs(unit_ind));      
            end

            [unit_num_wf,ind]  = sort(unit_num_wf', 'descend');
            unitIDs = unitIDs(ind);


            %% Plot PCA -----------------------------------------------------------
            h = figure;
            pos = get(gcf, 'Position');
            set(gcf, 'Position', [100 pos(2) pos(3)*1.5 pos(4)/1.5])

            pcaAxis = subplot(1,2,1);
                cla
                grid on 
                hold all
                xlabel('PC 1');
                ylabel('PC 2');
                title(PCA_Title);
            
            [coeff{1}, score{1}] = princomp(wf);
            PCA_2D_point_1 = score{1}(:,1);
            PCA_2D_point_2 = score{1}(:,2);

            h_pca = zeros([sum(unit_num_wf), 1]);
            for unit_ind = 1:number_units,
                indices = sortCode == unitIDs(unit_ind);

                h_pca(indices) = plot( repmat(PCA_2D_point_1(indices),[2 1]), ...
                                       repmat(PCA_2D_point_2(indices),[2 1]), ...
                                       'linestyle', '*');     
            end

            %% Plot Waveforms -----------------------------------------------------
            wfAxis = subplot(1,2,2);
                cla
                hold on
                grid on
                ylabel(YAxisLabel);
                xlabel(XAxisLabel);
                title(Waveform_Title);
                
            wf_handles = {};
            legend_handles = [];
            unitLabels = [];

            for unit_ind = 1:number_units,
                indices = find(sortCode == unitIDs(unit_ind));
                unitIDs_inOrder = unique(unitIDs);

                wf_handles{unit_ind} = plot( wf(indices,:)', ...
                                             'Color', get(h_pca(indices(1)), 'color'));     

                legend_handles(unit_ind) = wf_handles{unit_ind}(1);                      

                strUnitId = num2str(unitIDs(unit_ind));
                strUnitNumWf = num2str(unit_num_wf(unit_ind));
                unitLabels{unit_ind} = sprintf('%s - %s', strUnitId, strUnitNumWf);
            end

            legend('off')
            axis tight

            % Update Legend
            legend(legend_handles, unitLabels{:}, 'Location', 'NorthEastOutside');
            
        end
        
    end
    
    %% Hidden Methods
    methods (Access = private)
       
        % Format inputs to object
        function this = formatInputs(this)
            
            this.number_waveforms = length(this.time_stamps);
            this.wf_continuous = this.wf_continuous(:);
            
            if size(this.wf_snippets,2) ~= this.number_waveforms,
                this.wf_snippets = this.wf_snippets';
                if size(this.wf_snippets,2) ~= this.number_waveforms,
                    error('Data Dimension MisMatch!');
                end
            end
            
            this.number_samples = size(this.wf_snippets,1);         
        end        
        
    end
    
end

