
classdef MATLAB_KMeans < hsst.sorter
    
    properties (Constant)  
        sortMethodLabel = 'MATLAB-KM';
        defaultParameters = [1:6];
    end
    
    methods (Static)
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, ~, thresh, ...
                                                 ~, ~, ~, ...
                                                 sortparameter)
            
            % Handle Inputs
            if size(wf, 1) ~= length(ts),
                wf = wf';
                if size(wf, 1) ~= length(ts),
                    error('Time Stamp and Number of Waveform Mismatch');
                end
            end
            
            K = sortparameter;
            debug = false;
            showplots = false;
            V_thresh = thresh;
                        
            num_snips = size(wf, 1);
            num_samples = size(wf, 2);
                
            
%             %% Initalize U_K
%             snip_ind = datasample([1:num_snips], min(num_snips,1000), 'Replace',false);
% %             snip_ind = datasample([1:num_snips], max(min(num_snips,1000),num_snips*0.5), 'Replace',false);
% %             snip_ind = [1:num_snips];
%             est_mu_wf = wf(snip_ind,:);             
%             dist_matrix = squareform(pdist(est_mu_wf));
% 
%             val = max(max(dist_matrix));
%             [~, init_indices] = find(dist_matrix == val);
%             chosen_indices = 1:size(est_mu_wf,1) == init_indices(1);
% 
%             for k = 2:K,
%                 total_dist_sum = sum(dist_matrix(chosen_indices, :), 1);
%                 total_dist_sum(chosen_indices) = 0;
% 
%                 [~, next_index] = max(total_dist_sum);
%                 chosen_indices(next_index) = 1;
%             end
% 
%             [COEFF, score] = princomp(est_mu_wf);
%             mean_waveform = mean(est_mu_wf);
%     
%             low_d = score(chosen_indices, :);
%             low_d(:,3:end) = 0;
%             projected = bsxfun(@plus, low_d/COEFF, mean_waveform);
%             start_mu = projected;
            
                        
            %% Sort N times, pick lowest sort 
            
            warning('off', 'all');
            completed_successfully = [];
            completed_successfully(1) = 0;
            N = 5;
%             repmat(start_mu,[1 1 N])
            while all(~completed_successfully) & length(completed_successfully) < 5,
            try 
                [sortCode, ~, ~] = kmeans(wf, K, ...
                                          'start', 'cluster', ...
                                          'Replicates', N);
                completed_successfully(end+1) = 1;
            catch
                fprintf('KMeans Failure, will try again. ');
                completed_successfully(end+1) = 0;
            end
            end
            warning('on', 'all');
            
            if all(~completed_successfully),
                sortCode = ones([1 num_snips]);
            end
            
            if debug, 
                plotWaveforms(snips, sortCode, V_thresh, [])
            end
                
            property = [];
            sort_ids = sortCode;
            
            
        end
        
    end
    
end
