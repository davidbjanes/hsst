
classdef MATLAB_GMM_PCA < hsst.sorter
    
    properties (Constant)
        sortMethodLabel = 'MATLAB-GMM-PCA';
        defaultParameters = [1:6];
    end
    
    methods (Static)      
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, ~, ~, ...
                                                 ~, ~, ~, sortparameter)
            
            % Handle Inputs
            if size(wf, 1) ~= length(ts),
                wf = wf';
                if size(wf, 1) ~= length(ts),
                    error('Time Stamp and Number of Waveform Mismatch');
                end
            end
            
            K = sortparameter;
            
            num_snips = size(wf, 1);
            num_samples = size(wf, 2);

            
            %% PCA Decomp
            
            [coeff{1}, score{1}] = pca(wf);
            PCA_2D_point = score{1}(:,1:2);
            PCA_3D_point = score{1}(:,1:3);
            PCA_10D_point = score{1}(:,1:10);
            
            wf = PCA_3D_point;
            
            %% Sort
            % Gaussian Mixture Model of correlation histogram (Non
            % Deterministic : run 10 times and pick lowest NlogL value
            warning('off', 'all')
            N = 10;
            gmm = gmdistribution.fit(wf, K, ...
                                     'CovType', 'diagonal',...
                                     'Regularize', 10e-5, ...
                                     'Replicates', N);
            warning('on', 'all')
                                        
            wf_pdf = zeros([K num_snips]);
            for k = 1:K,
                wf_pdf(k,:) = mvnpdf( wf, gmm.mu(k,:), gmm.Sigma(:,:,k) );
            end
            
            [~, sortCode] = max(wf_pdf, [], 1);                           
            
            property = gmm.NlogL;
            sort_ids = sortCode;
            
        end
        
    end
    
end
