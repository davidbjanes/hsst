classdef DukeSort < hsst.sorter
    
    properties (Constant)
        sortMethodLabel = 'DukeSort';
        defaultParameters = [1e-3 1e-2 0.01 0.1 1 10 100];
    end
        
    methods (Static)
        
        function [sort_ids, wf, ts, property] = sort(wf, ts, Fs, thresh, ...
                                                     align_sample, dur, raw_data, ...
                                                     sortparameter)
            
            verbose = true;
            showplots = false;
            
            
            %% Parameters
            % basic parameters
            params.maxClusters  = 20;  % Upper bound on number of clusters
            params.initClusters = 10;  % How many clusters to use for initialization
            params.burnin       = 100; % Number of burn-in samples
            params.collect      = 500; % Number of collection samples
            params.space        = 1;   % Number of samples until next collection

            % mixture model parameters:
%             params.mixture_hyper    = 1; % lower is fewer clusters, higher is more clusters
            params.mixture_hyper    = sortparameter; 
            
            % PCA/Dictionary Parameters:
            params.numDictionary    = 3; % How many dictionary/ PCA terms to use.

            % Signal Noise Parameters, only used in dictionary learning models
            params.noiseVariance    = 1/30; % Roughly correct for 3sigma threshold detection
            params.lag1correlation  = 0.9;  % Noise lag-1 autocorrelation
            
            model = 1;
            
            try
                
                %% Load Data
                if length(ts) == size(wf, 2),
                    x = wf;
                else
                    x = wf';
                end

                xn = sqrt(size(x,2))./(norm(x)).*x;
                [u,e,v] = svd(x,'econ');
                s = v(:,1:params.numDictionary);
                s = s.*sqrt(size(x,2));

                [P, N] = size(x);
                tmp = (ones(P,1))*(1+params.lag1correlation^2);
                tmp(1) = 1;
                tmp(end) = 1;
                params.signalLambda = 1./params.noiseVariance*(params.lag1correlation*....
                    (diag(ones(P-1,1),1)+diag(ones(P-1,1),-1))+diag(tmp));

                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Model 1 is simply a GMM where we learn the number
                % of mixture components on the SVD decomposition. This
                % is called the Dirichlet Process Gaussian Mixture
                % model.  See:
                %
                % D. Carlson, J. Vogelstein, Q.Wu, W. Lian, M. Zhou,
                % C.R. Stoetzner, D. Kipke, D. Weber, D. Dunson
                % and L. Carin, "Sorting Electrophysiological Data
                % via Dictionary Learning  & Mixture Modeling", to
                % appear in IEEE TBME, 2013
                %
                % "A non-parametric Bayesian approach to spike sorting,"
                % F. Wood, S. Goldwater, and M. J. Black, EMBS, 2006.
                %
                % Initial fit is done with a MCMC chain.  The post-
                % processing is done with an EM algorithm.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if model == 1;
                    [DPGMMmap, logLikelihood] = DPGMM(s, params);

    %                 figure(1)
    %                 drawPCs(s,DPGMMmap.z);

    %                 figure(2)
    %                 drawWaveforms(x,DPGMMmap.z);

                    %  Post-processing for DPGMM
    %                 addpath('emgm')
                    z = DPGMMmap.z;
                    n = numel(z);
                    nz = full(sparse(z,ones(n,1),ones(n,1)));
                    [nzs,rendx] = sort(nz,'descend');
                    numClus = find(nzs==0,1,'first')-1;
                    z2 = z;
                    for c = 1:numClus
                        z2(z==rendx(c)) = c;
                    end
                    [zDPGMMPost] = emgm(s',z2(:)');
                    sortCode =  zDPGMMPost;
    %                 figure(3)
    %                 drawPCs(s,zDPGMMPost);
    %                 
    %                 figure(4)
    %                 drawWaveforms(x,zDPGMMPost);
                end


                %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Model 2 is simply a GMM where we learn the number
                % of mixture components in a novel way via a
                % heirarchical Gamma Process on the SVD decomposition.
                % This is called the Hierarchical Gamma Process
                % Gaussian Mixture model.  See:
                %
                % D. Carlson, J. Vogelstein, Q.Wu, W. Lian, M. Zhou,
                % C.R. Stoetzner, D. Kipke, D. Weber, D. Dunson
                % and L. Carin, "Sorting Electrophysiological Data
                % via Dictionary Learning  & Mixture Modeling", to
                % appear in IEEE TBME, 2013
                %
                % Initial fit is done with a MCMC chain.  The post-
                % processing is done with an EM algorithm.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if model==2;
                    [HGPmap, logLikelihood] = HGPGMM(s,params);

    %                 figure(1)
    %                 drawPCs(s,HGPmap.z);
    %                 
    %                 figure(2)
    %                 drawWaveforms(x,HGPmap.z);

                    [HGPpost, ellk] = emHGP(s,params,HGPmap);
                    zHGPpost = HGPpost.z;

                    sortCode =  zHGPpost;
    %                 figure(3)
    %                 drawPCs(s,HGPpost.z);
    %                 
    %                 figure(4)
    %                 drawWaveforms(x,HGPpost.z);
                end
                
            catch error_message
            
                disp(error_message)
                sortCode = zeros([1 length(ts)]);
            
            end
            
            
            %% Output
            sort_ids = sortCode;
            wf = wf;
            ts = ts;
            property = [];
            
        end
        
    end    
    
end