clear
addpath('data')
addpath('support')
addpath('algorithms')
addpath('OPASS')
dataset = 2;
model = 2;

%% Parameter Sets:
% basic parameters
params.maxClusters  = 20;  % Upper bound on number of clusters
params.initClusters = 10;  % How many clusters to use for initialization
params.burnin       = 100; % Number of burn-in samples
params.collect      = 500; % Number of collection samples
params.space        = 1;   % Number of samples until next collection

% mixture model parameters:
params.mixture_hyper    = 1; % lower is fewer clusters, higher is more clusters

% PCA/Dictionary Parameters:
params.numDictionary    = 3; % How many dictionary/ PCA terms to use.

% Signal Noise Parameters, only used in dictionary learning models
params.noiseVariance    = 1/30; % Roughly correct for 3sigma threshold detection
params.lag1correlation  = 0.9;  % Noise lag-1 autocorrelation

%% load some data
if dataset==1
    % hc1
    load HarrisData
    % try a single channel:
    x=squeeze(ec_spikes(1,:,:));
    xn=sqrt(size(x,2))./(norm(x)).*x;
    [u,e,v]=svd(x,'econ');
    s=v(:,1:params.numDictionary);
    s=s.*sqrt(size(x,2));
elseif dataset==2
    % 4 cluster test set based on pitt data
    load clusters
    x=Spike;
    xn=sqrt(size(x,2))./(norm(x)).*x;
    [u,e,v]=svd(x,'econ');
    s=v(:,1:params.numDictionary);
    s=s.*sqrt(size(x,2));
end
[P,N]=size(x);
tmp=(ones(P,1))*(1+params.lag1correlation^2);tmp(1)=1;tmp(end)=1;
params.signalLambda=1./params.noiseVariance*(params.lag1correlation*....
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if model==1;
    %%
    [DPGMMmap,logLikelihood]=DPGMM(s,params);
    %%
    figure(1)
    drawPCs(s,DPGMMmap.z);
    %
    figure(2)
    drawWaveforms(x,DPGMMmap.z);
    %%  Post-processing for DPGMM
    addpath('emgm')
    z=DPGMMmap.z;n=numel(z);
    nz=full(sparse(z,ones(n,1),ones(n,1)));
    [nzs,rendx]=sort(nz,'descend');
    numClus=find(nzs==0,1,'first')-1;
    z2=z;
    for c=1:numClus
        z2(z==rendx(c))=c;
    end
    [zDPGMMPost]=emgm(s',z2(:)');
    %%
    figure(3)
    drawPCs(s,zDPGMMPost);
    %
    figure(4)
    drawWaveforms(x,zDPGMMPost);
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
    %%
    [HGPmap,logLikelihood]=HGPGMM(s,params);
    %%
    figure(1)
    drawPCs(s,HGPmap.z);
    %
    figure(2)
    drawWaveforms(x,HGPmap.z);
    %%
    [HGPpost,ellk]=emHGP(s,params,HGPmap);
    %%
    figure(3)
    drawPCs(s,HGPpost.z);
    %
    figure(4)
    drawWaveforms(x,HGPpost.z);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Still editing...
%
% Model 3 is a low-rank GMM where we learn the number
% of mixture components and we infer the low-rank
% dictionary in lieu of using predefined features. This
% model is called Dirichlet Process-Dictionary Learning
% Gaussian Mixture Model.
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
if model==3
    %%
    [DLmap,logLikelihood]=DPDLGMM(xn,s,params);
%     x=xn;
%     DPDLGMM
    %%
    figure(1)
    drawPCs(s,DLmap.z);
%     %%
%     figure(1);
%     drawPCs(DLmap.s,DLmap.z);
    %%
    figure(2)
    drawWaveforms(x,DLmap.z);
    %%
    [DLpost,ellk]=emDLGMM(xn,params,DLmap);
    %%
    figure(3)
    drawPCs(s,DLpost.z);
    %
    figure(4)
    drawWaveforms(x,DLpost.z);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 4 implements a real-time streaming data
% algorithm that clusters in one pass over the data.  
% Results will not necessarily be as good as the off-
% line algorithms, but can be implemented in systems
% for real time control.
%
% This model is called fake-OPASS. (OPASS and OPASS-a,
% provide integrated spike sorting and detection. 
% OPASS-a allows spike waveforms to evolve over time.
%  
% Article currently under review.
% D. Carlson, J. Vogelstein, V. Rao, and L. Carin.  
% "Real-Time Inference for a Gamma Process Model of 
% Neural Spiking," under review at NIPS 2013.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if model==4
     %%
    params.alph=1e-2; % Chinese restaurant process parameter.
    params.kappa_0=.05; 
    params.nu_0=1;
    K=params.numDictionary;
    params.Phi_0=.05*eye(K);
    [z,ngam,muu,Lam,nu,kappa,Phi]=fake_opass(s',params);
    
    %% Plot non-trivial clusters;
    figure(1)
    drawPCs(s,z);
        %%
    figure(2)
    drawWaveforms(x,z);
    %%
        %%  Post-processing for DPGMM
    addpath('emgm')
    z;n=numel(z);
    nz=full(sparse(z,ones(n,1),ones(n,1)));
    [nzs,rendx]=sort(nz,'descend');
    numClus=find(nzs==0,1,'first')-1;
    z2=z;
    for c=1:numClus
        z2(z==rendx(c))=c;
    end
    [zfOPASSp]=emgm(s',z2(:)');
    %%
    figure(3)
    drawPCs(s,zfOPASSp);
    %
    figure(4)
    drawWaveforms(x,zfOPASSp);
end





