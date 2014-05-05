function [map,llk]=HGPGMM(s,params)
[N,P]=size(s);
%% Parameters
burnin=params.burnin;
collection=params.collect;
space=params.space;
maxiter=burnin+collection*space;
maxClusters = params.maxClusters;
initClusters = params.initClusters;
alpha_mixture_0=params.mixture_hyper;
alpha_mixture=gamrnd(alpha_mixture_0,1,maxClusters,1);
% Atomic parameters (Gaussian-Wishart)
muu_0=zeros(P,1);
kappa_0=1;
nu_0=2*P;
lamb_0=eye(P)./nu_0;lamb_0inv=inv(lamb_0);lamb_0c=chol(lamb_0);lamb_0ci=chol(lamb_0inv);
% Save space for Gaussian parameters s~N(muu,Lamb^-1)
clus_muu=zeros(maxClusters,P);
clus_lamb=cell(maxClusters,1);
% Set mixture weights
mixture_weights=dirichletrnd(alpha_mixture);
maxllk=-inf;
%% Variable Initialization:
z = kmeans(s,initClusters);
for clus=1:maxClusters
    if clus<=initClusters
        clus_muu(clus,:)=mean(s(z==clus,:));
        clus_lamb{clus}=inv(cov(s(z==clus,:)));
    else
        clus_muu(clus,:)=randn(1,P);
        clus_lamb{clus}=wishrnd(lamb_0,nu_0,lamb_0ci);
    end
end
%% MCMC sampler
for iter=1:maxiter
    %% sample which cluster each spike is from:
    z=sampleGMM(s,mixture_weights,clus_muu,clus_lamb);
    %% sample mixture weights:
    counts=full(sparse(z,ones(N,1),ones(N,1),maxClusters,1));
    [mixture_weights,alpha_mixture,HGPllk]=HGPsampleweights(counts,alpha_mixture,alpha_mixture_0);
    mixture_weights=mixture_weights./sum(mixture_weights);
    %% sample Gaussian parameters:
    [clus_muu,clus_lamb]=updateAllNormalWisharts(s,z,lamb_0,lamb_0inv,lamb_0c,nu_0,kappa_0,muu_0,maxClusters);
    %% Store if necessary
    if (iter>burnin) & mod(iter-burnin,space)==0
        % Calculate log-likelihood
        collectioniter=(iter-burnin)/space;
        llk(collectioniter)=...
            GMMloglikelihood(s,z,clus_muu,clus_lamb,mixture_weights)+...
            HGPllk;
        if llk(collectioniter)>maxllk
           maxllk=llk(collectioniter);
           map.llk=maxllk;
           map.z=z;
           map.clus_muu=clus_muu;
           map.clus_lamb=clus_lamb;
           map.mixture_weights=mixture_weights;
           map.alpha_mixture=alpha_mixture;
        end
    end
end

