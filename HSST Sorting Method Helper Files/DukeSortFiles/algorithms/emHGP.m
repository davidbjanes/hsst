function [map,llk]=emHGP(s,params,map)
[N,P]=size(s);
%% Parameters
burnin=params.burnin;
collection=params.collect;
space=params.space;
maxiter=burnin+collection*space;
maxClusters = params.maxClusters;
alpha_mixture_0=params.mixture_hyper;
alpha_mixture=gamrnd(alpha_mixture_0,1,maxClusters,1);
% Atomic parameters (Gaussian-Wishart)
muu_0=zeros(P,1);
kappa_0=1;
nu_0=2*P;
lamb_0=1./nu_0*eye(P);lamb_0inv=inv(lamb_0);lamb_0c=chol(lamb_0);lamb_0ci=chol(lamb_0inv);
% Save space for Gaussian parameters s~N(muu,Lamb^-1)
clus_muu=map.clus_muu;
clus_lamb=map.clus_lamb;
% Set mixture weights
mixture_weights=map.mixture_weights;
maxllk=map.llk;
HGPllk=0;
%% EM optimization
for iter=1:30
    %% sample which cluster each spike is from:
    [~,zdistribution]=sampleGMM(s,mixture_weights,clus_muu,clus_lamb);
    [~,z]=max(zdistribution);
    %% sample mixture weights:
    counts=sum(zdistribution,2);
    [mixture_weights,alpha_mixture,HGPllk]=HGPMAPweights(counts,alpha_mixture,alpha_mixture_0);
    mixture_weights=mixture_weights./sum(mixture_weights);
    %% sample Gaussian parameters:
    [clus_muu,clus_lamb]=MAPAllNormalWisharts(s,zdistribution,lamb_0,lamb_0inv,lamb_0c,nu_0,kappa_0,muu_0);
    %% Store if necessary
    % Calculate log-likelihood
    llk(iter)=...
        GMMloglikelihood(s,z,clus_muu,clus_lamb,mixture_weights)+...
        HGPllk;
%     if llk(iter)>maxllk
        maxllk=llk;
        map.z=z;
        map.clus_muu=clus_muu;
        map.clus_lamb=clus_lamb;
        map.mixture_weights=mixture_weights;
        map.alpha_mixture=alpha_mixture;
%     end
end

