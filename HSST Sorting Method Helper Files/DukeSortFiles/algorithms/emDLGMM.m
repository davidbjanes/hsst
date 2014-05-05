function [map,llk]=emDLGMM(x,params,map);
[N,P]=size(map.s);
T=size(x,1);
%% Parameters
burnin=params.burnin;
collection=params.collect;
space=params.space;
maxiter=burnin+collection*space;
maxClusters = params.maxClusters;
initClusters = params.initClusters;
alpha_mixture=params.mixture_hyper;
% Atomic parameters (Gaussian-Wishart)
muu_0=zeros(P,1);
kappa_0=1;
nu_0=2*P;
lamb_0=1*eye(P);lamb_0inv=inv(lamb_0);lamb_0c=chol(lamb_0);lamb_0ci=chol(lamb_0inv);
maxllk=-inf;
signalLambda=params.signalLambda;
maxllk=map.llk;
HGPllk=0;
clus_muu=map.clus_muu;
clus_lamb=map.clus_lamb;
mixture_weights=map.mixture_weights;
dictionary=map.dictionary;
s=map.s;
%% MCMC sampler
for iter=1:10
    %% sample which cluster each spike is from:
    [~,zdistribution]=sampleGMM(s,mixture_weights,clus_muu,clus_lamb);
    [~,z]=max(zdistribution);
    %% sample a new dictionary:
    [dictionary]=sampleDictionary(x,dictionary,s',signalLambda,1);
    %% sample the new factor scores s
    s=sampleFactorScores(x,dictionary,s,signalLambda,clus_muu,clus_lamb,z,1);
    %% sample mixture weights:
    weights=alpha_mixture+sum(zdistribution,2);
    mixture_weights=max(weights-1,0);
    mixture_weights=mixture_weights./sum(mixture_weights);
    %% sample Gaussian parameters:
    [clus_muu,clus_lamb]=MAPAllNormalWisharts(s,zdistribution,lamb_0,lamb_0inv,lamb_0c,nu_0,kappa_0,muu_0);
    %% Store if necessary
    llk(iter)=fitllk(x,dictionary,s',signalLambda)...
        +GMMloglikelihood(s,z,clus_muu,clus_lamb,mixture_weights);
        maxllk=llk;
        map.z=z;
        map.dictionary=dictionary;
        map.s=s;
        map.clus_muu=clus_muu;
        map.clus_lamb=clus_lamb;
        map.mixture_weights=mixture_weights;
        map.alpha_mixture=alpha_mixture;
end

