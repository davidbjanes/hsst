function [map,llk]=DPDLGMM(x,s,params);
[N,P]=size(s);
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
lamb_0=eye(P)./nu_0;lamb_0inv=inv(lamb_0);lamb_0c=chol(lamb_0);lamb_0ci=chol(lamb_0inv);
% Save space for Gaussian parameters s~N(muu,Lamb^-1)
clus_muu=zeros(maxClusters,P);
clus_lamb=cell(maxClusters,1);
full_lamb=clus_lamb;
% Set mixture weights
mixture_weights=dirichletrnd(repmat(alpha_mixture,maxClusters,1));
maxllk=-inf;
dictionary=(s\x')';
dictionary=bsxfun(@rdivide,dictionary,sqrt(sum(dictionary.^2)));
s=(dictionary\x)';
signalLambda=params.signalLambda;
signalSigma=inv(signalLambda);
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
    %% sample a dictionary:
%     [dictionary]=sampleDictionary(x,dictionary,s',signalLambda);  
%     s=s';
    [dictionary]=sampleFullDictionary(x,dictionary,s',signalLambda);
    %% sample the new factor scores s
    s=sampleFactorScores(x,dictionary,s,signalLambda,clus_muu,clus_lamb,z);
    %% Collapse variables;
    full_muu=clus_muu*dictionary';
    dLd=dictionary'*signalLambda*dictionary;
    for clus=1:maxClusters
        tmp=signalLambda-signalLambda*dictionary*((clus_lamb{clus}+dLd)\...
            dictionary'*signalLambda);
        full_lamb{clus}=.5*(tmp+tmp');
    end
    %% sample which cluster each spike is from:
%     z=sampleGMM(s,mixture_weights,clus_muu,clus_lamb);
    z=sampleGMM(x',mixture_weights,full_muu,full_lamb);
    %% sample mixture weights:
    weights=alpha_mixture+full(sparse(z,ones(N,1),ones(N,1),maxClusters,1));
    mixture_weights=dirichletrnd(weights);
    %% sample Gaussian parameters:
    [clus_muu,clus_lamb]=updateAllNormalWisharts(s,z,lamb_0,lamb_0inv,lamb_0c,nu_0,kappa_0,muu_0,maxClusters);
    %% Store if necessary
    if (iter>burnin) & mod(iter-burnin,space)==0
        % Calculate log-likelihood
        collectioniter=(iter-burnin)/space;
        llk(collectioniter)=fitllk(x,dictionary,s',signalLambda)...
            +GMMloglikelihood(s,z,clus_muu,clus_lamb,mixture_weights);
        if llk(collectioniter)>maxllk
           maxllk=llk(collectioniter);
           map.llk=maxllk;
           map.z=z;
           map.dictionary=dictionary;
           map.s=s;
           map.clus_muu=clus_muu;
           map.clus_lamb=clus_lamb;
           map.mixture_weights=mixture_weights;
        end
    end
end

