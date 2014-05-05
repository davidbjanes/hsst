function [z,Likelihood]=sampleGMM(s,mixture_weights,muu,lamb)
% Chooses which mixture component is associated with each data point and
% returns a entry z
[N,P]=size(s);
maxClus=numel(mixture_weights);
logLikelihood=zeros(maxClus,N);
for clus=1:maxClus
    lambchol=chol(lamb{clus});
    logLikelihood(clus,:)=log(mixture_weights(clus))+...
        sum(diag(lambchol))-.5*sum((lambchol*bsxfun(@minus,s,muu(clus,:))').^2);
end
logLikelihood=bsxfun(@minus,logLikelihood,max(logLikelihood));
Likelihood=exp(logLikelihood);
Likelihood=bsxfun(@rdivide,Likelihood,sum(Likelihood));
% z=Discrete(Likelihood);
z=sample_vector(Likelihood);


function z=Discrete(P)
% Samples a discrete distribution from (un)normalized probabilities P
N=size(P,2);
Pcum=cumsum(P);
z=sum(bsxfun(@gt,P,(rand(1,N).*Pcum(end,:))))+1;
if max(z)>20
   1; 
end
1;


function x = sample_vector(p,col)
% SAMPLE_VECTOR    
% From the LIGHTSPEED project: http://research.microsoft.com/en-us/um/people/minka/software/lightspeed/
%
% Sample from multiple categorical distributions.
% X = SAMPLE_VECTOR(P) returns a row vector of cols(P) integers, where
% X(I) = SAMPLE(P(:,I)).
%
% X = SAMPLE_VECTOR(P,COL) returns a row vector of length(COL) integers, where
% X(I) = SAMPLE(P(:,COL(I))).
% This is equivalent to SAMPLE_VECTOR(P(:,COL)), but faster.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

cdf = cumsum(p);
if nargin >= 2
  cdf = cdf(:,col);
end
if any(cdf(end,:) <= 0)
  error('distribution is all zeros');
end
[rows cols] = size(cdf);
u = rand(1,cols).*cdf(end,:);
x = sum(cdf < repmat(u,rows,1),1) + 1;