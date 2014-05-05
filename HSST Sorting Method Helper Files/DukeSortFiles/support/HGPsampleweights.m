function [weights,alpha,llk]=HGPsampleweights(nz,alpha,alpha_0)
% Sample from a HGP according to:
% M. Zhou and L. Carin, "Augment-and-conquer Negative Binomial Processes,"
% NIPS 2012
%
% dec 8/12/13
minval=1e-5;
c=1;
p=.5;
Ntables=CRTrnd(nz,alpha);
alpha=gamrnd(alpha_0+sum(Ntables,2),1./(c-log(1-p)));
alpha=max(alpha,minval);
weights=gamrnd(alpha+nz,p);
weights=max(weights,minval);
llk=sum(alpha*log(p))-sum(gammaln(alpha))+(alpha-1)'*log(weights)...
    -p/(1-p)*sum(weights)+(alpha_0-1)*sum(log(alpha))-sum(alpha);