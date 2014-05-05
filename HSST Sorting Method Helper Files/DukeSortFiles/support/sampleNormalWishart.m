function [muu,lamb]=sampleNormalWishart(x,muu_0,kappa_0,nu_0,lamb_0inv)
% function [muu,lamb]=sampleNormalWishart(x,muu_0,kappa_0,nu_0,lamb0inv)
%
% Samples from a normal-wishart distribution with prior
% NW(muu_0,kappa_0,nu_0,lamb0) and observed data x is p by n;
[P,N]=size(x);
xmean=mean(x,2);
kappa=kappa_0+N;
muup=(N*xmean+kappa_0*muu_0)/kappa;
nu=nu_0+N;
if N>1
    lambpinv=lamb_0inv+N*cov(x')+kappa_0*N/(kappa_0+N)*(xmean-muu_0)*(xmean-muu_0)';
else
    lambpinv=lamb_0inv+kappa_0*N/(kappa_0+N)*(xmean-muu_0)*(xmean-muu_0)';
end
lambp=inv(lambpinv);lambp=.5*(lambp+lambp');
lambpc=chol(lambp);
lamb=wishrnd(lambp,nu,lambpc);
muu=muup+chol(lamb)\randn(P,1)./sqrt(kappa);
1;