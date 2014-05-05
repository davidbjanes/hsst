function [weights,alpha,llk]=HGPMAPweights(nz,alpha,alpha_0)
% Sample from a HGP according to:
% M. Zhou and L. Carin, "Augment-and-conquer Negative Binomial Processes,"
% NIPS 2012
%
% dec 8/12/13
c=1;
p=.5;
minval=1e-6;
for iter=1:10
    Ntables=eCRT(nz,alpha);
    alpha=max(minval,(alpha_0+sum(Ntables,2)-1)./((c-log(1-p))));
end
weights=max(minval,(alpha+nz-1).*p);
llk=sum(alpha*log(p))-sum(gammaln(alpha))+(alpha-1)'*log(weights)...
    -p/(1-p)*sum(weights)+(alpha_0-1)*sum(log(alpha))-sum(alpha);


function t=eCRT(m,r)
t=zeros(size(m));
for n=1:numel(m)
    p=r(n)./(r(n):m(n)-1+r(n));
    t(n)=sum(p(:));
end