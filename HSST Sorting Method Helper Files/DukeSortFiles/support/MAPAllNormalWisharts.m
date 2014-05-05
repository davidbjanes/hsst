function [muu,lamb]=MAPAllNormalWisharts(s,z,lamb_0,lamb_0inv,lamb_0c,nu_0,kappa_0,muu_0);
maxClusters=size(z,1);
[~,zm]=max(z);
zs=sum(z,2);

[N,P]=size(s);
muu=zeros(maxClusters,P);
lamb=cell(maxClusters,1);
for clus=1:maxClusters
%     sclus=s(zm==clus,:);
    if sum(zm==clus)==0 % no data, sample from prior
        lamb{clus}=lamb_0;
        muu(clus,:)=muu_0;
    else
        [muu(clus,:),lamb{clus}]=...
            MAPNormalWishart(z(clus,:),s',muu_0,kappa_0,nu_0,lamb_0inv);
    end
end

function [muu,lamb]=MAPNormalWishart(z,x,muu_0,kappa_0,nu_0,lamb_0inv)
% function [muu,lamb]=sampleNormalWishart(x,muu_0,kappa_0,nu_0,lamb0inv)
%
% Gets the MAP from a normal-wishart distribution with prior
% NW(muu_0,kappa_0,nu_0,lamb0) and observed data x is p by n;
sz=sum(z);
[P,N]=size(x);
xmean=sum(bsxfun(@times,x,z(:)'),2)./sz;
kappa=kappa_0+sz;
muup=(sz*xmean+kappa_0*muu_0)./kappa;
nu=nu_0+sz;
xm=bsxfun(@minus,x,xmean);
xms=bsxfun(@times,x,z(:)');
covx=xms*xm'./sz;
lambpinv=lamb_0inv+sz*covx+kappa_0*sz/(kappa_0+sz)*(xmean-muu_0)*(xmean-muu_0)';
lambp=inv(lambpinv);lambp=.5*(lambp+lambp');
lambpc=chol(lambp);
% lamb=wishrnd(lambp,nu,lambpc);
lamb=lambp*nu;
muu=muup;
1;