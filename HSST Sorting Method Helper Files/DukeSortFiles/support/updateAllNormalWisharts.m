function [muu,lamb]=updateAllNormalWisharts(s,z,lamb_0,lamb_0inv,lamb_0c,nu_0,kappa_0,muu_0,maxClusters)
[N,P]=size(s);
muu=zeros(maxClusters,P);
lamb=cell(maxClusters,1);
for clus=1:maxClusters
    sclus=s(z==clus,:);
    if numel(sclus)==0 % no data, sample from prior
        lamb{clus}=wishrnd(lamb_0,nu_0,lamb_0c);
        muu(clus,:)=muu_0(:)'+(chol(lamb{clus})\randn(P,1))'./sqrt(kappa_0);
    else
        [muu(clus,:),lamb{clus}]=...
            sampleNormalWishart(sclus',muu_0,kappa_0,nu_0,lamb_0inv);
    end
end