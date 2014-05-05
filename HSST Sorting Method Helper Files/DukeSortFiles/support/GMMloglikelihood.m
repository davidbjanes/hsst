function llk=GMMloglikelihood(s,z,muu,lamb,mixture_weights);
[N,P]=size(s);
maxClus=numel(mixture_weights);
logLikelihood=zeros(1,N);
for clus=1:maxClus
    ndx=find(z==clus);
    if numel(ndx)==0
        continue
    end
    lambchol=chol(lamb{clus});
    logLikelihood(1,ndx)=...
        sum(diag(lambchol))-.5*sum((lambchol*bsxfun(@minus,s(ndx,:),muu(clus,:))').^2);
end
nz=sparse(z,ones(N,1),ones(N,1),maxClus,1);
llk=sum(logLikelihood)+nz'*mixture_weights;
