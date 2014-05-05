function s=sampleFactorScores(x,d,s,lamb,clus_muu,clus_lamb,z,map)
% Assumes a prior of d_k~N(0,1/P I_P)
m=1;
if nargin>7;if map>0;
        m=0;
end;end;
[P,N]=size(x);
K=size(d,2);
maxClus=size(clus_muu,1);
DLD=d'*lamb*d;
for clus=1:maxClus
    ndx=z==clus;
    n=sum(ndx);
    if n==0
        continue
    end
    prec=DLD+clus_lamb{clus};
    prec=prec+eye(K)*1e-4*max(prec(:));
    precmuu=bsxfun(@plus,(d'*lamb)*x(:,ndx),clus_lamb{clus}*clus_muu(clus,:)');
    cprec=chol(prec);
    s(ndx,:)=(cprec\(m*randn(K,n)+cprec'\precmuu))';
    if max(s(:))>100
       1; 
    end
end
1;
