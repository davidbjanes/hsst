function [d,s]=sampleFullDictionary(x,d,s,lamb,map)
% Assumes a prior of d_k~N(0,1/P I_P)
% Rescales dictionary to 1.
m=1;
if nargin>4;if map>0;
        m=0;
end;end;
[P,N]=size(x);
K=size(d,2);
ss=s*s';
precmuu=lamb*(x*s');
precmuu=precmuu(:);
prec=P*eye(K*P)+kron(lamb,ss);
cprec=chol(prec);
dt=cprec\(m*rand(P*K,1)+cprec'\precmuu);
d=reshape(dt,P,K);
d=bsxfun(@rdivide,d,sqrt(sum(d.^2)));
1;
% xm=x-d*s;
% for k=1:K
%     sk=s(k,:);
%     dk=d(:,k);
%     xmk=xm+dk*sk;
%     lambp=P*eye(P)+lamb*sum(sk.^2);
%     tmp=lamb*(xmk*sk');
%     chollamb=chol(lambp);
%     dk=chollamb\(m*rand(P,1)+chollamb'\tmp);
%     scale=norm(dk);
%     dk=dk./scale;
% %     sk=sk.*scale;
%     xm=xm-dk*sk;
%     d(:,k)=dk;
%     s(k,:)=sk;
% end
% 1;

