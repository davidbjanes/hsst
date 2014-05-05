function llk=fitllk(x,d,s,lamb)
[P,N]=size(x);
xm=x-d*s;
clamb=chol(lamb);
llk=N*log(sum(diag(lamb)))-.5*norm(clamb*xm).^2;
