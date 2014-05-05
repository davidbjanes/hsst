function likeli = LogNormal(X,mu,J)
[p,n] = size(X); Xmu = X - mu*ones(1,n); 
likeli = -0.5*p*log(2*pi) + sum(log(diag(chol(J)))) - 0.5*sum(Xmu.*(J*Xmu),1);
% likeli=logsumexp(likeli)-log(n);
return