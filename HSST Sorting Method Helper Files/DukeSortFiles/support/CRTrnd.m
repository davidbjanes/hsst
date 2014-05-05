function t=CRTrnd(m,r)
% Samples the number of tables used in a chinese restaurant process with
% total number of customers m and parameter r.
t=zeros(size(m));
for n=1:numel(m)
    p=r(n)./(r(n):m(n)-1+r(n));
    t(n)=sum(rand(m(n),1)<p(:));
end