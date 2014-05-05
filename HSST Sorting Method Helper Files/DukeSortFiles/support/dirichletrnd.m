function pii=dirichletrnd(weights)
% Samples from a dirichlet distribution.
% Assumes weights are columnwise, so if there are multiple columns this
% will sample multiple Dirichlet distributions.
pii=gamrnd(weights,1);
pii=bsxfun(@rdivide,pii,sum(pii,1));