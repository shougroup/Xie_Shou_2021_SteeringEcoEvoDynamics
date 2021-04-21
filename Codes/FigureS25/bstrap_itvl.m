function [lb, ub] = bstrap_itvl(x,y,f,n_bstraps,q)
dist = zeros(n_bstraps,1);
for i = 1 : n_bstraps
    idx = randsample(length(x), length(x), true);
    dist(i) = f(x(idx), y(idx));
end
lb = quantile(dist, q);
ub = quantile(dist, 1-q);
