function [result] = fastbinorv(N,p)
if isempty(N)
    error('you are calling fastbniorv on an empty vector. This is a waste of everyone''s time.');
end

[~,n] = size(N);
if n > 1
    error('N must be a column vector')
end

thresh_b2p = 30;
thresh_p2n = 800;
normcond = @(N,p) N > 9 * (1-p) ./ p & N > 9 * p ./ (1-p);
result = zeros(length(N),1);
if length(p) == 1
    p = ones(length(N),1) * p;
    result = N;
end
if and(min(N) > thresh_p2n , normcond(N,p))
    result = round(normrnd(N .* p, sqrt(N .* p .* (1 - p))));
else
    
    bino_ind = find(N > 0 & N < thresh_b2p);
    pois_ind = find(N >= thresh_b2p & (N < thresh_p2n | ~normcond(N,p)));
    norm_ind = find(N >= thresh_p2n & normcond(N,p));
    if ~isempty(bino_ind)
        randy = randi(10^8);
        result(bino_ind) = (mybinornd_mex(int32(randy),(N(bino_ind))', (p(bino_ind))'))';
    end
    result(pois_ind) = poissrnd(N(pois_ind) .* p(pois_ind));
    result(norm_ind) = round(normrnd(N(norm_ind) .* p(norm_ind), sqrt(N(norm_ind) .*...
        p(norm_ind) .* (1 - p(norm_ind)))));
end
if sum(sum(result<0 + result>N)) > 0
    error('hi')
end
result(result<0) = result(result<0)* 0;
result(result>N) = result(result>N) * 0 + N(result>N);
end