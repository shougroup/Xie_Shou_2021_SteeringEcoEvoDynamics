function y = mu_spec(n)
u=rand(n, 1);
y=zeros(size(u));

% % half of the mutations are null
% y(u<=0.5) = -1;
% 
% % half of the mutations are drawn from the distribution in Eq. xx
% idx = find(u > 0.5);
% if ~isempty(idx)
%     y(idx) = mu_factorDunham(length(idx));
% end

y(u <= 0.1) = -1;

idx = find(u > 0.9);
if ~isempty(idx)
    y(idx)=mu_factorDunham(length(idx));
end

