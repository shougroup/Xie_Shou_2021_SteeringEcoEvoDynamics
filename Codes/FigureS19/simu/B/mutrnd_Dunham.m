% mutation spectrum obtained from Dunham lab

function [fp_mut] = mutrnd_Dunham(params,fp_background)
sp0=params(1);
sn0=params(2);
g=params(3);
sp=sp0./(1+g*(fp_background-0.13)/0.13);
sn=sn0*(1+g*(fp_background-0.13)/0.13);
nc=sp+sn.*(1-exp(-1./sn));
divalue=sn.*(1-exp(-1./sn))./nc;

u=rand(size(fp_background));
delta_fp=zeros(size(fp_background));

idx=find(u<=divalue);
if ~isempty(idx)
    delta_fp(idx)=sn(idx).*log(u(idx).*(sp(idx)+sn(idx).*(1-exp(-1./sn(idx))))./sn(idx)+exp(-1./sn(idx)));
end
idx=find(u>divalue);
if ~isempty(idx)
    delta_fp(idx)=-sp(idx).*log((sp(idx)+sn(idx).*(1-exp(-1./sn(idx))))./sp(idx).*(1-u(idx)));
end
fp_mut = fp_background .* (1 + delta_fp);
fp_mut(fp_mut > 1) = 1;
end