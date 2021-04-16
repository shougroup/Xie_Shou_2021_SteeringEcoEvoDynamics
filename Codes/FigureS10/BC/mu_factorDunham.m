function y=mu_factorDunham(n)
s1=0.05;
s2=0.067;
u=rand(n,1);
y=zeros(size(u));
idx=find(u<=s2/(s1+s2));
y(idx)=log((s1+s2)*u(idx)/s2)*s2;
idx=find(u>s2/(s1+s2));
y(idx)=-s1*log(1-(u(idx)-s2/(s1+s2))*(s1+s2)/s1);
y(y<-1)=-1;










