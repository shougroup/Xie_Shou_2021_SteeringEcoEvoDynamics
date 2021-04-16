clear
% uncomment one of the below 3 blocks to generate the corresponding figures

% % "No spike"
% phiM0 = 0.7;

% % "30%-H spike"
% spike_frac = 0.3;
% phiM0 = 0.7*(1-spike_frac);

% "30%-M spike"
spike_frac = 0.3;
phiM0 = 0.7*(1-spike_frac)+spike_frac;

% if fp = 0.13, phi_SS = 0.717
phi_SS = 0.717;
% Parameters shown in Table 1
gMmax = 0.7;
gHmax = 0.3;
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR = 1/3;
K_HR = 1/3;
K_MB = 100/3*2;
d1 = gMmax*5e-3;
d2 = gHmax*5e-3;

% N0: BM_0, R(0)=1, T0 is the maturation time
N0 = 100;
R0 = 1;
T0 = 17;

options=odeset('RelTol',1e-6,'abstol',1e-10);

fp = (0.1: 0.005: 0.16);
P = zeros(size(fp));
for j = 1:length(fp)
    para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp(j); K_MR; K_HR; K_MB; d1; d2];
    [T,X]=ode15s(@(t,x) MHDynamics(t,x,para),[0 T0],[N0*phiM0;N0*(1-phiM0);0;R0;0],options);
    if ~isreal(X)
        error('imaginary value')
    elseif nnz(X<-1e-5)>0
        error('negative value')
    end
    P(j) = X(end, 5);    
end
P_fit = polyfit(fp, P, 1);

figure(1)
plot(fp, P, 'm', 'linewidth', 2)
hold on
plot([0 1], [0 1]*P_fit(1)+P_fit(2), 'k:', 'linewidth', 2)
hold off
axis([0.1 0.16 100 800])
title (['slope = ', num2str(round(P_fit(1)))])
xlabel('$\overline{f_{P}}(0)$', 'interpreter', 'latex')
ylabel('P(T)')
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches',...
    'ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])
   