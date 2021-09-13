% A: spike_frac = 0. For the parent cycle, uncomment line 7 and comment 8 and 9. For the offspring cycle,
% uncomment line 9 and comment 7 and 8.
% B: spike_frac = 0.3. For the parent cycle, uncomment line 7 and comment 8 and 9. For the offspring cycle,
% uncomment line 8 and comment 7 and 9.
clear
spike_frac = 0; 
% phiM0 = [0.65 0.95]*(1-spike_frac);
phiM0 = [0.32 0.36 0.41 0.47 0.51 0.56 0.62 0.66]; % offspring w/ spiking
% phiM0 = [0.5 0.56 0.6 0.67 0.74 0.79 0.86 0.93]; % offspring w/o spiking
fp = 0.13;
% Parameters shown in Table 1
gMmax = 0.7;
gHmax = 0.3;
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR = 1/3;
K_HR = 1/5;
K_MB = 100/3;
d1 = gMmax*5e-3;
d2 = gHmax*5e-3;

% N0: BM_0, R(0)=1, T0 is the maturation time
N0 = 100;
R0 = 1;
T0 = 17;

options=odeset('RelTol',1e-6,'abstol',1e-10);
para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; d1; d2];
% c = {'b', [0 0.6 0]};
c(1:8) = {'k'};
hold on
for j = 1:length(phiM0)
[T,X]=ode15s(@(t,x) MHDynamics(t,x,para),[0 T0],[N0*phiM0(j);N0*(1-phiM0(j));0;R0;0],options);
if ~isreal(X)
    error('imaginary value')
elseif nnz(X<-1e-5)>0
    error('negative value')
end
% 
plot(T, X(:, 1)./(X(:, 1)+X(:, 2)), 'color', c{j}, 'linewidth', 2)
% hold on
plot([0 T(end)], [phiM0(j) X(end, 1)./(X(end, 1)+X(end, 2))*(1-spike_frac)], 'o', 'color', c{j}, 'linewidth', 3)
end
plot([0 T0], [0.726 0.726], 'linewidth', 3, 'color', [1 0.5 0])
hold off
axis([0 T0 0 1])
xticks([])
yticks((0:0.2:1))
xlabel('t')
ylabel('\phi_{M}(t)')
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches',...
    'ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])
% print('HighSpike.svg', '-dsvg')   