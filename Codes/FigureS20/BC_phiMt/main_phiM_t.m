clear
% B: no spiking
spike_frac = 0;
phiM0 = [0.5 0.85]; 

% % C: 30%-H spiking
% spike_frac = 0.3;
% phiM0 = [0.35 0.6]; % 30% spiking

fp = 0.18;
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
T0 = 17/2;

options=odeset('RelTol',1e-6,'abstol',1e-10);
para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; d1; d2];
c = {'b', [0 0.6 0]};
for j = 1:2
[T,X]=ode15s(@(t,x) MHDynamics(t,x,para),[0 T0],[N0*phiM0(j);N0*(1-phiM0(j));0;R0;0],options);
if ~isreal(X)
    error('imaginary value')
elseif nnz(X<-1e-5)>0
    error('negative value')
end
plot(T, X(:, 1)./(X(:, 1)+X(:, 2)), 'color', c{j}, 'linewidth', 2)
hold on
% plot([0 T(end)], [phiM0(j) X(end, 1)./(X(end, 1)+X(end, 2))*(1-spike_frac)], 'o', 'color', c{j}, 'linewidth', 3)
end
plot([0 T0], [0.68 0.68], ':', 'linewidth', 3, 'color', [255 127 42]/255)
if spike_frac > 0
    plot([0 T0], [0.68 0.68]*(1-spike_frac), ':', 'linewidth', 3, 'color', [0 0.5 0.5])
end
hold off
axis([0 T0 0 1])
xticks([])
% yticks((0.1:0.2:0.7))
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches',...
    'ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])
% print('HighSpike.svg', '-dsvg')   