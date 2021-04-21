clear
spike_frac = 0.3;
fp = 0.13;
% Parameters obtained from the simulation
gMmax = 0.7;
gHmax = 0.3;
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR = 1;
K_HR = 0.56;
K_MB = 85.7;
d1 = gMmax*5e-3;
d2 = gHmax*5e-3;

% N0: BM_0, R(0)=1, T0 is the maturation time
N0 = 100;
R0 = 1;
T0 = 17;

options=odeset('RelTol',1e-6,'abstol',1e-10);
para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; d1; d2];
% figure('color', 'none')
phiM0 = [0.5 0.7]*(1-spike_frac);
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
    plot(0, phiM0(j) , 'x', 'color', c{j}, 'linewidth', 2, 'markersize', 16)
    plot( T(end), X(end, 1)./(X(end, 1)+X(end, 2))*(1-spike_frac), 'o', 'color', c{j}, 'linewidth', 3)
end
plot([0 T0], [0.717 0.717], ':', 'linewidth', 3, 'color', [0.5 0.5 0])
hold off
axis([0 T0 0.1 0.8])
xticks([])
yticks((0.1:0.2:0.7))
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches',...
    'ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])
% print('HighSpike.svg', '-dsvg')