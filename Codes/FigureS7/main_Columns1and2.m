clear
% uncomment one of the below 3 blocks to generate the corresponding figures

% "No spike"
phiM0 = [0.5 0.9];

% % "30%-H spike"
% spike_frac = 0.3;
% phiM0 = [0.5 0.9]*(1-spike_frac);

% % "30%-M spike"
% spike_frac = 0.3;
% phiM0 = [0.5 0.9]*(1-spike_frac)+spike_frac;

fp = 0.13;
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
para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; d1; d2];
c = {'b', [0 0.6 0]};
for j = 1:2
[T,X]=ode15s(@(t,x) MHDynamics(t,x,para),[0 T0],[N0*phiM0(j);N0*(1-phiM0(j));0;R0;0],options);
if ~isreal(X)
    error('imaginary value')
elseif nnz(X<-1e-5)>0
    error('negative value')
end

figure(1)
plot(T, X(:, 4), 'color', c{j}, 'linewidth', 2)
hold on

figure(2)
plot(T, X(:, 1)./(X(:, 1)+X(:, 2)), 'color', c{j}, 'linewidth', 2)
hold on
end

figure(1)
hold off
axis([0 T0 0 1])
xticks([])
xlabel('t')
ylabel('R(t)')
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches',...
    'ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])

figure(2)
plot([0 T0], [phi_SS phi_SS], ':', 'linewidth', 3, 'color', [1 0.5 0])
hold off
axis([0 T0 0 1])
xticks([])
xlabel('t')
ylabel('\phi_{M}(t)')
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches',...
    'ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])

 