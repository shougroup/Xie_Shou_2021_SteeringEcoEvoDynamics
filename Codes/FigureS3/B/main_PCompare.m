clear
ax = [400 1300 400 1300];
load('C2000/comm_all/newborns')
load('C2000/comm_all/P_all')
P_cal = zeros(1, 100);
d1 = 0.7*5e-3;
d2 = 0.3*5e-3;
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
R0 = 1;
T0 = 17;
options=odeset('RelTol',1e-6,'abstol',1e-10);
for i = 1:100
    M0 = sum(newborns(i).M_L);
    H0 = sum(newborns(i).H_L);
    fp = sum(newborns(i).M_L .* newborns(i).fp)/M0;
    gMmax = sum(newborns(i).M_L .* newborns(i).gM_max)/M0;
    gHmax = sum(newborns(i).H_L .* newborns(i).gH_max)/H0;
    K_MR = sum(newborns(i).M_L .* newborns(i).K_MR)/M0;
    K_MB = sum(newborns(i).M_L .* newborns(i).K_MB)/M0;
    K_HR = sum(newborns(i).H_L .* newborns(i).K_HR)/H0;
    para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; d1; d2];
    [T,X]=ode15s(@(t,x) MHDynamics(t,x,para),[0 T0],[M0;H0;0;R0;0],options);
    if ~isreal(X)
        error('imaginary value')
    elseif nnz(X<-1e-5)>0
        error('negative value')
    end
    P_cal(i) = X(end,5);
end
figure(1)
scatter(P_cal, P_all, 'k', 'linewidth', 2)
hold on
plot(ax(1:2), ax(1:2), 'b--', 'linewidth', 3)
hold off
axis(ax)
set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','ticklength',[0.04 0.04],'position',[1 0.7 3 3])%,'position',[1 0.7 3 3],'yticklabel',[])
xlabel('P(T) calculated')
ylabel('P(T) from simulation')
% print('PCompare.svg', '-dsvg')