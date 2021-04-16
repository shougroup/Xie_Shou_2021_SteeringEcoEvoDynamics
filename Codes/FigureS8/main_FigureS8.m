clear
Tm = (13:4:21);
% Tm = 17;
phiM0m = (0.05:0.05:0.95);
[Tmesh, phiM0_mesh] = meshgrid(Tm, phiM0m);
P_anc = zeros(size(Tmesh));
P_evo = zeros(size(Tmesh));
anc_fp = 0.13;
evo_fp = 0.38;
R0 = 1;
N0 = 100;
gM_max = 0.7;
gH_max = 0.3;
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR = 1/3;
K_HR = 1/5;
K_MB = 100/3;
dM = gM_max*5e-3;
dH = gH_max*5e-3;
options=odeset('RelTol',1e-6,'abstol',1e-10);
para = [gM_max; gH_max; c_BM; c_RM; c_RH; anc_fp; K_MR; K_HR; K_MB; dM; dH];
% solve differential equations
for i = 1:length(Tmesh(:))
    T0 = Tmesh(i);
    phiM0 = phiM0_mesh(i);
    [T, X] = ode15s(@(t,x) MHDynamics(t,x,para), [0 T0], [N0*phiM0; N0*(1-phiM0); 0; R0; 0],options);
    P_anc(i) = X(end, 5);
end
para = [gM_max; gH_max; c_BM; c_RM; c_RH; evo_fp; K_MR; K_HR; K_MB; dM; dH];
% solve differential equations
for i = 1:length(Tmesh(:))
    T0 = Tmesh(i);
    phiM0 = phiM0_mesh(i);
    [T, X] = ode15s(@(t,x) MHDynamics(t,x,para), [0 T0], [N0*phiM0; N0*(1-phiM0); 0; R0; 0],options);
    P_evo(i) = X(end, 5);
end
%%
for i = 1:length(Tm)
    figure(i)
    plot(phiM0m, P_anc(:, i), 'b', 'linewidth', 3)
    hold on
    plot(phiM0m, P_evo(:, i), 'r', 'linewidth', 3)
    hold off
    axis([0 1 0 3000])
    xlabel('\phi_{M}(0)')
    ylabel('P(T)')
    title(['T = ' num2str(Tm(i))])
    set(gca,'LineWidth',3,'FontSize',18,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'yticklabel',[])
end