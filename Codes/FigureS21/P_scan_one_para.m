% Calculate community function varying one phenotype while fixing the others
% clear
para_m = (0.05:0.01:0.95);
x_range = [0 1];
P = zeros(size(para_m));
% Fnum = 12;
cycle = 3000;
%% other phenotypes are taken from #cnum of the chosen communities from cycle in folder /Fnum
cnum = 1;
[para, BM0] = find_determinants_MutualSP(folder_name, cycle, cnum);
% %% other phenotypes are taken from the average of the chosen communities from cycle in folder /Fnum
% [para, BM0] = mean_determinants_ExploitSP(Fnum, cycle);
%%
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
dM = 0.7*5e-3;
dH = 0.3*5e-3;
% consts = [c_BM; c_RM; c_RH; dM; dH];
N0 = 100;
R0 = 1;
T0 = 17;
pcs = 1e-9;

options=odeset('RelTol',1e-6,'abstol',1e-10);

for j=1:length(para_m)
%     % varying phenotype is K_MB
%     para(9) = para_m(j);

%     varying determinant is phiM0
    BM0 = 100*[para_m(j), 1-para_m(j)];

%     para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; K_HA; dM; dH];
    % solve differential equations
    [T, X] = ode15s(@(t,x) MHMutualism(t,x,para), [0 T0], [BM0';0;R0;0],options);
    if ~isreal(X)
        error('imaginary value')
    elseif nnz(X < -1e-4)>0
        error('negative value')
    end
    P(j)= X(end, 5);
end
% %%
% hold on
% plot(para_m, P, 'linewidth',2)
% xlim(x_range)
% xlabel('phiM0')
% ylabel('P(T)')
% set(gca, 'box', 'on', 'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04,0.04])%,'yticklabel',[])
% % hold off
