clear

spike_frac = 0;
cycle = 24;
comm_type_num = 2;
N = 100;

% range and increment for f_P
x0=0.12;
dx=0.0005;
xt=0.15;
% range and increment for phi_M(0)
y0=0.5;
dy=0.01;
yt=0.9;
xsize = round((xt-x0)/dx+1);
ysize = round((yt-y0)/dy+1);

filename = 'LocalMap.mat';
if exist(filename, 'file')
    load(filename)
else
    gMmax = 0.7;
    gHmax = 0.3;
    c_BM = 1/3;
    c_RM = 1e-4;
    c_RH = 1e-4;
    K_MR = 1/3;
    K_HR = 1/5;
    K_MB = 100/3;
    dM = gMmax*5e-3;
    dH = gHmax*5e-3;
    
    % N0: BM_0, R(0)=1, T0 is the maturation time
    N0 = 100;
    R0 = 1;
    T0 = 20;
    
    options=odeset('RelTol',1e-6,'abstol',1e-10);
    P = zeros(xsize,ysize); % community function at different f_P and phi_M(0)
    MfracDiff = zeros(xsize,ysize); % offspring phi_M(0) - parent phi_M(0) at different f_P and phi_M(0), no spiking
    m = 0; % counter for values of f_P
    % calculate P(T) at different f_P and phi_M(0)
    for i = x0 : dx : xt
        m = m+1; % counter for values of phi_M(0)
        n = 0;
        fp = i;
        for j = y0:dy:yt
            % parameters for differential equations
            para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; dM; dH];
            % solve differential equations
            [T, X] = ode15s(@(t,x) MHDynamics(t,x,para), [0 T0], [N0*j; N0*(1-j); 0; R0; 0],options);
            if ~isreal(X)
                error('imaginary value')
            elseif nnz(X < -1e-5)>0
                error('negative value')
            end
            n = n+1;
            P(m, n)= X(end, 5);
            MfracDiff(m, n) = X(end, 1)/(X(end, 1)+X(end, 2)) - j;
        end
    end
    save(filename, 'P', 'MfracDiff')
end

%% plot the Newborns and selected Newborns on the landscape
load(['C' num2str(cycle) '/comm_all/newborns'])
load(['C' num2str(cycle) '/comm_selected'])
NB_rseed = [newborns.rseed];
NB_pnum = [newborns.parentnum];
fp0 = zeros(100, 1);
phiM0 = zeros(100, 1);
for i = 1:100
    fp0(i) = sum(newborns(i).fp .* newborns(i).M_L)/sum(newborns(i).M_L);
    phiM0(i) = sum(newborns(i).M_L)/(sum(newborns(i).M_L)+sum(newborns(i).H_L));
end
sel_idx = zeros(length(comm_selected), 1);
for i = 1:length(comm_selected)
    temp_idx = find(NB_rseed == comm_selected(i).rseed & NB_pnum == comm_selected(i).parentnum);
    if length(temp_idx) == 1
        sel_idx(i) = temp_idx;
    else
        error('cannot find NB of Adult %d in Cycle %d ', cnum, cycle)
    end
end
%%
MfracT = MfracDiff + ones(xsize, 1) * (y0 : dy : yt);
MfracDiff_adjusted = MfracT * (1-spike_frac) - ones(xsize, 1) * (y0 : dy : yt);
figure(2)
contour((x0:dx:xt), (y0:dy:yt), P', (min(P(:)) : 50 : max(P(:))), 'linewidth', 3);
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04,0.04])%,'yticklabel',[])
cm = colormap('gray');
c = 1-cm;
colormap(c)
hold on
contour((x0:dx:xt),(y0:dy:yt), MfracDiff_adjusted', [0 0], 'color',[1 0.2 0], 'linewidth',3)
scatter(fp0, phiM0,18,'k','filled')
scatter(fp0(sel_idx), phiM0(sel_idx), 'g', 'linewidth', 2)
hold off
axis([x0 xt y0 yt])
xlabel('$\overline{f_{P}}(0)$', 'interpreter', 'latex')
ylabel('\phi_{M}(0)')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])


