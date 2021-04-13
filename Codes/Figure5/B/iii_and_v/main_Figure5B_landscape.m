clear

spike_frac = 0.3;
comm_type_num = 2;
N = 100;
r = 12;

% range and increment for fp
x0 = 0.13;
dx = 1e-4;
xt = 0.16;
% range and increment for B10frac
y0 = 0.3;
dy = 0.005;
yt = 0.7;
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
    T0 = 17;
    
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

fp0_all = zeros(N, 1);
Mfrac0_all = zeros(N, 1);
folder_name = ['RepeatData/R' num2str(r)];
load([folder_name '/newborns'])
load([folder_name '/adults'])
P_all = [adults.P];
[~, idx] = sort(P_all, 'descend');
for i = 1:N
    fp0_all(i) = sum(newborns(i).fp.*newborns(i).M_L)/sum(newborns(i).M_L);
    Mfrac0_all(i) = sum(newborns(i).M_L)/sum(newborns(i).M_L+newborns(i).H_L);
end
fp0_sel = fp0_all(idx(1:comm_type_num));
Mfrac0_sel = Mfrac0_all(idx(1:comm_type_num));
%%
MfracDiff_spike = (MfracDiff + ones(xsize, 1) * (y0:dy:yt))*(1-spike_frac)- ones(xsize, 1) * (y0:dy:yt);
figure(1)
contour((x0:dx:xt), (y0:dy:yt), P', (650:50:1.2e3), 'linewidth',2)
cm = colormap('gray');
c = 1-cm;
colormap(c)
caxis([650 1.2e3])
hold on
contour((x0:dx:xt), (y0:dy:yt), MfracDiff_spike', [0 0], 'color', [0 0.5 0.5], 'linewidth', 2)
plot(fp0_all, Mfrac0_all, 'k.', 'markersize', 18)
plot(fp0_sel, Mfrac0_sel, 'go', 'linewidth',2)
hold off
% colorbar
axis([x0 xt y0 yt])
xlabel('$\overline{f_{P}}(0)$', 'interpreter', 'latex')
ylabel('\phi_{M}(0)')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])


