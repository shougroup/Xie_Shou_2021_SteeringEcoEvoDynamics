clear
% Percentage of Newborn biomass replaced with H biomass. For the orange attractor, spike_frac = 0.
% For the teal attractor corresponding to the 30%-H spiking strategy, spike_frac = 0.3.
spike_frac = 0;
% Parameters shown in Table 1
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

% range and increment for f_P
x0=0.01;
dx=0.01;
xt=0.99;
% range and increment for phi_M(0)
y0=0.01;
dy=0.01;
yt=0.99;

xsize=round((xt-x0)/dx+1);
ysize=round((yt-y0)/dy+1);
%% Calculate the map if it hasn't been done. Otherwise just load the data in FullMap.mat. 
if ~exist('FullMap.mat')
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
    save('FullMap', 'P', 'MfracDiff')
else
    load('FullMap')
end
%%
% plot 3D community function landscape in Figure 2C
figure(1)
contour3((x0:dx:xt), (y0:dy:yt), P', (min(P(:)) : 200 : max(P(:))),'linewidth',3)
cm = colormap('gray');
c = 1-cm;
colormap(c)
ylim([0.01 0.99])
xlim([0.01 0.99])
set(gca,'LineWidth',3,'FontSize',18,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'yticklabel',[])
grid off
xlabel('f_P')
ylabel('phi_M(0)')
zlabel('P(T)')
%%
% plot 2D community function landscape with attractors
MfracT = MfracDiff + ones(xsize, 1) * (y0 : dy : yt);
MfracDiff_adjusted = MfracT * (1-spike_frac) - ones(xsize, 1) * (y0 : dy : yt);
figure(2)
contour((x0:dx:xt), (y0:dy:yt), P', (min(P(:)) : 240 : max(P(:))), 'linewidth', 3);
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04,0.04])%,'yticklabel',[])
cm = colormap('gray');
c = 1-cm;
colormap(c)
hold on
contour((x0:dx:xt),(y0:dy:yt), MfracDiff_adjusted', [0 0], 'color',[1 0.2 0], 'linewidth',3)
hold off
ylabel('phi_M(0)')
xlabel('f_P')