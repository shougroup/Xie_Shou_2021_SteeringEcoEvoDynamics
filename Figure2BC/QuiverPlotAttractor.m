clear
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
%%
if ~exist('FullMap.mat')
P = zeros(xsize,ysize); % community function at different f_P and phi_M(0)
MfracDiff = zeros(xsize,ysize); % offspring phi_M(0) - parent phi_M(0) at different f_P and phi_M(0)

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
        MfracDiff(m,n) = X(end, 1)/(X(end, 1)+X(end, 2))* - j;
    end
end
save('FullMap', 'P', 'MfracDiff')
else
    load('FullMap')
end

%%

st1 = 0.1; % phi_M(0) starts from 0.1
st2 = 0.9; % phi_M(0) starts from 0.9
fp_sample1 = (0.1 : 0.2 : 0.3)'; % fp values for the quivers below the attractor
fp_sample2 = (0.1 : 0.2 : 0.9)'; % fp values for the quivers above the attractor
xidx1 = round(fp_sample1/0.01); % index of fp_sample1 to corresponding to matrix MfracDiff
xidx2 = round(fp_sample2/0.01); % index of fp_sample2 to corresponding to matrix MfracDiff
yidx1 = round(st1/0.01)*ones(size(fp_sample1)); % index of st1 to corresponding to matrix MfracDiff
yidx2 = round(st2/0.01)*ones(size(fp_sample2)); % index of st2 to corresponding to matrix MfracDiff

%% generate the quiver plot of Figure 2B 
figure(1)
% Plot the attractor
contour((x0:dx:xt), (y0:dy:yt), MfracDiff',[0 0],'color',[1 0.2 0],'linewidth',3)
hold on
set(gca,'LineWidth',3, 'FontSize',18, 'FontName','Arial', 'fontweight','bold', 'units','inches', 'position',[1 1 3 3], 'ticklength',[0.04 0.04])%,'yticklabel',[])
% for all fp values in fp_sample1, phiM0 start from st1. The corresponding phiM(T)-phiM(0) can be found from the
% corresponding element of MfracDiff. The next arrow then starts from phiM(T) of the previous cycle. Draw 3 consecutive arrows.  
for i = 1:3
    lidx = sub2ind(size(MfracDiff), xidx1, yidx1); % find the corresponding index 
    quiver(fp_sample1, yidx1*0.01, zeros(size(fp_sample1)), MfracDiff(lidx), 0, 'o','markersize',2,'color','k','linewidth',2,'showarrowhead','on','markerfacecolor','k')
    % the new starting point of a new arrow is the end of the current arrow
    yidx1=yidx1+round(MfracDiff(lidx)/0.01);
end

% for all fp values in fp_sample2, phiM0 start from st2. The corresponding phiM(T)-phiM(0) can be found from the
% corresponding element of MfracDiff. The next arrow then starts from phiM(T) of the previous cycle. Draw 3 consecutive arrows.  
for i=1:7
    lidx = sub2ind(size(MfracDiff), xidx2, yidx2);  % find the corresponding index 
    quiver(fp_sample2, yidx2*0.01, zeros(size(fp_sample2)), MfracDiff(lidx), 0, 'o','markersize',2,'color','k','linewidth',2,'showarrowhead','on','markerfacecolor','k')
    yidx2 = yidx2+round(MfracDiff(lidx)/0.01); % the new starting point of a new arrow is the end of the current arrow
    yidx2(yidx2==0) = 1; % bound the index corresponding to phiM(0) 
end

hold off

ylabel('phiM(0)')
xlabel('fp')


