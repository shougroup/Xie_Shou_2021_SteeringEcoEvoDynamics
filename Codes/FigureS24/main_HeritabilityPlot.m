clear
load('C549/comm_all/adults')
P_parent = [adults.P];
P_off = zeros(100, 1);
P_errbar_upper = zeros(100, 1);
P_errbar_lower = zeros(100, 1);
for i =1:100
    load(['Data/adults' num2str(i)])
    P_temp = [adults_off.P];
    P_temp = P_temp(P_temp > 1e-9);
    P_off(i) = mean(P_temp);
    P_errbar_upper(i) = quantile(P_temp, 0.75); % 75% quantile of offspring community function
    P_errbar_lower(i) = quantile(P_temp, 0.25); % 75% quantile of offspring community function
end
%%
figure(1)
% plot(P_parent,P_off,'k.','markersize',20)
errorbar(P_parent,P_off,P_errbar_lower-P_off,P_errbar_upper-P_off,'ko')
axis([300 1000 300 1000])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
% print('Heritability_errorbar.pdf', '-dpdf')  

%% plot heritability of the next cycle
load('C550/comm_all/adults')
P_parent = [adults.P];
P_off = zeros(100, 1);
P_errbar_upper = zeros(100, 1);
P_errbar_lower = zeros(100, 1);
for i =1:100
    load(['C551/adults' num2str(i)])
    P_temp = [adults_off.P];
    P_temp = P_temp(P_temp > 1e-9);
    P_off(i) = mean(P_temp);
    P_errbar_upper(i) = quantile(P_temp, 0.75); % 75% quantile of offspring community function
    P_errbar_lower(i) = quantile(P_temp, 0.25); % 75% quantile of offspring community function
end
%%
figure(2)
% plot(P_parent,P_off,'k.','markersize',20)
errorbar(P_parent,P_off,P_errbar_lower-P_off,P_errbar_upper-P_off,'ko')
axis([300 1000 300 1000])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])