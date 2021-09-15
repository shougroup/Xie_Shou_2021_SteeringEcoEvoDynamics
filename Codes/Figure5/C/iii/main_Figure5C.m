% plot Figure 5C(ii).
clear

load('C549/comm_all/adults')
P_parent = [adults.P];
P_off = zeros(100, 1);
P_errbar_upper = zeros(100, 1);
P_errbar_lower = zeros(100, 1);

for i =1:100
    load(['HeriData/adults' num2str(i)]) % simulation results for Adults offspring from the ith parent
    P_off_temp = [adults_off.P];
    P_off_temp = P_off_temp(P_off_temp>1e-5);
    P_off(i) = mean(P_off_temp); % average offspring community function
    P_errbar_upper(i) = quantile(P_off_temp, 0.75); % 75% quantile of offspring community function
    P_errbar_lower(i) = quantile(P_off_temp, 0.25); % 75% quantile of offspring community functionend
end
%%
pfit = polyfit(P_parent, P_off, 1);
figure(1)
errorbar(P_parent,P_off,P_errbar_lower-P_off,P_errbar_upper-P_off,'ko')
hold on
plot([800 1100], [800 1100]*pfit(1)+pfit(2), '--', 'color', [0 0.5 0], 'linewidth', 2)
hold off
axis([800 1100 800 1100])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
xlabel('Parent function')
ylabel('Offspring function')
title(['Heritability (slope) = ' num2str(pfit(1),3)])
    