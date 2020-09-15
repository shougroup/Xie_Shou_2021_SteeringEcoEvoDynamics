% plot Figure 5C.
clear
load('R1/adults')
P_parent = [adults.P];
P_off = zeros(100, 1);
P_std = zeros(100, 1);
for i =1:100
    load(['Data/adults' num2str(i)]) % simulation results for Adults offspring from the ith parent
    P_off(i) = mean([adults_off.P]); % average offspring community function
    P_std(i) = std([adults_off.P], 1); % standard diviation of offspring community function
end
%%
figure(1)
errorbar(P_parent,P_off,P_std,'ko')
axis([800 1100 800 1100])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])

    