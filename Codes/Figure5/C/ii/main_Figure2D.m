% Plot the determinants, ave. fp(0) and phiM(0) of offspring vs. parent
clear
load('C549/comm_all/newborns')
fp0_parent = zeros(100, 1);
fp0_off = zeros(100, 1);
fp0_errbar_upper = zeros(100, 1);
fp0_errbar_lower = zeros(100, 1);
phiM0_parent = zeros(100, 1);
phiM0_off = zeros(100, 1);
phiM0_errbar_upper = zeros(100, 1);
phiM0_errbar_lower = zeros(100, 1);
for i =1:100
    fp0_parent(i) = sum(newborns(i).fp .* newborns(i).M_L)/sum(newborns(i).M_L);
    phiM0_parent(i)= sum(newborns(i).M_L) / (sum(newborns(i).M_L)+sum(newborns(i).H_L));
end

for i = 1:100
    load(['HeriData/newborns' num2str(i)])
    l = length(newborns);
    temp_fp = zeros(l, 1);
    temp_phi = zeros(l, 1);
    for j = 1:l
        temp_fp(j) = sum(newborns(j).fp .* newborns(j).M_L)/sum(newborns(j).M_L); 
        temp_phi(j) = sum(newborns(j).M_L) / (sum(newborns(j).M_L)+sum(newborns(j).H_L));
    end
    fp0_off(i) = mean(temp_fp);
    fp0_errbar_upper(i) = quantile(temp_fp, 0.75);
    fp0_errbar_lower(i) = quantile(temp_fp, 0.25);
    phiM0_off(i) = mean(temp_phi);
    phiM0_errbar_upper(i) = quantile(temp_phi, 0.75);
    phiM0_errbar_lower(i) = quantile(temp_phi, 0.25);
end

%% Plot w errorbar
figure(1)
errorbar(fp0_parent, fp0_off, fp0_errbar_lower-fp0_off, fp0_errbar_upper-fp0_off, 'ko')
% plot(fp0_parent,fp0_off,'k.','markersize',20)
axis([0.13 0.15 0.13 0.15])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
xlabel('$\overline{f_{P}}(0)$ of parents', 'interpreter', 'latex') 
ylabel('$\overline{f_{P}}(0)$ of offspring', 'interpreter', 'latex')

figure(2)
errorbar(phiM0_parent, phiM0_off, phiM0_errbar_lower-phiM0_off, phiM0_errbar_upper-phiM0_off,'ko')
% plot(phiM0_parent,phiM0_off,'k.','markersize',20)
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
axis([0.5 0.9 0.5 0.9])
xlabel('\phi_{M}(0) of parents') 
ylabel('\phi_{M}(0) of offspring')