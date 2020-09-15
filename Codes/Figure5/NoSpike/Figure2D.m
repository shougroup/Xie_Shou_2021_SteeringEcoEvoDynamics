% Plot the determinants, ave. fp(0) and phiM(0) of offspring vs. parent
clear
load('C549/comm_all/newborns')
fp0_parent = zeros(100, 1);
fp0_off = zeros(100, 1);
fp0_std_off = zeros(100, 1);
phiM0_parent = zeros(100, 1);
phiM0_off = zeros(100, 1);
phiM0_std_off = zeros(100, 1);
for i =1:100
    fp0_parent(i) = sum(newborns(i).fp .* newborns(i).M_L)/sum(newborns(i).M_L);
    phiM0_parent(i)= sum(newborns(i).M_L) / (sum(newborns(i).M_L)+sum(newborns(i).H_L));
end

for i = 1:100
    load(['Data/newborns' num2str(i)])
    l = length(newborns);
    temp_fp = zeros(l, 1);
    temp_phi = zeros(l, 1);
    for j = 1:l
        temp_fp(j) = sum(newborns(j).fp .* newborns(j).M_L)/sum(newborns(j).M_L); 
        temp_phi(j) = sum(newborns(j).M_L) / (sum(newborns(j).M_L)+sum(newborns(j).H_L));
    end
    fp0_off(i) = mean(temp_fp);
    fp0_std_off(i) = std(temp_fp, 1);
    phiM0_off(i) = mean(temp_phi);
    phiM0_std_off(i) = std(temp_phi, 1);
end

%% Plot w errorbar
figure(1)
errorbar(fp0_parent,fp0_off,fp0_std_off,'ko')
% plot(fp0_parent,fp0_off,'k.','markersize',20)
axis([0.13 0.15 0.13 0.15])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
xlabel('Determinant of parents') 
ylabel('Determinant of offspring')

figure(2)
errorbar(phiM0_parent,phiM0_off,phiM0_std_off,'ko')
% plot(phiM0_parent,phiM0_off,'k.','markersize',20)
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
axis([0.5 0.9 0.5 0.9])
xlabel('Determinant of parents') 
ylabel('Determinant of offspring')