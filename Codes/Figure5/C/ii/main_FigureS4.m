% Plot the determinants, ave. fp(0) and phiM(0) of offspring vs. parent
clear
load('C549/comm_all/newborns')
fp0_parent = zeros(100, 1);
fp0_off = zeros(100, 200);
fp0_off_norm = zeros(100, 200);
fp0_std_off = zeros(100, 1);
phiM0_parent = zeros(100, 1);
phiM0_off = zeros(100, 200);
phiM0_off_norm = zeros(100, 200);
phiM0_std_off = zeros(100, 1);
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
    fp0_off(i, 1:l) = temp_fp;
    fp0_off_norm(i, 1:l) = (temp_fp - mean(temp_fp))/std(temp_fp, 1);
    fp0_std_off(i) = std(temp_fp, 1);
    phiM0_off(i, 1:l) = temp_phi;
    phiM0_off_norm(i, 1:l) = (temp_phi - mean(temp_phi))/std(temp_phi, 1);
    phiM0_std_off(i) = std(temp_phi, 1);
end

%%
figure(1)
histogram(fp0_off_norm(abs(fp0_off_norm(:)) > 1e-9), (-5:0.25:5), 'normalization', 'pdf', 'linewidth', 1, 'facecolor', 'none')
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.02 0.02])%,'xticklabel',[])
ylabel('PDF','FontSize',16,'FontName','Arial','fontweight','bold');
xlim([-4.9 4.9])
xlabel('Normalized $\overline{f_{P}}(0)$', 'interpreter', 'latex')
ylabel('PDF')
figure(2)
histogram(phiM0_off_norm(abs(phiM0_off_norm(:)) > 1e-9),(-5:0.25:5), 'normalization', 'pdf', 'linewidth', 1, 'facecolor', 'none')
hold on
x = (-5:0.1:5);
plot(x, normpdf(x), 'b', 'linewidth', 2)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.02 0.02])%,'xticklabel',[])
ylabel('PDF','FontSize',16,'FontName','Arial','fontweight','bold');
hold off
xlim([-4.9 4.9])
xlabel('normalized \phi_{M}(0)')
ylabel('PDF')

