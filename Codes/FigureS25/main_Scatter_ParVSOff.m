% clear
c = 780;
load(['C' num2str(c-1) '/comm_all/newborns'])
load(['C' num2str(c-1) '/comm_all/P_all'])
load(['C' num2str(c-1) '/comm_all/Pn'])
load(['C' num2str(c-1) '/ParResults'])
load(['C' num2str(c) '/OffResults'])

off_rep = 6;
sl = 4;
test_rep = length(heri_par_idx);
[P_temp, I] = sort(Pn + P_all, 'descend');
P_par1 = P_all(I(1:test_rep));
P_par_m1 = P_temp(1:test_rep);
nb_sorted = newborns(I);
M0_par1 = zeros(1, test_rep);
H0_par1 = zeros(1, test_rep);
for i = 1:test_rep
    M0_par1(i) = sum(nb_sorted(i).M_L);
    H0_par1(i) = sum(nb_sorted(i).H_L);
end
P_par = P_par + Pn_par;
for i = 1:sl
    [P_temp, I] = sort(P_par(i, :), 'descend'); 
    M0_par(i, :) = M0_par(i, I);
    H0_par(i, :) = H0_par(i, I);
    P_par(i, :) = P_par(i, I);
%     Pn_par(i, :) = Pn_par(i, I);
end
M0_par = [M0_par1; M0_par(:, 1:test_rep)];
H0_par = [H0_par1; H0_par(:, 1:test_rep)];
% P_par_m = [P_par_m1; P_par(:, 1:test_rep)];
P_par = [P_par_m1; P_par(:, 1:test_rep)];
phiM0_par = M0_par./(M0_par+H0_par);
phiM0_off = M0_off./(M0_off+H0_off);
P_off = P_off + Pn_off;
%%
szs = [0.5 0.8 0.5 0.8; 0.3 0.7 0.3 0.7; 0.1 0.4 0.1 0.4;];
titles = {'', 'B: 30%-H spiking', 'A: 60%-H spiking'};
close all
for k = [3 2]
    x = phiM0_par(k, :);
    y = mean(phiM0_off(1+(k-1)*off_rep:k*off_rep, :));
    pfit = polyfit(x, y, 1);
    figure
    scatter(phiM0_par(k, :), mean(phiM0_off(1+(k-1)*off_rep:k*off_rep, :)), 36, 'k', 'filled')
    hold on
    plot(szs(k, 1:2), szs(k, 1:2)*pfit(1)+pfit(2), 'r--', 'linewidth', 2)
    hold off
    axis(szs(k, :))
    xlabel('\phi_{M}(0) of parent')
    ylabel('ave. \phi_{M}(0) of offspring')
    title(titles{k})
    set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])
%     print([num2str(k) 'phiM0.svg'], '-dsvg')
end
%%
szs = [0 650 0 650];
% close all
for k = [3 2]
    x = P_par(k, :);
    y = nanmedian(P_off(1+(k-1)*off_rep:k*off_rep, :));
    pfit = polyfit(x, y, 1);
    figure
scatter(P_par(k, :), nanmedian(P_off(1+(k-1)*off_rep:k*off_rep, :)), 36, 'k', 'filled')
    axis(szs)
    hold on
    plot(szs(1:2), szs(1:2)*pfit(1)+pfit(2), 'r--', 'linewidth', 2)
    hold off
    xlabel('P(T) of Parent')
    ylabel('ave. P(T) of offspring')
    title(titles{k})
    set(gca,'LineWidth',3,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','ticklength',[0.04 0.04],'position',[1 0.7 3 3], 'color', 'none')%,'position',[1 0.7 3 3],'yticklabel',[])
%     print([num2str(k) 'PT.svg'], '-dsvg')
end

