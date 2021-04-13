clear
clear
comm_type_num = 2;
comm_rep_num = 50;
rp = 100;
fp0_all = zeros(rp, comm_type_num*comm_rep_num);
fp0_sel = zeros(rp, comm_type_num);
fpT_sel = zeros(rp, comm_type_num);
for j = 1:rp
    folder_name = ['RepeatData/R' num2str(j)];
    load([folder_name '/newborns'])
    load([folder_name '/adults'])
    P_all(j, :) = [adults.P];
    [~, idx] = sort(P_all(j, :), 'descend');
    for i=1:comm_type_num*comm_rep_num
        fp0_all(j,i) = sum(newborns(i).fp.*newborns(i).M_L)/sum(newborns(i).M_L);
    end
    fp0_sel(j, :) = fp0_all(j, idx(1:comm_type_num));
    for i = 1:comm_type_num
        a_idx = idx(i);
        fpT_sel(j, i) = sum(adults(a_idx).fp .* adults(a_idx).M_L)/sum(adults(a_idx).M_L);
    end
end
%%
dfp0 = mean(fp0_sel, 2)-mean(fp0_all, 2);
dfp_total = mean(fpT_sel, 2)-mean(fp0_all, 2);
dfp_sel = fpT_sel-fp0_sel;
%%
figure(1)
grp = [zeros(length(dfp_total), 1); ones(length(dfp0), 1); 2*ones(length(dfp_sel(:)), 1)];
h = boxplot([dfp_total; dfp0; dfp_sel(:)], grp, 'labels',{'total','inter','intra'},...
    'symbol', 'o', 'orientation', 'horizontal');
xlim([-6 6]*1e-3)
xticks((-4:2:4)*1e-3)
set(h, 'LineStyle', '-', 'linewidth', 2)
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','ticklength',[0.04,0.04],'units','inches','position',[1 1 3 3])
