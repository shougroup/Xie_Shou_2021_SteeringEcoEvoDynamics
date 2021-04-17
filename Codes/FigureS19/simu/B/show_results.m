experiments = {'pipette_m1w100pu62bh', 'pipette_m1w100pu62no_bh', ...
    'pipette_m1w100pu62bh_nonull', 'pipette_m1w100pu62nobh_nonull'};

num_exp = length(experiments);
exp_struct = {};
f = figure;
for i = 1 : num_exp
    exp_struct{i} = load(strcat('pooled_results/',experiments{i}));
end
%%
set(0,'DefaultLineLineWidth',1.5,'DefaultTextFontWeight','bold',...
    'DefaultTextFontSize',14,'DefaultAxesFontWeight','bold',...
    'DefaultAxesFontSize',14, 'DefaultLineMarkerSize',10);

hold on
for i = 1 : num_exp
    thr_max = 2735.5 * exp_struct{i}.rc.multiplier;
    plot(exp_struct{i}.pdt_chosen_avg / thr_max,'LineWidth',2)
end
xlabel('cycle number')
ylabel('\langleP(T)_{chosen}\rangle / P(T)_{max}')
legend('pick top 10', 'pick top ~2 (standard)','pick top 10, no null','pick top~2 no null','Location','northwest')
title('product over time')