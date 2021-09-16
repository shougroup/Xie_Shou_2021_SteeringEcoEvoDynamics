clear
close all

spike_num = 5; % number of candidate spiking strategies
C = 2000; % number of cycles
print_flag = false; % whether or not print the figures.
rc = true; % whether or not heritability check was performed.
tile_flag = false; % tile_flag = true, figures are presented in tiles; otherwise, in individual figures
f_num = [1 2 3];
r_num = 3;
cl = {'k', 'c', [0.8 0.8 0.8]};

fp = {};
P = {};
K_MR = {};
K_MB = {};
K_HR = {};
b_Mmax = {};
b_Hmax = {};
if ~exist('PlotData','dir')
    mkdir('PlotData')
end
for i = 1:r_num
    filename=['PlotData/Data' num2str(f_num(i)) '.mat'];
    if ~exist(filename, 'file')
        DataCollect_SP(i, C);
    end
    load(filename)
    fp{i} = fp_commmean;
    P{i} = P_commmean;
    K_MR{i} = K_MR_commmean;
    K_MB{i} = K_MB_commmean;
    K_HR{i} = K_HR_commmean;
    b_Mmax{i} = g_Mmax_commmean;
    b_Hmax{i} = g_Hmax_commmean;
end
fc = 0; % counter for figures
%%
if tile_flag == true
    tiledlayout(2,4)
end
% plot P(T) from 3 runs
if tile_flag == true
    nexttile
else
    fc = fc+1;
    figure(fc)
end
hold on
for i = 1:r_num
    plot((1:length(P{i})),P{i},'color',cl{i},'Linewidth',1);
end
hold off
axis([1 C 0 3000])
xticks((1e3:1e3:C))
set(gca,'box','on','LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('P(T)','FontSize',16,'FontName','Arial','fontweight','bold');

%%
% plot fp from 3 runs
if tile_flag == true
    nexttile
else
    fc = fc+1;
    figure(fc)
end
hold on
for i = 1:r_num
    plot((1:length(fp{i})),fp{i},'color',cl{i},'Linewidth',1);
end
hold off
xticks((1e3:1e3:C))
axis([1 C 0 0.5])
set(gca,'box','on','LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('f_P','FontSize',16,'FontName','Arial','fontweight','bold');
%%
if tile_flag == true
    nexttile
else
    fc = fc+1;
    figure(fc)
end
hold on
for i = 1:r_num
    plot((1:length(b_Mmax{i})),b_Mmax{i},'color',cl{i},'Linewidth',1);
    plot((1:length(b_Hmax{i})),b_Hmax{i},'color',cl{i},'Linewidth',1);
end
plot([1 C],[0.7 0.7],'g--','Linewidth',1.5)
plot([1 C],[0.3 0.3],'g--','Linewidth',1.5)
hold off
xticks((1e3:1e3:C))
axis([1 C 0 0.85])
set(gca,'box','on','LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('b_{M,Hmax}','FontSize',16,'FontName','Arial','fontweight','bold');

%%
% plot 1/K_MR from 3 runs
if tile_flag == true
    nexttile
else
    fc = fc+1;
    figure(fc)
end
hold on
for i = 1:r_num
    plot((1:length(K_MR{i})),1./K_MR{i},'color',cl{i},'Linewidth',1);
end
plot([1 C],[3 3],'g--','Linewidth',1.5)
hold off
xticks((1e3:1e3:C))
axis([1 C 0 4])
% yticks([0 0.1 0.2 0.3 0.4])
% yticklabels({'0','1','2','3','4' })
set(gca,'box','on','LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('1/K_{MR}','FontSize',16,'FontName','Arial','fontweight','bold');

% plot 1/K_MB from 3 runs
if tile_flag == true
    nexttile
else
    fc = fc+1;
    figure(fc)
end
hold on
for i = 1:r_num
    plot((1:length(K_MB{i})),1./K_MB{i},'color',cl{i},'Linewidth',1);
end
plot([1 C],[0.030 0.030],'g--','Linewidth',1.5)
hold off
xticks((1e3:1e3:C))
yticks([0 0.006 0.012 0.018 0.024 0.03])
yticklabels({'0','6','12','18','24','30'})
axis([1 C 0 0.033])
set(gca,'box','on','LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('1/K_{MB}','FontSize',16,'FontName','Arial','fontweight','bold');

% plot K_HR from 3 runs
if tile_flag == true
    nexttile
else
    fc = fc+1;
    figure(fc)
end
hold on
for i = 1:r_num
    plot((1:length(K_HR{i})),1./K_HR{i},'color',cl{i},'Linewidth',1);
end
plot([1 C],[5 5],'g--','Linewidth',1.5)
hold off
axis([1 C 0 6])
xticks((1e3:1e3:C))
% yticks([0 0.1 0.2 0.3 0.4 0.5 0.6])
% yticklabels({'0','1','2','3','4','5','6' })
set(gca,'box','on','LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('1/K_{HR}','FontSize',16,'FontName','Arial','fontweight','bold');

%%
if rc == true
    if tile_flag == true
        nexttile
    else
        fc = fc+1;
        figure(fc)
    end
    for j = 1:r_num
        filename = ['PlotData/Check' num2str(f_num(j)) '/CheckSummary.mat'];
        if ~exist(filename, 'file')
            HeriCheckCycles(i, spike_num);
        end
        load(filename)
        frac = zeros(C, 1);
        check_cycle_m = [0; check_cycle_m(:)];
        for i = 1:length(check_cycle_m)-1
            frac(check_cycle_m(i)+1:check_cycle_m(i+1)) = spike_before_m(i, 1);
        end
        frac(check_cycle_m(end):C) = spike_after_m(end, 1);
        plot((1:C), frac, 'color', cl{j})
        hold on
    end
    hold off
    xticks((1e3:1e3:C))
    yticks([-0.6 -0.3 0 0.3 0.6])
    yticklabels({'60%-M','30%-M','0','30%-H','60%-H' })
    axis([1 C -0.7 0.7])
    set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04],'yticklabel',[])
    xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
    ylabel('Spikeing fraction','FontSize',16,'FontName','Arial','fontweight','bold');
end