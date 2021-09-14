clear
close all

C = 300; % number of cycles
Fnum = 1; % number of simulation replicate

cl = {'k','c',[0.8 0.8 0.8]}; % color of the line
if ~exist('PlotData','dir')
    mkdir('PlotData')
end
for i = 1:length(Fnum)
    filename=['PlotData/Data' num2str(Fnum(i)) '.mat'];
    if ~exist(filename, 'file')
        DataCollect(i, C);
    end
    load(filename)
    figure(1)
    plot((1:C),fp0_commmean(1:C), 'color', cl{i}, 'Linewidth', 1);
    hold on
    figure(2)
    plot((1:C), P_commmean(1:C), 'color', cl{i}, 'Linewidth',1);
    hold on
end

%%
figure(1)
plot((1:C),fp0_commmean(1:C),'k','Linewidth',1);
hold on
plot([-40 C],[0.41 0.41],'k--','Linewidth',1.5)
hold off
axis([-40 C 0 0.5])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 1 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('f_P(0)','FontSize',16,'FontName','Arial','fontweight','bold');


figure(2)
M0frac_commmean(1) = 0.6;
plot((1:C),M0frac_commmean(1:C),'k','Linewidth',1);
hold on
plot([-40 C],[0.54 0.54],'k--','Linewidth',1.5)
hold off
axis([-40 C 0.5 0.8])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 1 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('ϕ_M(0)','FontSize',16,'FontName','Arial','fontweight','bold');


figure(3)
plot((1:C),P_commmean(1:C),'k','Linewidth',1);
hold on
plot([-40 C],[2735.5 2735.5],'k--','Linewidth',1.5)
hold off
axis([-40 C 0 3000])
set(gca,'LineWidth',1,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 1 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('P(T))','FontSize',16,'FontName','Arial','fontweight','bold');
