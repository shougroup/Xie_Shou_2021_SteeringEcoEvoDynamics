clear
close all

C = 3; % number of cycles
Fnum = (1:3); % number of simulation replicate

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

% plot evolutionary dynamics of fp
figure(1)
plot([1 C],[0.41 0.41],'k--','Linewidth',1.5) % optimal fp
hold off
axis([1 C 0 0.5])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('f_P','FontSize',16,'FontName','Arial','fontweight','bold');

% plot evolutionary dynamics of P(T)
figure(2)
plot([1 C],[2735.5 2735.5], 'k--','Linewidth', 1.5) % maximal community function P(T)
hold off
axis([1 C 0 3000])
set(gca,'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04 0.04])%,'xticklabel',[])
xlabel('Cycles','FontSize',16,'FontName','Arial','fontweight','bold');
ylabel('P(T)','FontSize',16,'FontName','Arial','fontweight','bold');

