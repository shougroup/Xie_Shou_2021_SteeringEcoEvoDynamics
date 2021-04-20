clear
hold on
% directory of simulations w/  spiking
fnames1 = {'A/1', 'A/2', 'A/3'}; 
% directory of simulations w/o spiking
fnames2 = {'B/1', 'B/2', 'B/3'}; 

for Fnum = 1:3
    folder_name = fnames1{Fnum};
    P_scan_one_para;
    phiMss = phiM_SteadyState(para);
    plot(para_m, P,'r','linewidth',2)
    plot([phiMss phiMss], [0 3e3],'r--','linewidth',1);
    xticks((0:0.2:1))
    xlim(x_range)
    xlabel('\phi_{M}(0)')
    ylabel('P(T)')
    set(gca, 'box', 'on', 'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04,0.04])%,'yticklabel',[])
end

for Fnum = 1:3
    folder_name = fnames2{Fnum};
    P_scan_one_para;
    phiMss = phiM_SteadyState(para);
    plot(para_m, P,'k','linewidth',2)
    plot([phiMss phiMss], [0 3e3],'k--','linewidth',1);
%     xlim(x_range)
%     xlabel('phiM0')
%     ylabel('P(T)')
%     set(gca, 'box', 'on', 'LineWidth',2,'FontSize',16,'FontName','Arial','fontweight','bold','units','inches','position',[1 1 3 3],'ticklength',[0.04,0.04])%,'yticklabel',[])
end
hold off
