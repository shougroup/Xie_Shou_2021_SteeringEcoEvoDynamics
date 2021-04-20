% Collect info about cycles when heritability is checked. It copies ParResults and
% OffResults to the folder. 
% newborns: newborns of the parents under current spiking strategy.
% P_all and Pn: Comm function and noise of the parents under current spiking strategy.
function HeriCheckCycles(Fnum, s)
fname1=['PlotData/Check' num2str(Fnum)];
if ~exist(fname1,'dir')
    mkdir(fname1)
end
load([num2str(Fnum) '/comm_all/check_cycle_m'])
check_counter = length(check_cycle_m);
spike_before_m = zeros(check_counter, s);
spike_after_m = zeros(check_counter, s);
lb_m = zeros(check_counter, s);
ub_m = zeros(check_counter, s);
corr_mean = zeros(check_counter, s);
for i = 1:check_counter
    n = check_cycle_m(i);
    fname2 = [fname1 '/C' num2str(n)];
    if ~exist(fname2,'dir')
        mkdir(fname2)
    end 
    load([num2str(Fnum) '/C' num2str(n-1) '/ParResults'])
    load([num2str(Fnum) '/C' num2str(n) '/OffResults'])
    copyfile([num2str(Fnum) '/C' num2str(n-1) '/ParResults.mat'], fname2)
    copyfile([num2str(Fnum) '/C' num2str(n) '/OffResults.mat'], fname2)
    copyfile([num2str(Fnum) '/C' num2str(n-1) '/comm_all'], fname2)
    spike_before_m(i, :) = spike_all;
    corr_mean(i, :) = heri';
    lb_m(i, :) = lb';
    ub_m(i, :) = ub';
    load([num2str(Fnum) '/C' num2str(n) '/spike_all'])
    spike_after_m(i, :) = spike_all;
end
%%
save([fname1 '/CheckSummary.mat'], 'check_cycle_m', 'spike_before_m', 'spike_after_m',...
    'lb_m', 'ub_m', 'corr_mean')
    
    