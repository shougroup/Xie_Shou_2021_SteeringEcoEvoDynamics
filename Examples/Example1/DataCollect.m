% Process the raw simulation data for plotting evolutionary dynamics.
% Fnum: name of the directory (a number) that contains simulation output
% C: number of cycles.
function DataCollect(Fnum, C)
dataname = ['PlotData/Data' num2str(Fnum)];
pcs = 1e-9; % precision of zero
%%
fp0_commmean = zeros(C,1); % inter-community average of fp
M0frac_commmean = zeros(C,1); % inter-community average of \phi_M(0)
MTfrac_commmean = zeros(C,1); % inter-community average of \phi_M(T)
P_commmean = zeros(C,1); % inter-community average of community function P(T)

for n = 1:C
    load([num2str(Fnum) '/C' num2str(n) '/comm_all/newborns.mat']);
    load([num2str(Fnum) '/C' num2str(n) '/comm_selected.mat']);
    rseed_NB = [newborns.rseed];
    pnum_NB = [newborns.parentnum];
    comm_type_num = length(comm_selected);
    M0 = zeros(comm_type_num, 1); % total M biomass in selected Newborns
    H0 = zeros(comm_type_num, 1); % total H biomass in selected Newborns
    MT = zeros(comm_type_num, 1); % total M biomass in selected Adults
    HT = zeros(comm_type_num, 1); % total H biomass in selected Adults
    fp0_mean = zeros(comm_type_num, 1); % intra-community average fp of selected Newborns
    for i = 1:comm_type_num
        M0(i) = comm_selected(i).M_t(1);
        H0(i) = comm_selected(i).H_t(1);
        MT(i) = comm_selected(i).M_t(end);
        HT(i) = comm_selected(i).H_t(end);
        % find selected Newborns based on rseed and parentnum
        c_num = find((rseed_NB == comm_selected(i).rseed) & (pnum_NB == comm_selected(i).parentnum));
        if length(c_num) ~= 1
            error('error in finding selected Newborns in cycle %d', n)
        end
        if sum(newborns(c_num).M_L) < pcs
            fp0_mean(i) = 0;
        else
            fp0_mean(i)=sum(newborns(c_num).fp.*newborns(c_num).M_L)/sum(newborns(c_num).M_L);
        end
    end
    fp0_commmean(n) = mean(fp0_mean);
    M0frac_commmean(n) = mean(M0./(M0+H0));
    MTfrac_commmean(n) = mean(MT./(MT+HT));
    P_commmean(n) = mean([comm_selected.P]);
end

save(dataname,'fp0_commmean','P_commmean','M0frac_commmean','MTfrac_commmean')