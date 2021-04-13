function simu_one_cycle(n, spike_frac)
m = n-1;
comm_type_num = 2;
mut_rate = 1e-2;
% Parameters shown in Table 1
% evolved parameters shown in Table 1
gM_max = 0.7; % max growth rate of M
gH_max = 0.3; % max growth rate of H
dM = gM_max*5e-3; % death rate of M
dH = gH_max*5e-3; % death rate of H
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR = 1/3;
K_HR = 1/5;
K_MB = 100/3;
fp_Bound = 1; % fp is between 0 and 1

% N0: BM_0, R(0)=1, T0 is the maturation time
N0 = 100;
R0 = 1;
T0 = 17;

N = 100; % number of communities within a cycle
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4; % maximal number of cells in the community
t_bin = 0.05; % time step in the simulation
pcs=1e-15; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps

% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;

% structure for communities
comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));

const_struct=struct('t_binnum',t_binnum, 'max_popul', max_popul, 'comm_rep_num', comm_rep_num,...
    'comm_type_num', comm_type_num, 'pcs', pcs, 'BM_target', BM_target, 't_bin', t_bin,...
    'fp_Bound',  fp_Bound, 'R0',  R0, 'dM',  dM, 'dH',  dH, ...
    'gM_max', gM_max, 'gH_max', gH_max, 'K_MR',  K_MR, 'K_HR', K_HR, 'K_MB', K_MB, ...
    'c_BM', c_BM, 'c_RM', c_RM, 'c_RH', c_RH, 'mut_rate', mut_rate);

filename = ['C' num2str(m) '/comm_selected.mat'];
load(filename);
if exist('HeriData/distrng.mat', 'file')
    % use saved random number generator to obtain the same results
    load('HeriData/distrng.mat')
    rng(distrng)
else
    rng('shuffle')
    distrng = rng;
    save('HeriData/distrng', 'distrng')
end
folder_name = ['C' num2str(n) '/comm_all'];
if ~isfolder(folder_name)
    mkdir(folder_name)
end
comm_all(1:comm_type_num*comm_rep_num,1) = comm_struct;
rep_counter = 0;
% reproduce chosen Adult communities.
for i = 1 : comm_type_num
    if rep_counter >= comm_type_num*comm_rep_num
        break
    end
    dil_factor = floor((comm_selected(i).M_t(t_binnum+1)+comm_selected(i).H_t(t_binnum+1))/BM_target/(1-spike_frac));
    if dil_factor == 0
        continue
    end
    rep_num_temp = min(dil_factor,comm_rep_num);
    comm_all_idx = min(comm_type_num*comm_rep_num,rep_counter+rep_num_temp);
    comm_all(rep_counter+1:comm_all_idx) = fixBM0_spike(comm_selected(i), comm_struct, const_struct, dil_factor, BM_target*spike_frac, rep_counter, i);
    rep_counter = rep_counter+rep_num_temp;
end
% assign random number seed to Newborns of the next cycle
rseed = randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
for ri = 1 : comm_type_num*comm_rep_num
    comm_all(ri).rseed = rseed(ri);
end
newborns = comm_all;
save([folder_name '/newborns.mat'],'newborns');
parfor rep = 1:N
    if sum(newborns(rep).M_L) + sum(newborns(rep).H_L) < pcs
        error('An Newborn is empty')
    end
    comm_all(rep) = simu_one_comm(newborns(rep), comm_struct, const_struct);
end
adults = comm_all;
save([folder_name '/adults.mat'], 'adults')

[~, idx] = sort([adults.P], 'descend');
comm_selected = adults(idx(1:comm_type_num));
save(['C' num2str(n) '/comm_selected.mat'], 'comm_selected')


