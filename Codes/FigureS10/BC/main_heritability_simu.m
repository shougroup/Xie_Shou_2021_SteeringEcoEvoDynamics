% Generate data for parent-offspring comparison. For every parent community, reproduce
% into as many offspring communities as possible.

clear

spike_frac = 0; % fraction of H pure culture spiked in
n = 24;
comm_type_num = 1;

% Parameters shown in Table 1
% evolved parameters shown in Table 1
gM_max = 0.7; % max growth rate of M
gH_max = 0.3; % max growth rate of H
dM = gM_max*5e-3; % death rate of M
dH = gH_max*5e-3; % death rate of H
fp_start = 0.13; % fp at the beginning of selection
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
T0 = 20;
mut_rate = 1e-2;
N = 100;
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

% load the parent Adults
filename = ['C' num2str(n) '/comm_all/adults.mat'];
if ~exist(filename, 'file')
    load(['C' num2str(n) '/comm_all/newborns.mat'])
    adults(1:N, 1) = comm_struct;
    % rep is the index of communities within one cycle
    parfor rep = 1:N
        if sum(newborns(rep).M_L) + sum(newborns(rep).H_L) < pcs
            error('An Newborn is empty')
        end
        adults(rep) = simu_one_comm(newborns(rep), comm_struct, const_struct);
    end
    save(filename, 'adults')
else
    load(filename);
end
% load('R1/adults.mat')
if ~isfolder('HeriData')
    mkdir('HeriData')
end
filename = 'HeriData/distrng_m.mat';
if exist(filename, 'file')
    % set the randome number generator to obtain identical results
    load(filename)
else
    rng('shuffle')
    distrng_m = randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
    save(filename, 'distrng_m')
end

for j = 1:100
    rep_counter = 0;
    % reproduce each Adult communities.
    comm_selected = adults(j);
    rng(distrng_m(j),'twister');
    dil_factor = floor((comm_selected.M_t(t_binnum+1)+comm_selected.H_t(t_binnum+1))/BM_target/(1-spike_frac));
    if dil_factor == 0
        continue
    end
    comm_all = fixBM0_spike(comm_selected, comm_struct, const_struct, dil_factor, BM_target*spike_frac, rep_counter, j);
    % assign random number seed to Newborns of the next cycle
    rseed = randi(2^32-1, dil_factor,1,'uint32');
    for ri = 1 : length(comm_all)
        comm_all(ri).rseed = rseed(ri);
    end
    newborns = comm_all;
    % save the selected communities and the current state of the random number generator
    save(['HeriData/newborns' num2str(j)],'newborns');
    % rep is the index of communities within one cycle
    parfor rep = 1:length(comm_all)
        if sum(newborns(rep).M_L) + sum(newborns(rep).H_L) < pcs
            error('An Newborn is empty')
        end
        comm_all(rep) = simu_one_comm(newborns(rep), comm_struct, const_struct);
    end
    adults_off = comm_all;
    save(['HeriData/adults' num2str(j)],'adults_off');
end

