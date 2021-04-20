function SPMutualParaInit
% Set the improvement threshold that would trigger a heritability check
% If the improvement over the last recent_cycles is less than inc_thre, check heritability 
s_check.min_check_int = 6000;
s_check.inc_thre = 0;
s.spike_initial = [0 0.3 0.6 -0.3 -0.6];
s.spike_clone_num = 5;
spike_clone_num = s.spike_clone_num;
% number of offspring communities from each parent community to test
% heritability
s.off_rep_max = 6;
s.n_bstrap = 1e3; % number of bootstraping for heritability
% 1-q is the confidence interval of the heritability
s.q = 0.05;
% when checking the heritability, only use up to the top 50 parent Adults
s.test_rep_max = 100;
s.HeriSwitch = 1;
% upper bound for gH_max. for gH_max_Bound = 0.3, set V = 1.
s.gH_max_Bound = 0.3;
% the factor for amount of R(0)
s.V = 1;
% minimal number of Adults allowed to reproduce. comm_type_num = 1 for the
% top-dog strategy, comm_type_num = n for the top n% strategy.
s.comm_type_num = 10;
% mutation rate corresponding to effective mutation rate of 2e-3
% to turn off the mutation, set mut_rate=0.
s.mut_rate = 1e-2; % mutation rate corresponding to effective mutation rate of 2e-3
s.Pn_sig = 50; % std of the measurement noise.
s.N = 100; % number of communities within a cycle
% BM_target is the target biomass, T0 is the maturation time
s.BM_target = 100;
s.T0 = 20;
% comm_type_num * comm_rep_num = number of communities within one cycle.
s.comm_rep_num = s.N/s.comm_type_num; % maximal number of offspring community from one Adult.
s.max_popul = 1e4*s.V; % maximal number of cells in an Adult
max_popul = s.max_popul;
s.nb_popul = s.BM_target*2; % maximal number of cells in a Newborn
nb_popul = s.nb_popul;
s.t_bin = 0.05; % time step in the simulation
s.pcs = 1e-9; % precision constant
s.t_binnum = int16(s.T0/s.t_bin); % number of time steps
t_binnum = s.t_binnum;
% ancestral or evolved parameters shown in Table 1
s.gM_max_start = 0.7/1.2; % max growth rate of M
s.gH_max_start = 0.3/1.2; % max growth rate of H
s.dM = 3.5e-3; % death rate of M
s.dH = 1.5e-3; % death rate of H
s.c_BM = 1/3;
s.c_RM = 1e-4;
s.c_RH = 1e-4;

s.fp_start = 0.1; % fp at the beginning of selection
s.K_MR_start = 1 ;
s.K_HR_start = 1 ;
s.K_MB_start = 500/3 ;
s.B0_start = 200/3 ;

% evolutionary bounds of the phenotypes
s.gM_max_Bound = 0.7;
s.fp_Bound = 1; % fp is between 0 and 1
s.K_MB_Bound = 100/3 ;
s.K_MR_Bound = 1/3 ;
s.K_HR_Bound = 0.2 ;
s.B0_Bound = 1000/3;
% if K_MR, K_HR, K_MB is larger than K_singular value, they are effectively
% infinity
s.K_singular = 1e5;

% R(0)=1 for each V
s.R0 = s.V;
s.sl = length(s.spike_initial)-1; % number of spiking fraction to be tested

% structure of an Adult community
s.comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'gM_max', zeros(max_popul,1), 'K_MB', zeros(max_popul,1), 'K_MR', zeros(max_popul,1),...
    'gH_max', zeros(max_popul,1), 'K_HR', zeros(max_popul,1), 'B0', zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
% structure of a Newborn community
s.newborn_struct = struct('M_L',zeros(nb_popul,1),'H_L',zeros(nb_popul,1),'fp',zeros(nb_popul,1),...
    'gM_max', zeros(nb_popul,1), 'K_MB', zeros(nb_popul,1), 'K_MR', zeros(nb_popul,1),...
    'gH_max', zeros(nb_popul,1), 'K_HR', zeros(nb_popul,1), 'B0', zeros(nb_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
% structure of simulation constants
s.const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',s.comm_rep_num,...
    'comm_type_num',s.comm_type_num,'pcs',s.pcs,'BM_target',s.BM_target, 't_bin', s.t_bin,...
    'gM_max_Bound', s.gM_max_Bound, 'gH_max_Bound', s.gH_max_Bound, 'fp_Bound', s.fp_Bound, ...
    'K_MB_Bound', s.K_MB_Bound,'K_MR_Bound', s.K_MR_Bound, 'B0_Bound', s.B0_Bound, ...
    'K_HR_Bound', s.K_HR_Bound, 'R0', s.R0, 'K_singular', s.K_singular, 'dM', s.dM, 'dH', s.dH, ...
    'c_BM', s.c_BM, 'c_RM', s.c_RM, 'c_RH', s.c_RH, 'mut_rate', s.mut_rate);
s.H_isolates_struct = struct('gH_max', -ones(spike_clone_num, 1), 'K_HR', -ones(spike_clone_num, 1),...
    'B0', -ones(spike_clone_num, 1), 'H_L', -ones(spike_clone_num, 1));
s.M_isolates_struct = struct('gM_max', -ones(spike_clone_num, 1), 'K_MR', -ones(spike_clone_num, 1),...
                'M_L', -ones(spike_clone_num, 1), 'fp', -ones(spike_clone_num, 1),...
                'K_MB', -ones(spike_clone_num, 1));
save('SPMutualPara.mat', '-struct', 's')
save('CheckDecPara.mat', '-struct', 's_check')