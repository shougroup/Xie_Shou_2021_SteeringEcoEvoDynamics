clear
% dbstop if error
% HeriSwitch=0: random, HeriSwitch=1: highest heri, HeriSwitch=-1: lowest heri
HeriSwitch = int8(1);
% number of previous cycles
C_prev = 777;
spike_initial = [0 0.3 0.6 -0.3 -0.6];
spike_clone_num = 10;
C = 780; % total number of cycles
% number of offspring communities from each parent community to test
% heritability
off_rep_max = 6;
n_bstrap = 1e3; % number of bootstraping for heritability
% 1-q is the confidence interval of the heritability
q = 0.05;
% when checking the heritability, only use up to the top 50 parent Adults
test_rep_max = 100;
% check heritability every check_period cycles.
check_period = 100;
% upper bound for gH_max. for gH_max_Bound = 0.3, set V = 1.
gH_max_Bound = 0.3;
% the factor for amount of R(0)
V = 1;
% minimal number of Adults allowed to reproduce. comm_type_num = 1 for the
% top-dog strategy, comm_type_num = n for the top n% strategy.
comm_type_num = 10;
% mutation rate corresponding to effective mutation rate of 2e-3
% to turn off the mutation, set mut_rate=0.
mut_rate = 1e-2; % mutation rate corresponding to effective mutation rate of 2e-3
Pn_sig = 50; % std of the measurement noise.
N = 100; % number of communities within a cycle
% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;
T0 = 17;
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4*V; % maximal number of cells in an Adult
nb_popul = BM_target*2; % maximal number of cells in a Newborn
t_bin = 0.05; % time step in the simulation
pcs = 1e-9; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps

% ancestral or evolved parameters shown in Table 1
gM_max_start = 0.7/1.2; % max growth rate of M
gH_max_start = 0.3/1.2; % max growth rate of H
dM = 3.5e-3; % death rate of M
dH = 1.5e-3; % death rate of H
fp_start = 0.1; % fp at the beginning of selection
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
K_MR_start = 1 ;
K_HR_start = 1 ;
K_MB_start = 500/3 ;

% evolutionary bounds of the phenotypes
gM_max_Bound = 0.7;
fp_Bound = 1; % fp is between 0 and 1
K_MB_Bound = 100/3 ;
K_MR_Bound = 1/3 ;
K_HR_Bound = 0.2 ;


% R(0)=1 for each V
R0 = V;
% if K_MR, K_HR, K_MB is larger than K_singular value, they are effectively
% infinity
K_singular = 1e3;
% structure of an Adult community
comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'gM_max', zeros(max_popul,1), 'K_MB', zeros(max_popul,1), 'K_MR', zeros(max_popul,1),...
    'gH_max', zeros(max_popul,1), 'K_HR', zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
% structure of a Newborn community
newborn_struct = struct('M_L',zeros(nb_popul,1),'H_L',zeros(nb_popul,1),'fp',zeros(nb_popul,1),...
    'gM_max', zeros(nb_popul,1), 'K_MB', zeros(nb_popul,1), 'K_MR', zeros(nb_popul,1),...
    'gH_max', zeros(nb_popul,1), 'K_HR', zeros(nb_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));
% structure of simulation constants
const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target, 't_bin', t_bin,...
    'gM_max_Bound', gM_max_Bound, 'gH_max_Bound', gH_max_Bound,...
    'fp_Bound', fp_Bound, 'K_MB_Bound', K_MB_Bound,'K_MR_Bound', K_MR_Bound, ...
    'K_HR_Bound', K_HR_Bound, 'R0', R0, 'K_singular', K_singular, 'dM', dM, 'dH', dH, ...
    'c_BM', c_BM, 'c_RM', c_RM, 'c_RH', c_RH, 'mut_rate', mut_rate);
%%
% rng('shuffle');
load('C778/distrng.mat')
rng(distrng)
% if the simulation continues from a previous one, pick up
if C_prev > 0
    load(['C' num2str(C_prev+1) '/comm_all/newborns'])
    % H_isolates_in stores H cells isolated from C_prev cycle needed for spiking
    load(['C' num2str(C_prev) '/H_isolates_in'])
    load(['C' num2str(C_prev) '/M_isolates_in'])
    % check_cycle_m stores all the number of cycles when heritability was checked
    if exist('comm_all/check_cycle_m.mat', 'file')
        load('comm_all/check_cycle_m')
        check_counter = length(check_cycle_m);
        if check_counter ==0
            spike_frac = 0; % fraction of H biomass to be spiked in current scheme
            % test the following spiking fractions of H biomass
            spike_test = spike_initial(2:end);
            spike_all = [spike_frac spike_test];
            check_cycle = C_prev+check_period;
        else
            % heritability is checked during the check_cycle
            check_cycle = max(C_prev+3, check_cycle_m(end)+check_period);
            %             load(['C' num2str(check_cycle_m(end)) '/spike_all']);
            load('comm_all/spike_all')
            spike_frac = spike_all(1);
            spike_test = spike_all(2:end);
        end
    else
        spike_frac = 0; % fraction of H biomass to be spiked in current scheme
        % test the following spiking fractions of H biomass
        spike_test = spike_initial(2:end);
        spike_all = [spike_frac spike_test];
        check_cycle_m = [];
        check_cycle = C_prev+check_period;
        check_counter = 0;
    end
    sl = length(spike_test); % number of spiking fraction to be tested
else
    % initialize the simulation
    fp_init = zeros(nb_popul, 1);
    gM_max_init = zeros(nb_popul, 1);
    K_MB_init = zeros(nb_popul, 1);
    K_MR_init = zeros(nb_popul, 1);
    gH_max_init = zeros(nb_popul, 1);
    M_L_init = zeros(nb_popul, 1);
    K_HR_init = zeros(nb_popul, 1);
    H_L_init=zeros(nb_popul, 1);
    % death_probability=zeros(max_popul,1);
    M_counter = 60 ;
    H_counter = BM_target - M_counter;
    M_L_init(1 : M_counter) = 1;
    H_L_init(1 : H_counter) = 1;
    fp_init(1 : M_counter) = fp_start;
    gM_max_init(1 : M_counter) = gM_max_start;
    K_MB_init(1 : M_counter) = K_MB_start;
    K_MR_init(1 : M_counter) = K_MR_start;
    gH_max_init(1 : H_counter)=gH_max_start;
    K_HR_init(1 : H_counter) = K_HR_start;
    
    newborns(1, 1:N)...
        =struct('M_L', M_L_init,'H_L', H_L_init,'fp', fp_init,...
        'gM_max', gM_max_init, 'K_MB', K_MB_init, 'K_MR', K_MR_init,...
        'gH_max', gH_max_init, 'K_HR', K_HR_init,...
        'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
        'P',0,'parentnum',1,'rseed',uint32(0));
    rseed=randi(2^32-1,N,1,'uint32');
    for ri=1:N
        newborns(ri).rseed=rseed(ri);
    end
    if ~exist('C1/comm_all', 'dir')
        mkdir('C1/comm_all')
    end
    save('C1/comm_all/newborns', 'newborns')
    H_isolates_in = struct('gH_max', gH_max_start*ones(spike_clone_num, 1),...
        'K_HR', K_HR_start*ones(spike_clone_num, 1),...
        'H_L', 1/log(2)*ones(spike_clone_num, 1));
    M_isolates_in = struct('gM_max', gM_max_start*ones(spike_clone_num, 1),...
        'K_MR', K_MR_start*ones(spike_clone_num, 1), 'M_L', 1/log(2)*ones(spike_clone_num, 1),...
        'fp', fp_start*ones(spike_clone_num, 1), 'K_MB', K_MB_start*ones(spike_clone_num, 1));
    spike_frac = 0; % fraction of H pure culture spiked in
    spike_test = spike_initial(2:end);
    sl = length(spike_test);
    spike_all = [spike_frac spike_test];
    check_cycle_m = [];
    check_cycle = check_period;
    check_counter = 0;
end
% randomize the random number seed
n = C_prev+1;
check_flag = false;
% flag_m = zeros(C,2);
while n <= C
    %     flag_m(n, :) = [n uint8(check_flag)];
    % create a folder Cn to save the results of the nth cycle
    folder_name1 = ['C' num2str(n)];
    if ~exist(folder_name1, 'dir')
        mkdir(folder_name1)
    end
    folder_name2=['C' num2str(n) '/comm_all'];
    if ~exist(folder_name2, 'dir')
        mkdir(folder_name2)
    end
    comm_all(1,1:N) = comm_struct;
    % rep is the index of communities within one cycle
    parfor rep = 1:N
        if sum(newborns(rep).M_L) + sum(newborns(rep).H_L) < pcs
            error('An Newborn is empty')
        end
        comm_all(rep) = simu_one_comm(newborns(rep), comm_struct, const_struct);
    end
    newborns(1, 1:N) = newborn_struct;
    distrng=rng;
    if Pn_sig>pcs
        % Pn is the measurement noise in P(T)
        Pn = normrnd(0, Pn_sig, size([comm_all.P]));
        save([folder_name2 '/Pn'],'Pn');
    else
        Pn = 0;
    end
    P_all = [comm_all.P];
    save([folder_name2 '/P_all'], 'P_all');
    [~, I] = sort([comm_all.P]+Pn, 'descend');
    if max(P_all) < pcs
        break
    end
    % I=randperm(N);
    % in cycle check_cycle-2, parent communities are generated if there are enough cells
    if abs(n - check_cycle + 2) < 0.1 && check_cycle <= C
        sel_counter = 0;
        rep_counter = 0;
        % rep_num_M is a matrix for the numbers of newborns for selection
        % and tests for heritability under different spiking fraction
        rep_num_M = zeros(N, sl+1); % N = 100, sl+1 = 3, # of substitution conditions;
        for i = 1 : N
            if rep_counter >= N
                break
            end
            % BM is the biomass of the Adult
            BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
            % test_rep_temp is the total number of test comms that the Adult can
            % contribute alternative (not current) substitution conditions for
            % heritability checks
            test_rep_temp = floor((BM - comm_rep_num*BM_target*(1-abs(spike_frac)))...
                /BM_target/sum(1-abs(spike_test)));
            % if the number of test comms is too small, postpone the test cycle by 1
            if test_rep_temp < round(comm_rep_num/5)
                check_cycle = check_cycle+10;
                check_flag = false;
                break
            else
                rep_num_M(i, 1) = comm_rep_num;
                rep_num_M(i, 2:end) = min(test_rep_temp, comm_rep_num);
                check_flag = true;
            end
            sel_counter = sel_counter+1;
            rep_counter = rep_counter+comm_rep_num;
        end
        % if there are at least N/5 comms to test heritability
        if check_flag == true
            comm_selected = comm_all(I(1:sel_counter));
            % rem_num_M has sel_counter rows for sel_counter chosen Adults
            rep_num_M = rep_num_M(1:sel_counter, :);
            test_rep_total = sum(rep_num_M(:, end));
            % newborns_par is the structure array for parent Newborns
            newborns_par(1:sl, 1:test_rep_total) = newborn_struct;
            % H isolates that might be used to spike during the next cycle
            H_isolates_out(1, 1:sel_counter) ...
                = struct('gH_max', -ones(spike_clone_num, 1), 'K_HR', -ones(spike_clone_num, 1),...
                'H_L', -ones(spike_clone_num, 1));
            % M isolates that might be used to spike during the next cycle
            M_isolates_out(1, 1:sel_counter) ...
                = struct('gM_max', -ones(spike_clone_num, 1), 'K_MR', -ones(spike_clone_num, 1),...
                'M_L', -ones(spike_clone_num, 1), 'fp', -ones(spike_clone_num, 1),...
                'K_MB', -ones(spike_clone_num, 1));
            rseed = randi(2^32-1, sl, test_rep_total, 'uint32');
            comm_par_counter = 0;
            rep_counter = 0;
            for i = 1:sel_counter
                if abs(rep_num_M(i, 1) - comm_rep_num) > pcs
                    error('possible errors in determining newborn rep num')
                end
                % reproduce the ith Adult into rep_num_M(i, 1) Newborns for selection,
                % rep_num_M(i, 2) Newborns for testing heritability when spiking
                % spike_test(1) H biomass, rep_num_M(i, 3) Newborns for testing
                % heritability when spiking spike_test(2) H biomass, and so on
                % comm_temp is a row structure array with sum(rep_num_M(i, :)) elements
                [comm_temp, H_isolates_out(i), M_isolates_out(i)] = pipette_SpikeMix_SPHM(comm_selected(i), newborn_struct, ...
                    const_struct, [spike_frac spike_test], rep_num_M(i, :), ...
                    H_isolates_in(comm_selected(i).parentnum), M_isolates_in(comm_selected(i).parentnum),...
                    spike_clone_num, i);
                newborns(rep_counter+1:rep_counter+rep_num_M(i, 1)) ...
                    = comm_temp(1:rep_num_M(i, 1));
                rep_counter = rep_counter + rep_num_M(i, 1);
                % newborns_par(j, :) contains parent Newborns with jth substitution fraction
                for j = 1:sl
                    newborns_par(j, comm_par_counter+1:comm_par_counter+rep_num_M(i,end)) ...
                        = comm_temp(comm_rep_num+1+(j-1)*rep_num_M(i,end) : comm_rep_num+j*rep_num_M(i,end));
                end
                comm_par_counter = comm_par_counter+rep_num_M(i,end);
            end
            %             % If an Adult doesn't have H to be isolated for spiking for the next cycle,
            %             % fill the spot with H isolates from the Adult with the highest P(T). If none
            %             % of have H to be isolated for spiking, use those H isolates from the previous
            %             % cycle.
            %             H_isolates_out = clean_HM_isolates(H_isolates_out, H_isolates_in);
            H_isolates_in = H_isolates_out;
            %             M_isolates_out = clean_HM_isolates(M_isolates_out, M_isolates_in);
            M_isolates_in = M_isolates_out;
            % assign random number seeds to each parent Newborn.
            for i = 1:sl
                for j = 1:test_rep_total
                    newborns_par(i,j).rseed = rseed(i,j);
                end
            end
        end
        % in cycle check_cycle-1, offspring communities are generated if there are enough cells
    elseif abs(n - check_cycle + 1) < 0.1
        if check_flag==0
            error('n = %d, check_flag=0, inconsistent', n)
        end
        test_rep = min(test_rep_total, test_rep_max);
        % H isolates that might be used to spike during the next cycle
        sel_counter = 0;
        rep_counter = 0;
        rep_num_C = zeros(N, 1);
        heri_par_idx = [];
        heri_par_counter = 0;
        % Calculate how many Adults are required to check the heritability. If more than
        % test_rep are needed, postpond heritability check by 5 cycles.
        for i = 1 : N
            if rep_counter >= N
                break
            end
            BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
            dil_factor = floor(BM/BM_target/(1-abs(spike_frac)));
            if dil_factor == 0
                continue
            end
            rep_num_temp = min([dil_factor-off_rep_max, comm_rep_num, N-rep_counter]);
            % if a chosen Adult produces less than test_rep_num Newborns, skip this Adult
            % for heritability check
            if rep_num_temp > 0
                heri_par_counter = heri_par_counter+1;
                heri_par_idx(heri_par_counter) = i;
                rep_num_C(i) = rep_num_temp+off_rep_max;
            else
                rep_num_temp = min([dil_factor, comm_rep_num, N-rep_counter]);
                rep_num_C(i) = rep_num_temp;
            end
            sel_counter = sel_counter+1;
            rep_counter = rep_counter+rep_num_temp;
        end
        % if more than test_rep Adults are required to check heritability, postpond the
        % check by 10 cycles
        if sel_counter > test_rep
            check_cycle = check_cycle+10;
            check_flag = false;
            newborns(1, 1:N) = newborn_struct;
            rep_num_C = zeros(N, 1);
        else
            check_counter = check_counter+1;
            check_cycle_m(check_counter) = check_cycle;
            save('check_cycle_m','check_cycle_m')
            H_isolates_out(1, 1:sel_counter) ...
                = struct('gH_max', -ones(spike_clone_num, 1), 'K_HR', -ones(spike_clone_num, 1),...
                'H_L', -ones(spike_clone_num, 1));
            % M isolates that might be used to spike during the next cycle
            M_isolates_out(1, 1:sel_counter) ...
                = struct('gM_max', -ones(spike_clone_num, 1), 'K_MR', -ones(spike_clone_num, 1),...
                'M_L', -ones(spike_clone_num, 1), 'fp', -ones(spike_clone_num, 1),...
                'K_MB', -ones(spike_clone_num, 1));
            % newborns_off is the structure array for offspring Newborns
            newborns_off(1:off_rep_max*(sl+1), 1:test_rep) = newborn_struct;
            comm_selected = comm_all(I(1:sel_counter));
            rep_num_C = rep_num_C(1:sel_counter);
            % P_par is the array for P(T) of parents
            P_par = zeros(size(newborns_par));
            P_par_sorted = zeros(size(newborns_par));
            M0_par = zeros(size(newborns_par));
            H0_par = zeros(size(newborns_par));
            % adults_par is the structure array for parent Adults
            adults_par(1:sl, 1:test_rep_total) = comm_struct;
            parfor i = 1:test_rep_total*sl
                adults_par(i) = simu_one_comm(newborns_par(i), comm_struct, const_struct);
                P_par(i) = adults_par(i).P;
                M0_par(i) = adults_par(i).M_t(1);
                H0_par(i) = adults_par(i).H_t(1);
            end
            save([folder_name1 '/ParResults'],'P_par','M0_par','H0_par','spike_all')
            if Pn_sig>pcs
                Pn_par = normrnd(0, Pn_sig, size(P_par));
                save([folder_name1 '/ParResults'], 'Pn_par', '-append')
                P_par = P_par + Pn_par;
            end
            % sort the parent Adults according to their functions
            for i = 1:sl
                [P_par_sorted(i,:), I_par] = sort(P_par(i,:), 'descend');
                adults_par(i,:) = adults_par(i, I_par);
            end
            % only the top test_rep parent Adults reproduce
            adults_par = adults_par(:,1:test_rep);
            rseed = randi(2^32-1,off_rep_max*(sl+1), test_rep, 'uint32');
            % reproduce Adults of the current substitution ratio into Newborns for
            % selection and offspring Newborns for checking heritability
            rep_counter = 0;
            for i = 1 : sel_counter
                BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
                dil_factor = floor(BM/BM_target/(1-abs(spike_frac)));
                if dil_factor == 0
                    continue
                end
                rep_num_temp = min([dil_factor-off_rep_max, comm_rep_num, N-rep_counter]);
                if rep_num_temp > 0
                    [comm_temp, H_isolates_out(i), M_isolates_out(i)] = pipette_SpikeMix_SPHM(comm_all(I(i)), newborn_struct, ...
                        const_struct, spike_frac, rep_num_C(i), ...
                        H_isolates_in(comm_all(I(i)).parentnum), ...
                        M_isolates_in(comm_all(I(i)).parentnum), spike_clone_num, i);
                    newborns(rep_counter+1:rep_counter+rep_num_temp) ...
                        = comm_temp(1:rep_num_temp);
                    newborns_off(1 : off_rep_max, i) ...
                        = transpose(comm_temp(rep_num_temp+1:end));
                else
                    rep_num_temp = min([dil_factor, comm_rep_num, N-rep_counter]);
                    [newborns(rep_counter+1:rep_counter+rep_num_temp), H_isolates_out(i), M_isolates_out(i)]...
                        = pipette_SpikeMix_SPHM(comm_all(I(i)), newborn_struct, ...
                        const_struct, spike_frac, rep_num_C(i), ...
                        H_isolates_in(comm_all(I(i)).parentnum), M_isolates_in(comm_all(I(i)).parentnum), ...
                        spike_clone_num, i);
                end
                rep_counter = rep_counter+rep_num_temp;
            end
            for i = sel_counter+1:test_rep
                BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
                dil_factor = floor(BM/BM_target/(1-abs(spike_frac)));
                off_rep_num = min(off_rep_max, dil_factor);
                if dil_factor > 0
                    [comm_temp, ~, ~]...
                        = pipette_SpikeMix_SPHM(comm_all(I(i)), newborn_struct, ...
                        const_struct, spike_frac, off_rep_num, H_isolates_in(comm_all(I(i)).parentnum),...
                        M_isolates_in(comm_all(I(i)).parentnum), spike_clone_num, i);
                    newborns_off(1 : off_rep_num, i) = transpose(comm_temp);
                    heri_par_counter = heri_par_counter+1;
                    heri_par_idx(heri_par_counter) = i;
                end
            end
            % reproduce each parent Adult into up to test_rep_num offspring Newborns
            for i = 1:sl % i is the index for substitution fraction
                for j = 1:test_rep
                    BM = adults_par(i,j).M_t(end) + adults_par(i,j).H_t(end);
                    dil_factor = floor(BM/BM_target/(1-abs(spike_test(i))));
                    off_rep_num = min(off_rep_max, dil_factor);
                    [comm_temp, ~, ~] = pipette_SpikeMix_SPHM(adults_par(i,j), newborn_struct,...
                        const_struct, spike_test(i), off_rep_num,...
                        H_isolates_in(adults_par(i,j).parentnum),...
                        M_isolates_in(adults_par(i,j).parentnum),...
                        spike_clone_num, j);
                    newborns_off(i*off_rep_max+1:i*off_rep_max+off_rep_num, j) ...
                        = transpose(comm_temp);
                end
            end
            %             H_isolates_out = clean_HM_isolates(H_isolates_out, H_isolates_in);
            H_isolates_in = H_isolates_out;
            %             M_isolates_out = clean_HM_isolates(M_isolates_out, M_isolates_in);
            M_isolates_in = M_isolates_out;
            for i = 1:off_rep_max * (sl+1) * test_rep
                newborns_off(i).rseed = rseed(i);
            end
            clear newborns_par adults_par
        end
    elseif abs(n - check_cycle) < 0.1
        if check_flag == 0
            error('n = %d, check_flag=0, inconsistent', n)
        end
        % P_off is the array for P(T) of each offspring Adult
        P_off = zeros(size(newborns_off));
        M0_off = zeros(size(newborns_off));
        H0_off = zeros(size(newborns_off));
        parfor i=1 : test_rep*off_rep_max*(sl+1)
            if sum(newborns_off(i).M_L)+sum(newborns_off(i).H_L) < pcs
                P_off(i) = NaN;
                M0_off(i) = 0;
                H0_off(i) = 0;
            else
                comm_temp = simu_one_comm(newborns_off(i), comm_struct, const_struct);
                P_off(i) = comm_temp.P;
                if P_off(i) <0
                    error('negative P_off')
                end
                M0_off(i) = comm_temp.M_t(1);
                H0_off(i) = comm_temp.H_t(1);
            end
        end
        save([folder_name1 '/OffResults'],'P_off','M0_off','H0_off')
        clear newborns_off
        load(['C' num2str(check_cycle-1) '/comm_all/P_all'])
        if Pn_sig>pcs
            Pn_off = normrnd(0, Pn_sig, size(P_off));
            save([folder_name1 '/OffResults'], 'Pn_off', '-append')
            P_off = P_off + Pn_off;
            load(['C' num2str(check_cycle-1) '/comm_all/Pn'])
            P_sorted = sort(P_all+Pn, 'descend');
        else
            P_sorted = sort(P_all, 'descend');
        end
        check_flag = 0;
        % heri stores the heritability under each substitution fraction
        heri = zeros(sl+1, 1);
        % ub stores the upper confidence interval under each substitution fraction
        ub = zeros(sl+1, 1);
        % lb stores the lower confidence interval under each substitution fraction
        lb = zeros(sl+1, 1);
        % heri stores the p_value of heritability under each substitution fraction
        %         p_val = zeros(sl+1, 1);
        % heritability is defined as the Spearman correlation coefficient between P(T) of
        % parent Adults and average P(T) among offspring Adults from each lineage
        heri(1) = slope_func((P_sorted(heri_par_idx))',...
            (nanmedian(P_off(1:off_rep_max, heri_par_idx)))');
        % the confidence interval is estimated from bootstraping
        [lb(1), ub(1)] = bstrap_itvl((P_sorted(heri_par_idx))',...
            (nanmedian(P_off(1:off_rep_max, heri_par_idx)))', @slope_func, n_bstrap, q);
        for i = 1:sl
            heri(i+1) = slope_func(P_par_sorted(i, 1:test_rep)', ...
                (nanmedian(P_off(i*off_rep_max+1:(i+1)*off_rep_max, 1:test_rep)))');
            [lb(i+1), ub(i+1)] = bstrap_itvl(P_par_sorted(i, 1:test_rep)',...
                (nanmedian(P_off(i*off_rep_max+1:(i+1)*off_rep_max, 1:test_rep)))', @slope_func, n_bstrap, q);
        end
        clear P_par_sorted
        
        if HeriSwitch == 1
            [heri_sorted, I_heri] = sort(heri(2:end),'descend');
            spike_idx = [];
            for i = 1:sl
                % if the heritability of alternative spiking ratio is larger than
                % the 95% confidence interval, switch.
                if heri_sorted(i) > ub(1)
                    
                    %             % if the lb of the heritability of alternative spiking ratio is larger than
                    %             % the ub of the heritability of the current spiking ratio, switch.
                    %             if lb(I_heri(i)+1) > ub(1)
                    spike_idx = I_heri(i)+1;
                    break
                end
            end
            if ~isempty(spike_idx)
                spike_frac = spike_all(spike_idx);
                spike_all(spike_idx) = [];
                spike_test = spike_all;
                spike_all = [spike_frac spike_test];
            end
        elseif HeriSwitch == -1
            [heri_sorted, I_heri] = sort(heri(2:end),'ascend');
            spike_idx = [];
            for i = 1:sl
                if heri_sorted(i) < lb(1)
                    spike_idx = I_heri(i)+1;
                    break
                end
            end
            if ~isempty(spike_idx)
                spike_frac = spike_all(spike_idx);
                spike_all(spike_idx) = [];
                spike_test = spike_all;
                spike_all = [spike_frac spike_test];
            end
        elseif HeriSwitch == 0
            spike_all = spike_all(randperm(sl+1));
            spike_frac = spike_all(1);
            spike_test = spike_all(2:end);
        else
            error('HeriSwitch value not valid')
        end
        
        save([folder_name1 '/OffResults'],'heri','lb','ub','heri_par_idx','-append')
        save([folder_name1 '/spike_all'], 'spike_all')
        clear heri_par_idx
        if check_cycle + check_period <= C
            check_cycle = check_cycle + check_period;
        end
    end
    % Reproducing chosen Adults when heritability is not checked
    if check_flag ==0
        if exist('comm_selected', 'var')>0 || max([newborns.rseed])>0
            error('Communities already reproduced')
        end
        rep_counter = 0;
        sel_counter = 0;
        % H isolates that might be used to spike during the next cycle
        H_isolates_out(1, 1:N) ...
            = struct('gH_max', -ones(spike_clone_num, 1), 'K_HR', -ones(spike_clone_num, 1),...
            'H_L', -ones(spike_clone_num, 1));
        % M isolates that might be used to spike during the next cycle
        M_isolates_out(1, 1:N) ...
            = struct('gM_max', -ones(spike_clone_num, 1), 'K_MR', -ones(spike_clone_num, 1),...
            'M_L', -ones(spike_clone_num, 1), 'fp', -ones(spike_clone_num, 1),...
            'K_MB', -ones(spike_clone_num, 1));
        for i = 1 : N
            if rep_counter >= N
                break
            end
            dil_factor = floor((comm_all(I(i)).M_t(end)+comm_all(I(i)).H_t(end))/BM_target/(1-abs(spike_frac)));
            if dil_factor == 0
                continue
            end
            rep_num_temp = min([dil_factor, comm_rep_num, N-rep_counter]);
            [newborns(rep_counter+1:rep_counter+rep_num_temp), H_isolates_out(i), M_isolates_out(i)]...
                = pipette_SpikeMix_SPHM(comm_all(I(i)), newborn_struct, ...
                const_struct, spike_frac, rep_num_temp, H_isolates_in(comm_all(I(i)).parentnum),...
                M_isolates_in(comm_all(I(i)).parentnum), spike_clone_num, i);
            sel_counter=sel_counter+1;
            rep_counter=rep_counter+rep_num_temp;
        end
        %         H_isolates_out = clean_HM_isolates(H_isolates_out, H_isolates_in);
        H_isolates_in = H_isolates_out(1:sel_counter);
        %         M_isolates_out = clean_HM_isolates(M_isolates_out, M_isolates_in);
        M_isolates_in = M_isolates_out(1:sel_counter);
        comm_selected = comm_all(I(1:sel_counter));
    end
    % assign random number seed to Newborns of the next cycle
    rseed = randi(2^32-1, N, 1, 'uint32');
    for ri = 1 : N
        newborns(ri).rseed = rseed(ri);
    end
    clear H_isolates_out M_isolates_out
    % save the selected communities and the current state of the random number generator
    save([folder_name1 '/comm_selected'],'comm_selected');
    save([folder_name1 '/H_isolates_in'],'H_isolates_in');
    save([folder_name1 '/M_isolates_in'],'M_isolates_in');
%     save([folder_name1 '/distrng'],'distrng');
    clear comm_selected
    n = n+1;
    % save Newborns for the next cycle
    folder_name3=['C' num2str(n) '/comm_all'];
    if ~exist(folder_name3, 'dir')
        mkdir(folder_name3)
    end
    save([folder_name3 '/newborns'],'newborns');
end
% save Newborns of the (C+1)th cycle
if ~exist('comm_all','dir')
    mkdir('comm_all')
end
% save('comm_all/newborns', 'newborns')
% save('comm_all/check_cycle_m','check_cycle_m')
% save('comm_all/spike_all','spike_all')
% save('comm_all/flag_m','flag_m')
