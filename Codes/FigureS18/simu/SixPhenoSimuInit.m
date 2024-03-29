function [spike_all, check_cycle_prev, H_isolates_in, M_isolates_in, newborns,...
     P_sel_dynamics, check_cycle_m]...
    = SixPhenoSimuInit(C_prev)
load('SixPhenoPara')
% if the simulation continues from a previous one, pick up
if C_prev > 0
    load(['C' num2str(C_prev+1) '/comm_all/newborns'])
    % H_isolates_in stores H cells isolated from C_prev cycle needed for spiking
    load(['C' num2str(C_prev) '/H_isolates_in'])
    load(['C' num2str(C_prev) '/M_isolates_in'])
    load('comm_all/P_sel_dynamics')
    % check_cycle_m stores all the number of cycles when heritability was checked
    if exist('comm_all/check_cycle_m.mat', 'file')
        load('comm_all/check_cycle_m')
        check_counter = length(check_cycle_m);
        if check_counter ==0
            spike_all = spike_initial;
            check_cycle_prev = 0;
        else
            %             load(['C' num2str(check_cycle_m(end)) '/spike_all']);
            load('comm_all/spike_all')
            check_cycle_prev = check_cycle_m(end);
        end
    else
        spike_all = spike_initial;
        check_cycle_m = [];
        check_cycle_prev = 0;
    end
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
    if isfolder('C1')
        load('C1/comm_all/newborns.mat')
        rseed = [newborns.rseed];
    else
        rng('shuffle')
        rseed=randi(2^32-1,N,1,'uint32');
    end
    for ri=1:N
        newborns(ri).rseed=rseed(ri);
    end
    if ~exist('comm_all','dir')
        mkdir('comm_all')
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
    spike_all = spike_initial;
    check_cycle_m = [];
    P_sel_dynamics = [];
    check_cycle_prev = 0;
end