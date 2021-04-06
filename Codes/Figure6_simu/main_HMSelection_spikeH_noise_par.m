% Simulating community selection of H-M community in the simple scenario
% where only fp mutates. Dynamicsare plotted in Figure 6 and Figure 3A.
clear
C = 3; % total number of cycles
% comm_type_num = n for the top n% strategy.
comm_type_num = 2; 
Pn_sig = 100; % std of the measurement noise. In case of no noise, Pn_sig = 0. 
% reproducing method
repro_method = @fixBM0_spike;
% repro_method = @cell_sort; % reproduce through cell sorting, no spiking
spike_frac = 0; % fraction of H pure culture spiked in

T0 = 17; % Maturation time
mut_rate = 1e-2; % mut_rate = 1e-2 corresponding to effective mutation rate of 2e-3,
N = 100; % number of communities within a cycle
% comm_type_num * comm_rep_num = number of communities within one cycle.
comm_rep_num = N/comm_type_num; % maximal number of offspring community from one Adult.
max_popul = 1e4; % maximal number of cells in the community
t_bin = 0.05; % time step in the simulation
pcs=1e-15; % precision constant
t_binnum = int16(T0/t_bin); % number of time steps


% BM_target is the target biomass, T0 is the maturation time
BM_target = 100;

% R(0)=1, 
R0 = 1;
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
fp_max = 1; % fp is between 0 and 1

% initialize for the first cycle
M0=30; % number of M cells in a Newborn
H0=BM_target-M0; % number of H cells in a Newborn

% structure for communities
comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
    'M_t',zeros(t_binnum+1,1),'H_t',zeros(t_binnum+1,1),'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),...
    'P',0,'parentnum',0,'rseed',uint32(0));

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target);

% Initialize the first cycle
% column vector that stores fp of each M cell
fp = zeros(max_popul,1);
% column vector that stores biomass of each M cell
M_L = zeros(max_popul,1);
% column vector that stores biomass of each H cell
H_L = zeros(max_popul,1);
% number of M cells
M_counter = M0;
% number of H cells
H_counter = H0;

M_L(1 : M0) = 1;
H_L(1 : H0) = 1;
fp(1 : M0) = fp_start;
M_t=zeros(t_binnum+1, 1);
H_t=zeros(t_binnum+1, 1);
comm_all(1:comm_type_num*comm_rep_num,1)=struct('M_L',M_L,'H_L',H_L,'fp',fp,'M_t',M_t,...
    'H_t',H_t,'R',zeros(t_binnum+1,1),'B',zeros(t_binnum+1,1),'P',0,'parentnum',0,'rseed',uint32(0));

% shuffle the random number state
rng('shuffle');
rseed=randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');

for ri=1:comm_type_num*comm_rep_num
    comm_all(ri).rseed=rseed(ri);
end


% perform a total C cycles
for n = 1 : C
    % create a folder Cn to save the results of the nth cycle
    folder_name1 = ['C' num2str(n)];
    if ~exist(folder_name1,'dir')
        mkdir(folder_name1)
    end
    newborns = comm_all;
    % rep is the index of communities within one cycle
    parfor rep = 1 : comm_type_num*comm_rep_num
        % seed the random number generator
        rnodeseed = comm_all(rep).rseed;
        rng(rnodeseed,'twister');
        comm_rep = comm_struct;
        % copy the data from the structure with Newborn's configuration
        M_L = comm_all(rep).M_L;
        H_L = comm_all(rep).H_L;
        fp = comm_all(rep).fp;
        M_t = zeros(t_binnum+1, 1);
        H_t = zeros(t_binnum+1, 1);
        %  temporary column vectors for M and H's biomass
        M_L_Temp = zeros(max_popul,1);
        H_L_Temp = zeros(max_popul,1);
        % column vectors to store B and R's concentrations at each time
        % step
        B = zeros(t_binnum+1, 1);
        R = zeros(t_binnum+1, 1);
        P = 0;
        M_counter = nnz(M_L>pcs);
        H_counter = nnz(H_L>pcs);
        M_t(1) = sum(M_L);
        H_t(1) = sum(H_L);
        R(1) = R0;
        % perform simulation through t_binnum time steps
        for dt = 2 : t_binnum+1
            % solve for the dynamics of B and R within the time step
            para = [gM_max; gH_max; c_BM; c_RM; c_RH; K_MR; K_HR; K_MB; sum(M_L); sum(H_L)];
            fhandle = @(t, y) chem_conc_jacobian(t, y, para);
            options = odeset('Jacobian',fhandle,'RelTol',1e-5);
            [tx,y] = ode23s(@(t,y) chem_conc(t,y,para),[0 t_bin],[B(dt-1); R(dt-1)], options);
            if ~isreal(y)
                error('imaginary value')
            end
           
            BN = y(:,1)/K_MB;
            RN = y(:,2)/K_MR;
            KHKM = K_HR/K_MR;
            % calculate the instantaneous growth rates for M and H
            g_M = BN./(BN+RN).*RN./(RN+1)+RN./(BN+RN).*BN./(BN+1);
            g_H = RN./(RN+KHKM);
            
            % calculate the biomass for each M and H cells at the end of
            % the time step
            gM_exp = trapz(tx,g_M)*gM_max;
            gH_exp = trapz(tx,g_H)*gH_max;
            M_L_Temp(1:M_counter) = exp((1-fp(1:M_counter))*gM_exp).*M_L(1:M_counter);
            H_L_Temp(1:H_counter) = exp(gH_exp)*H_L(1:H_counter);
            % calculate the amount of Product at the end of the time step
            P = sum(fp(1:M_counter).*(M_L_Temp(1:M_counter)-M_L(1:M_counter))./(1-fp(1:M_counter)))+P;
            % update the biomass for M and H cells
            M_L(1:M_counter) = M_L_Temp(1:M_counter);
            H_L(1:H_counter) = H_L_Temp(1:H_counter);
            % register the concentration of B and R at the end of the time
            % step
            B(dt) = BN(end)*K_MB;
            R(dt) = RN(end)*K_MR;
            
            % calculate the death probability for each M cell and eliminate
            % those that die within this time step, update the total
            % biomass of M
            death_probability = rand(max_popul,1);
            M_L(death_probability < dM * t_bin) = 0;
            fp( death_probability < dM * t_bin) = 0;
            M_t(dt) = sum(M_L);
            
            % calculate the death probability for each H cell and eliminate
            % those that die within this time step, update the total
            % biomass of H
            death_probability = rand(max_popul,1);
            H_L( death_probability < dH * t_bin) = 0;
            H_t(dt) = sum(H_L);
            
            % if the biomass of a M cell is greater than 2, it divides
            div_idx = find(M_L >= 2);
            div_length = length(div_idx);
            % draw mutation effect from mu_spec, and change the value of fp
            % of both daughter cells
            if div_length > 0
                mut_multiplier1 = mu_spec(div_length).*double(rand(div_length,1) <= mut_rate)+1;
                mut_multiplier2 = mu_spec(div_length).*double(rand(div_length,1) <= mut_rate)+1;
                fp(M_counter + 1 : M_counter + div_length) = fp(div_idx) .* mut_multiplier1;
                fp(div_idx) = fp(div_idx) .* mut_multiplier2;
                M_L(M_counter + 1 : M_counter + div_length) = M_L(div_idx) / 2;
                M_L(div_idx) = M_L(div_idx) / 2;
                M_counter = M_counter+div_length;
            end
            fp(fp>fp_max) = fp_max;
            
            % if the biomass of a H cell is greater than 2, it divides
            div_idx = find(H_L >= 2);
            div_length = length(div_idx);
            if div_length>0
                H_L(H_counter+1:H_counter+div_length) = H_L(div_idx)/2;
                H_L(div_idx) = H_L(div_idx)/2;
                H_counter = H_counter+div_length;
            end
        end
        % eliminate dead cells and update the list of fp and biomass for
        % both M and H cells
        M_counter = nnz(M_L > pcs);
        comm_rep.M_L(1 : M_counter) = M_L(M_L > pcs);
        H_counter = nnz(H_L > pcs);
        comm_rep.H_L(1 : H_counter) = H_L(H_L > pcs);
        comm_rep.fp(1 : M_counter) = fp(M_L > pcs);
        comm_rep.M_t = M_t;
        comm_rep.H_t = H_t;
        comm_rep.B = B;
        comm_rep.R = R;
        comm_rep.P = P;
        comm_rep.parentnum = comm_all(rep).parentnum;
        comm_rep.rseed = comm_all(rep).rseed;
        % return the result to the structure array
        comm_all(rep) = comm_rep;
    end
    % create a folder to store the data for all communities
    folder_name2 = ['C' num2str(n) '/comm_all'];
    if ~exist(folder_name2,'dir')
        mkdir(folder_name2)
    end
    adults = comm_all;
%     save([folder_name2 '/adults'],'adults');
    save([folder_name2 '/newborns'],'newborns');
    % save the current state of the random number generator
    distrng = rng;
    % add stochastic noise to the community function
    if Pn_sig > pcs
        Pn = normrnd(0, Pn_sig, size([comm_all.P]));
    else
        Pn = 0;
    end
    % sort the communities in descending order of their functions
    [~, I] = sort([comm_all.P] + Pn,'descend');
    % I = randperm(comm_type_num*comm_rep_num);
    comm_all_sorted = comm_all(I);
    % counter for the Newborn in the next cycle
    rep_counter = 0;
    % counter for the selected Adults
    sel_counter = 0;
    comm_selected(1:comm_type_num*comm_rep_num,1) = comm_struct;
    % reproduce chosen Adult communities.
    for i = 1 : comm_type_num*comm_rep_num
        if rep_counter >= comm_type_num*comm_rep_num
            break
        end
        dil_factor = floor((comm_all_sorted(i).M_t(t_binnum+1)+comm_all_sorted(i).H_t(t_binnum+1))/BM_target/(1-spike_frac));
        if dil_factor == 0
            continue
        end
        rep_num_temp = min(dil_factor,comm_rep_num);
        comm_all_idx = min(comm_type_num*comm_rep_num,rep_counter+rep_num_temp);
        comm_all(rep_counter+1:comm_all_idx) = repro_method(comm_all_sorted(i), comm_struct, const_struct, dil_factor, BM_target*spike_frac, rep_counter, i);
        sel_counter = sel_counter+1;
        comm_selected(sel_counter) = comm_all_sorted(i);
        rep_counter = rep_counter+rep_num_temp;
    end
    % the selected Adult communities
    comm_selected = comm_selected(1:sel_counter);
    % assign random number seed to Newborns of the next cycle
    rseed = randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
    for ri = 1 : comm_type_num*comm_rep_num
        comm_all(ri).rseed = rseed(ri);
    end
    % save the selected communities and the current state of the random number generator 
    save([folder_name1 '/comm_selected'],'comm_selected');
    save([folder_name1 '/distrng'],'distrng');
    if Pn_sig > pcs
        save([folder_name1 '/Pn'],'Pn');
    end
end
% % save Newborns of the (C+1)th cycle 
% if ~exist('comm_all','dir')
%     mkdir('comm_all')
% end
% for i = 1:comm_type_num
%     comm = comm_all((i-1)*comm_rep_num+1:i*comm_rep_num);
%     save(['comm_all/comm' num2str(i*comm_rep_num)],'comm');
% end
