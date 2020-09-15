% Generate data for parent-offspring comparison. For every parent community, reproduce
% into as many offspring communities as possible. 

clear

spike_frac = 0; % fraction of H pure culture spiked in
n = 549;
comm_type_num = 1;
mut_rate = 1e-2;

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
fp_max = 1; % fp is between 0 and 1

% N0: BM_0, R(0)=1, T0 is the maturation time
N0 = 100;
R0 = 1;
T0 = 17;

N = 1000; % to allow parent community to reproduce as many offspring as possible 
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

const_struct=struct('t_binnum',t_binnum,'max_popul',max_popul,'comm_rep_num',comm_rep_num,...
    'comm_type_num',comm_type_num,'pcs',pcs,'BM_target',BM_target);

filename = ['C' num2str(n) '/comm_all/adults.mat'];
load(filename);
% load('R1/adults.mat')
if ~exist('Data','dir')
    mkdir('Data')
end
% rng('shuffle')
% distrng_m = randi(2^32-1,comm_type_num*comm_rep_num,1,'uint32');
% save('Data/distrng_m', 'distrng_m')

% set the randome number generator to obtain identical results
load('distrng_m')
rng(distrng_m)
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
    for ri = 1 : dil_factor 
        comm_all(ri).rseed = rseed(ri);
    end
    newborns = comm_all;
    % save the selected communities and the current state of the random number generator 
    save(['Data/newborns' num2str(j)],'newborns');
    parfor rep = 1 : length(comm_all)
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
    adults_off = comm_all;
    save(['Data/adults' num2str(j)],'adults_off');
end

