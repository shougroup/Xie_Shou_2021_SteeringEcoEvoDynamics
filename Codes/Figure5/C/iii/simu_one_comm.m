function comm_rep = simu_one_comm(newborn, comm_struct, const_struct)
rnodeseed = newborn.rseed;
rng(rnodeseed, 'twister');
t_binnum = const_struct.t_binnum;
max_popul = const_struct.max_popul;
pcs = const_struct.pcs;
R0 = const_struct.R0;
t_bin = const_struct.t_bin;
fp_Bound = const_struct.fp_Bound;
c_BM = const_struct.c_BM;
c_RM = const_struct.c_RM;
c_RH = const_struct.c_RH;
mut_rate = const_struct.mut_rate;
dM = const_struct.dM;
dH  = const_struct.dH;
gH_max = const_struct.gH_max ;
gM_max = const_struct.gM_max;
K_MB = const_struct.K_MB ;
K_MR = const_struct.K_MR ;
K_HR = const_struct.K_HR ;


M_counter = nnz(newborn.M_L>pcs);
H_counter = nnz(newborn.H_L>pcs);

M_L = zeros(max_popul, 1);
H_L = zeros(max_popul, 1);
fp = zeros(max_popul, 1);

% copy the data from the structure with Newborn's configuration
M_L(1:M_counter) = newborn.M_L(1:M_counter);
H_L(1:H_counter) = newborn.H_L(1:H_counter);
fp(1:M_counter) = newborn.fp(1:M_counter);
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
%         p0=sum(fp.*M_L./(1-fp));

M_t(1) = sum(M_L);
H_t(1) = sum(H_L);
R(1) = R0;

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
    fp(fp>fp_Bound) = fp_Bound;
    
    % if the biomass of a H cell is greater than 2, it divides
    div_idx = find(H_L >= 2);
    div_length = length(div_idx);
    if div_length>0
        H_L(H_counter+1:H_counter+div_length) = H_L(div_idx)/2;
        H_L(div_idx) = H_L(div_idx)/2;
        H_counter = H_counter+div_length;
    end
end
comm_rep = comm_struct;
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
comm_rep.parentnum = newborn.parentnum;
comm_rep.rseed = newborn.rseed;