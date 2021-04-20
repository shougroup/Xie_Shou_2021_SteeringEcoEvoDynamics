function comm_rep = simu_one_Mutual(newborn, comm_struct, const_struct)
rnodeseed = newborn.rseed;
rng(rnodeseed, 'twister');
t_binnum = const_struct.t_binnum;
max_popul = const_struct.max_popul;
pcs = const_struct.pcs;
R0 = const_struct.R0;
t_bin = const_struct.t_bin;
gM_max_Bound = const_struct.gM_max_Bound;
gH_max_Bound = const_struct.gH_max_Bound;
fp_Bound = const_struct.fp_Bound;
K_MB_Bound = const_struct.K_MB_Bound;
K_MR_Bound = const_struct.K_MR_Bound;
K_HR_Bound = const_struct.K_HR_Bound;
B0_Bound = const_struct.B0_Bound;
K_singular = const_struct.K_singular;
c_BM = const_struct.c_BM;
c_RM = const_struct.c_RM;
c_RH = const_struct.c_RH;
mut_rate = const_struct.mut_rate;
dM = const_struct.dM;
dH  = const_struct.dH;

M_counter = nnz(newborn.M_L>pcs);
H_counter = nnz(newborn.H_L>pcs);

M_L = zeros(max_popul, 1);
H_L = zeros(max_popul, 1);
fp = zeros(max_popul, 1);
gM_max = zeros(max_popul, 1);
gH_max = zeros(max_popul, 1);
K_MB = zeros(max_popul, 1);
K_MR = zeros(max_popul, 1);
K_HR = zeros(max_popul, 1);
B0 = zeros(max_popul, 1);
% copy the data from the structure with Newborn's configuration
M_L(1:M_counter) = newborn.M_L(1:M_counter);
H_L(1:H_counter) = newborn.H_L(1:H_counter);
fp(1:M_counter) = newborn.fp(1:M_counter);
gM_max(1:M_counter) = newborn.gM_max(1:M_counter);
K_MB(1:M_counter) = newborn.K_MB(1:M_counter);
K_MR(1:M_counter) = newborn.K_MR(1:M_counter);
gH_max(1:H_counter) = newborn.gH_max(1:H_counter);
K_HR(1:H_counter) = newborn.K_HR(1:H_counter);
B0(1:H_counter) = newborn.B0(1:H_counter);
M_t = zeros(t_binnum+1, 1);
H_t = zeros(t_binnum+1, 1);
%  temporary column vectors for M and H's biomass
M_LTemp = zeros(max_popul,1);
H_LTemp = zeros(max_popul,1);
% column vectors to store B and R's concentrations at each time
% step
B = zeros(t_binnum+1, 1);
R = zeros(t_binnum+1, 1);
P = 0;
%         p0=sum(fp.*M_L./(1-fp));

M_t(1) = sum(M_L);
H_t(1) = sum(H_L);
R(1) = R0;
% perform simulation through t_binnum time steps
for dt = 2 : t_binnum+1
    gM_maxM_L = gM_max(1:M_counter) .* M_L(1:M_counter);
    gH_maxH_L = gH_max(1:H_counter) .* H_L(1:H_counter);
    paras=struct('gM_maxM_L',gM_maxM_L,'gH_maxH_L',gH_maxH_L,...
        'c_BM', c_BM, 'c_RM', c_RM, 'c_RH', c_RH,...
        'K_MB', K_MB(1:M_counter),'K_MR',K_MR(1:M_counter),...
        'K_HR',K_HR(1:H_counter), 'B0', B0(1:H_counter));
    fhandle=@(t,y) chem_Mutual_jacobian(t, y , paras);
    options=odeset('Jacobian',fhandle,'RelTol',1e-5);
    [tx,y]=ode15s(@(t,y) chem_Mutual(t, y, paras), (0:1e-2:t_bin), [B(dt-1); R(dt-1)], options);
    if ~isreal(y)
        error('imaginary value')
    end
    % matrice where the row number is the number of cells and colomn
    % number is the B(t) and R(t) within the time step
    BN_mat = (1./K_MB(1:M_counter)) * y(:,1)';
    RN_M_mat = (1./K_MR(1:M_counter)) * y(:,2)';
    M_coef = BN_mat./(RN_M_mat+BN_mat).*RN_M_mat./(RN_M_mat+1) ...
        + RN_M_mat./(RN_M_mat+BN_mat).*BN_mat./(BN_mat+1);
    gMdt = trapz(tx, M_coef, 2) .* gM_max(1:M_counter) .* (1 - fp(1:M_counter));
    clear M_coef BN_mat RN_M_mat
    RN_H_mat = (1./K_HR(1:H_counter)) * y(:,2)';
    H_coef = RN_H_mat ./(RN_H_mat+1) .* exp(-1./B0(1:H_counter) * y(:,1)');
    gHdt = trapz(tx, H_coef, 2) .* gH_max(1:H_counter);
    clear H_coef RN_H_mat
    
    M_LTemp(1:M_counter) = exp(gMdt) .* M_L((1:M_counter));
    H_LTemp(1:H_counter)=exp(gHdt) .* H_L(1:H_counter);
    P = sum(fp(1:M_counter) .* (M_LTemp(1:M_counter) - M_L(1:M_counter))...
        ./(1-fp(1:M_counter)) ) + P;
    M_L(1:M_counter) = M_LTemp(1:M_counter);
    H_L(1:H_counter) = H_LTemp(1:H_counter);
    B(dt) = y(end, 1);
    R(dt) = y(end, 2);
    
    death_probability = rand(max_popul,1);
    M_L(death_probability < dM * t_bin) = 0;
    fp( death_probability < dM * t_bin) = 0;
    gM_max( death_probability < dM * t_bin) = 0;
    M_t(dt) = sum(M_L);
    
    death_probability = rand(max_popul,1);
    H_L( death_probability < dH * t_bin) = 0;
    gH_max( death_probability < dH * t_bin) = 0;
    H_t(dt) = sum(H_L);
    
    div_idx = find(M_L >= 2);
    div_length = length(div_idx);
    if div_length > 0
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        fp(M_counter+1 : M_counter+div_length)...
            = fp(div_idx) .* mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        K_MB(M_counter+1 : M_counter+div_length)...
            = K_MB(div_idx) ./mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        K_MR(M_counter+1 : M_counter+div_length)...
            = K_MR(div_idx) ./mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        gM_max(M_counter+1 : M_counter+div_length)...
            =gM_max(div_idx) .* mut_multiplier;
        
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        fp(div_idx) = fp(div_idx) .* mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        K_MB(div_idx)=K_MB(div_idx) ./mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        K_MR(div_idx)=K_MR(div_idx) ./mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        gM_max(div_idx)=gM_max(div_idx) .* mut_multiplier;
        
        M_L(M_counter+1:M_counter+div_length) = M_L(div_idx)/2;
        M_L(div_idx) = M_L(div_idx)/2;
        M_counter = M_counter + div_length;
        
        % if the phenotypes exceed the bounds, cap them at the bounds
        fp(fp > fp_Bound) = fp_Bound;
        K_MB((K_MB < K_MB_Bound)&(K_MB > pcs)) = K_MB_Bound;
        K_MR((K_MR < K_MR_Bound)&(K_MR > pcs)) = K_MR_Bound;
        gM_max(gM_max > gM_max_Bound) = gM_max_Bound;
    end
    
    div_idx = find(H_L >= 2);
    div_length = length(div_idx);
    if div_length > 0
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        gH_max(H_counter+1 : H_counter+div_length)...
            =gH_max(div_idx) .* mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        K_HR(H_counter+1 : H_counter+div_length)...
            = K_HR(div_idx) ./mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        B0(H_counter+1 : H_counter+div_length)...
            = B0(div_idx) .* mut_multiplier;
        
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        gH_max(div_idx) = gH_max(div_idx) .* mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        K_HR(div_idx) = K_HR(div_idx) ./mut_multiplier;
        mut_multiplier = mu_spec(div_length) .* double(rand(div_length,1) <= mut_rate) + 1;
        mut_multiplier = max(mut_multiplier, pcs);
        B0(div_idx) = B0(div_idx) .* mut_multiplier;
        
        H_L(H_counter+1 : H_counter+div_length) = H_L(div_idx)/2;
        H_L(div_idx) = H_L(div_idx)/2;
        H_counter = H_counter + div_length;
        
        gH_max(gH_max > gH_max_Bound) = gH_max_Bound;
        K_HR((K_HR < K_HR_Bound)&(K_HR > pcs)) = K_HR_Bound;
        B0(B0 > B0_Bound) = B0_Bound;
    end
    
end
% find the indice of the cells that are viable (positive growth
% rates and K smaller than the singular value)
temp1 = find((gM_max > pcs) & (K_MB < K_singular) & (K_MR < K_singular));
temp2 = find((gH_max > pcs) & (K_HR < K_singular));
M_t(end) = sum(M_L(temp1));
H_t(end) = sum(H_L(temp2));

comm_rep = comm_struct;
M_counter = length(temp1);
comm_rep.M_L(1:M_counter) = M_L(temp1);
H_counter = length(temp2);
comm_rep.H_L(1:H_counter) = H_L(temp2);

comm_rep.fp(1:M_counter) = fp(temp1);
comm_rep.gM_max(1:M_counter) = gM_max(temp1);
comm_rep.K_MB(1:M_counter) = K_MB(temp1);
comm_rep.K_MR(1:M_counter) = K_MR(temp1);
comm_rep.gH_max(1:H_counter) = gH_max(temp2);
comm_rep.K_HR(1:H_counter) = K_HR(temp2);
comm_rep.B0(1:H_counter) = B0(temp2);

comm_rep.M_t = M_t;
comm_rep.H_t = H_t;
comm_rep.B = B;
comm_rep.R = R;
comm_rep.P = P;
comm_rep.parentnum = newborn.parentnum;
comm_rep.rseed = newborn.rseed;