% reproduce Adult communities into Newborns through cell sorting, so that 
% the biomass in the Newborns are fixed BM(0)=BM_target and
% phi_M(0)=phi_M(T) of the parent Adults
% comm_select: Adults to be reproduced
% const_struct: a structure to pass constants
% dil_factor: dilution factor
% rep_counter: current number of Newborn communities.
% parentnum: the rank of the Adult community
function comm_rep = cell_sort(comm_selected, comm_struct, const_struct,...
    dil_factor, ~, rep_counter, parentnum)
% maximal number of offspring community from one Adult
comm_rep_num = const_struct.comm_rep_num;
% minimal number of Adults allowed to reproduce
comm_type_num = const_struct.comm_type_num;
BM_target = const_struct.BM_target;
% comm_struct=struct('M_L',zeros(max_popul,1),'H_L',zeros(max_popul,1),'fp',zeros(max_popul,1),...
%     'M_t',zeros(t_binnum,1),'H_t',zeros(t_binnum,1),'R',zeros(t_binnum,1),'B',zeros(t_binnum,1),...
%     'P',0,'parentnum',0,'rseed',uint32(0));
comm_rep(1 : comm_rep_num,1) = comm_struct;

M_counter = nnz(comm_selected.M_L);
H_counter = nnz(comm_selected.H_L);
M_L = comm_selected.M_L(1:M_counter);
H_L = comm_selected.H_L(1:H_counter);
fp = comm_selected.fp(1:M_counter);

% the fraction of M's biomass in the Adult
Mfrac=sum(M_L)/(sum(M_L)+sum(H_L));
% randomize the H and M cells
rand_idx1 = randperm(M_counter);
rand_idx2 = randperm(H_counter);
M_L_rand = M_L(rand_idx1);
fp_rand = fp(rand_idx1);
H_L_rand = H_L(rand_idx2);

% partition indice for M cells
par_idx1 = [0; zeros(dil_factor, 1)];
% partition M cells into groups with biomass closest to without exceeding
% BM_target*phi_M(T) 
for i = 1 : dil_factor
    M_L_remain = M_L_rand(par_idx1(i)+1 : M_counter);
    M_remain_accu = cumsum(M_L_remain);
    [~, idx_temp] = min( abs(M_remain_accu - BM_target*Mfrac) );
    if M_remain_accu(idx_temp) - BM_target*Mfrac > 0
        par_idx1(i+1) = idx_temp - 1 + par_idx1(i);
    else
        par_idx1(i+1) = idx_temp + par_idx1(i);
    end
end
% partition H cells into groups with biomass closest to without exceeding
% BM_target*(1-phi_M(T))
par_idx2=[0; zeros(dil_factor, 1)];
for i = 1 : dil_factor
    H_L_remain = H_L_rand(par_idx2(i)+1 : H_counter);
    H_remain_accu = cumsum(H_L_remain);
    [~, idx_temp] = min( abs(H_remain_accu-BM_target*(1-Mfrac)) );
    if H_remain_accu(idx_temp) - BM_target*(1-Mfrac) > 0
        par_idx2(i+1) = idx_temp - 1 + par_idx2(i);
    else
        par_idx2(i+1) = idx_temp + par_idx2(i);
    end
end

% if the Adult can generate more Newborns than comm_rep_num, keep only
% comm_rep_num
if dil_factor >= comm_rep_num
    for i = 1 : comm_rep_num
        M_num=par_idx1(i+1)-par_idx1(i);
        H_num=par_idx2(i+1)-par_idx2(i);
        if M_num>=1
            comm_rep(i).M_L(1:M_num) = M_L_rand(par_idx1(i)+1:par_idx1(i+1));
            comm_rep(i).fp(1:M_num) = fp_rand(par_idx1(i)+1:par_idx1(i+1));
        end
        if H_num>=1
            comm_rep(i).H_L(1:H_num)= H_L_rand(par_idx2(i)+1:par_idx2(i+1));
        end
        comm_rep(i).parentnum=parentnum;
    end
    % keep only comm_type_num*comm_rep_num Newborns
    if rep_counter+comm_rep_num > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1 : comm_type_num*comm_rep_num-rep_counter);
    end
% if the Adult generates less Newborns than comm_rep_num, keep them all
else
    for i = 1:dil_factor
        M_num = par_idx1(i+1)-par_idx1(i);
        H_num = par_idx2(i+1)-par_idx2(i);
        if M_num>=1
            comm_rep(i).M_L(1:M_num) = M_L_rand(par_idx1(i)+1:par_idx1(i+1));
            comm_rep(i).fp(1:M_num) = fp_rand(par_idx1(i)+1:par_idx1(i+1));
        end
        if H_num>=1
            comm_rep(i).H_L(1:H_num)= H_L_rand(par_idx2(i)+1:par_idx2(i+1));
        end
        comm_rep(i).parentnum = parentnum;
    end
    % keep only comm_type_num*comm_rep_num Newborns
    if rep_counter + dil_factor > comm_type_num*comm_rep_num
        comm_rep = comm_rep(1 : comm_type_num*comm_rep_num-rep_counter);
    else
        comm_rep = comm_rep(1 : dil_factor);
    end
end
    
    
