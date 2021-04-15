function comm_newborn = mix_pipette(comm_selected, comm_struct, const_struct)
comm_rep_num = const_struct.comm_rep_num;
comm_type_num = const_struct.comm_type_num;
max_popul = const_struct.max_popul;
M_L = zeros(max_popul*comm_type_num,1);
H_L = zeros(max_popul*comm_type_num,1);
fp = zeros(max_popul*comm_type_num,1);
M_counter = 0;
H_counter =0;
pcs = const_struct.pcs;
BM_target = const_struct.BM_target;
for i = 1: comm_type_num
    temp_idx = comm_selected(i).M_L > pcs;
    temp_l = nnz(temp_idx);
    M_L(M_counter+1 : M_counter+temp_l) = comm_selected(i).M_L(temp_idx);
    fp(M_counter+1 : M_counter+temp_l) = comm_selected(i).fp(temp_idx);
    M_counter = M_counter + temp_l;
    temp_idx = comm_selected(i).H_L > pcs;
    temp_l = nnz(temp_idx);
    H_L(H_counter+1 : H_counter+temp_l) = comm_selected(i).H_L(temp_idx);
    H_counter = H_counter + temp_l;
end  
M_L = M_L(1:M_counter);
H_L = H_L(1:H_counter);
dil_factor = floor((sum(M_L) + sum(H_L))/BM_target); 
rand_temp1 = randi(dil_factor, M_counter, 1);
rand_temp2 = randi(dil_factor, H_counter, 1);
comm_newborn(1 : comm_rep_num*comm_type_num) = comm_struct;
for i = 1 : comm_rep_num*comm_type_num
    temp_idx = find(rand_temp1 == i);
    M_num = length(temp_idx);
    if ~isempty(temp_idx)
        comm_newborn(i).M_L(1:M_num) = M_L(temp_idx);
        comm_newborn(i).fp(1:M_num) = fp(temp_idx);
    end
    temp_idx = find(rand_temp2 == i);
    H_num = length(temp_idx);
    if ~isempty(temp_idx)
        comm_newborn(i).H_L(1:H_num) = H_L(temp_idx);
    end
end

