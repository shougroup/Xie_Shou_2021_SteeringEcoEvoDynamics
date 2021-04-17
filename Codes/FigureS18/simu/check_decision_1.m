function [check_flag, rep_num_M, heri_par_idx] ...
    = check_decision_1(comm_all, I, spike_all, rep_num_M, N, comm_type_num, BM_target, off_rep_max, test_rep_max)
if size(rep_num_M, 2) <= 1
    error('previous cycle doesnt seem like check_cycle - 1')
end
spike_frac = spike_all(1);
comm_rep_num = N/comm_type_num;
check_flag = true;
test_rep_total2 = sum(rep_num_M(:, end)); % number of Adults under each alternative spiking strategy
test_rep_total1 = min(test_rep_total2, test_rep_max); % number of Adults used to estimate heritability
sel_counter = 0;
rep_counter = 0;
rep_num_M = zeros(N, 1);
heri_par_idx = [];
heri_par_counter = 0;
% Calculate how many Adults are required to check the heritability. If more than
% test_rep are needed, postpond heritability check by 10 cycles.
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
    % if a chosen Adult produces less than off_rep_max Newborns, skip this Adult
    % for heritability check
    if rep_num_temp > 0
        heri_par_counter = heri_par_counter+1;
        heri_par_idx(heri_par_counter) = i;
        rep_num_M(i) = rep_num_temp+off_rep_max;
    else
        rep_num_temp = min([dil_factor, comm_rep_num, N-rep_counter]);
        rep_num_M(i) = rep_num_temp;
    end
    sel_counter = sel_counter+1;
    rep_counter = rep_counter+rep_num_temp;
end
rep_num_M = rep_num_M(1:sel_counter);
% if more than test_rep Adults are required to check heritability, postpond the
% check by 10 cycles
if sel_counter > test_rep_total1
    check_flag = false;
    heri_par_idx = [];
    rep_num_M = [];
end

