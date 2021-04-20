function rep_num_M = rep_num_cal(comm_all, I, spike_all, N, comm_rep_num, BM_target)
rep_num_M = zeros(N, 1);
rep_counter = 0;
sel_counter = 0;
spike_frac = spike_all(1);
for i = 1 : N
    if rep_counter >= N
        break
    end
    dil_factor = floor((comm_all(I(i)).M_t(end)+comm_all(I(i)).H_t(end))/BM_target/(1-abs(spike_frac)));
    if dil_factor == 0
        error('# %d of chosen Adult of cycle %d is nearly empty', i, n)
    end
    rep_num_M(i) = min([dil_factor, comm_rep_num, N-rep_counter]);
    sel_counter = sel_counter+1;
    rep_counter = rep_counter+rep_num_M(i);
end
rep_num_M = rep_num_M(1 : sel_counter);
