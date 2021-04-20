function [check_flag, rep_num_M] = check_decision_2(comm_all, I, spike_all, N, comm_type_num, BM_target)
sl = length(spike_all) - 1;
spike_frac = spike_all(1);
spike_test = spike_all(2:end);
comm_rep_num = N/comm_type_num;
check_flag = true;
% rep_num_M is a matrix for the numbers of newborns for selection
% and tests for heritability under different spiking fraction
rep_num_M = zeros(comm_type_num, sl+1); % N = 100, sl+1 = 3, # of substitution conditions;
for i = 1 : comm_type_num
    % BM is the biomass of the Adult
    BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
    % test_rep_temp is the total number of test comms that the Adult can
    % contribute alternative (not current) substitution conditions for
    % heritability checks
    test_rep_temp = floor((BM - comm_rep_num*BM_target*(1-abs(spike_frac)))...
        /BM_target/sum(1-abs(spike_test)));
    % if the number of test comms is too small, postpone the test cycle by 1
    if test_rep_temp < round(comm_rep_num/5)
        check_flag = false;
        break
    else
        rep_num_M(i, 1) = comm_rep_num;
        rep_num_M(i, 2:end) = min(test_rep_temp, comm_rep_num);
    end
end


