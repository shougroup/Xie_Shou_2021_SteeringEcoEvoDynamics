function check_flag = decision_func(P_sel_dynamics, n, C, check_cycle_prev)
load('CheckDecPara')
check_flag = false;
% Heritability will be checked if the following conditions are met:
% * at least min_check_int cycles have passed since last check or the beginning of simulation
% * at least 2 cycles remains before the end of the simulation
% * over the last min_check_int cycles, the improvement in function is less than inc_thre/cycle
if n - check_cycle_prev > min_check_int && C-n > 2 && length(P_sel_dynamics) > min_check_int
    P_dynamics = P_sel_dynamics(n - min_check_int : n-1);
    % calculate the improvement rate with linear regression
    imp_rate = [ones(min_check_int, 1) (1:min_check_int)']\P_dynamics;
    if imp_rate(2) < inc_thre
        check_flag = true;
    end
end
end
    
