function check_flag = decision_func(P_sel_dynamics, n, C, check_cycle_prev)
load('CheckDecPara')
check_flag = false;
% check heritability every check_period cycles
if n - check_cycle_prev +2 >= check_period && C-n > 2 
        check_flag = true;
end

    
