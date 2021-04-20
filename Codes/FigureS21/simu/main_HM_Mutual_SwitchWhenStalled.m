clear
% dbstop if error
% number of previous cycles
C_prev = 0;
C = 3000; % total number of cycles
if C_prev == 0
    SPMutualParaInit
end
load('SPMutualPara')
%%
[spike_all, check_cycle_prev, H_isolates_in, M_isolates_in, ...
    newborns, P_sel_dynamics, check_cycle_m]...
    = SPMutualSimuInit(C_prev);
rng('shuffle'); % randomize the random number seed
n = C_prev+1;
check_flag = false;
% spike_frac = spike_all(1);
% spike_test = spike_all(2 : end);
% flag_m = zeros(C,2);
while n <= C
    %     flag_m(n, :) = [n uint8(check_flag)];
    % create a folder Cn to save the results of the nth cycle
    folder_name1 = ['C' num2str(n)];
    if ~exist(folder_name1, 'dir')
        mkdir(folder_name1)
    end
    folder_name2 = ['C' num2str(n) '/comm_all'];
    if ~exist(folder_name2, 'dir')
        mkdir(folder_name2)
    end
    comm_all(1,1:N) = comm_struct;
    % rep is the index of communities within one cycle
    parfor rep = 1:N
        if sum(newborns(rep).M_L) + sum(newborns(rep).H_L) < pcs
            error('An Newborn is empty')
        end
        comm_all(rep) = simu_one_Mutual(newborns(rep), comm_struct, const_struct);
    end
    newborns(1, 1:N) = newborn_struct;
    distrng=rng;
    if Pn_sig>pcs
        % Pn is the measurement noise in P(T)
        Pn = normrnd(0, Pn_sig, size([comm_all.P]));
        save([folder_name2 '/Pn'],'Pn');
    else
        Pn = 0;
    end
    P_all = [comm_all.P];
    save([folder_name2 '/P_all'], 'P_all');
    [~, I] = sort([comm_all.P]+Pn, 'descend');
    if max(P_all) < pcs
        break
    end
    if check_flag == false && decision_func(P_sel_dynamics, n, C, check_cycle_prev) == true
        [check_flag, rep_num_M]...
            = check_decision_2(comm_all, I, spike_all, N, comm_type_num, BM_target);
        if check_flag == false
            check_cycle_prev = n; 
            check_cycle = -1;
        else
            check_cycle = n+2;
        end
    end
    % in cycle check_cycle-2, parent communities are generated if there are enough cells
    if check_flag == true && abs(n - check_cycle + 2) < 0.1
        sel_counter = size(rep_num_M, 1);
        test_rep_total2 = sum(rep_num_M(:, end)); % total number of parent comms in check_cycle-2
        comm_selected = comm_all(I(1:sel_counter));
        % newborns_par is the structure array for parent Newborns
        newborns_par(1:sl, 1:test_rep_total2) = newborn_struct;
        % H isolates that might be used to spike during the next cycle
        H_isolates_out(1, 1:sel_counter) = H_isolates_struct;
        % M isolates that might be used to spike during the next cycle
        M_isolates_out(1, 1:sel_counter) = M_isolates_struct;
        rseed = randi(2^32-1, sl, test_rep_total2, 'uint32');
        comm_par_counter = 0;
        rep_counter = 0;
        for i = 1:sel_counter
            % reproduce the ith Adult into rep_num_M(i, 1) Newborns for selection,
            % rep_num_M(i, 2) Newborns for testing heritability when spiking
            % spike_test(1) H biomass, rep_num_M(i, 3) Newborns for testing
            % heritability when spiking spike_test(2) H biomass, and so on
            % comm_temp is a row structure array with sum(rep_num_M(i, :)) elements
            [comm_temp, H_isolates_out(i), M_isolates_out(i)] = pipette_SpikeMix_Mutual(comm_selected(i), newborn_struct, ...
                const_struct, spike_all, rep_num_M(i, :), ...
                H_isolates_in(comm_selected(i).parentnum), M_isolates_in(comm_selected(i).parentnum),...
                spike_clone_num, i);
            newborns(rep_counter+1:rep_counter+rep_num_M(i, 1)) ...
                = comm_temp(1:rep_num_M(i, 1));
            rep_counter = rep_counter + rep_num_M(i, 1);
            % newborns_par(j, :) contains parent Newborns with jth substitution fraction
            for j = 1:sl
                newborns_par(j, comm_par_counter+1:comm_par_counter+rep_num_M(i,end)) ...
                    = comm_temp(comm_rep_num+1+(j-1)*rep_num_M(i,end) : comm_rep_num+j*rep_num_M(i,end));
            end
            comm_par_counter = comm_par_counter+rep_num_M(i,end);
        end
        P_sel_dynamics(n, 1) = mean([comm_selected.P]);
        % If an Adult doesn't have H to be isolated for spiking for the next cycle, fill the spot with H isolates 
        % from the Adult with the highest P(T). If none of have H to be isolated for spiking, use those H isolates 
        % from the previous cycle.
        H_isolates_in = H_isolates_out;
        M_isolates_in = M_isolates_out;
        % assign random number seeds to each parent Newborn.
        for i = 1:sl
            for j = 1:test_rep_total2
                newborns_par(i,j).rseed = rseed(i,j);
            end
        end
        % in cycle check_cycle-1, offspring communities are generated if there are enough cells
    elseif check_flag == true && abs(n - check_cycle + 1) < 0.1
        [check_flag, rep_num_M, heri_par_idx] ...
            = check_decision_1(comm_all, I, spike_all, rep_num_M, N, comm_type_num, BM_target, off_rep_max, test_rep_max);
        if check_flag == true
            heri_par_counter = length(heri_par_idx);
            test_rep_total1 = min(test_rep_total2, test_rep_max); % number of ordered parent comms decided in cycle check_cycle-1
            sel_counter = size(rep_num_M, 1);
            check_cycle_m = [check_cycle_m; check_cycle];
            save('comm_all/check_cycle_m', 'check_cycle_m')
            H_isolates_out(1, 1:sel_counter) = H_isolates_struct;
            M_isolates_out(1, 1:sel_counter) = M_isolates_struct;
            % newborns_off is the structure array for offspring Newborns
            newborns_off(1:off_rep_max*(sl+1), 1:test_rep_total1) = newborn_struct;
            comm_selected = comm_all(I(1:sel_counter));
            % P_par is the array for P(T) of parents
            P_par = zeros(size(newborns_par));
            P_par_sorted = zeros(size(newborns_par));
            M0_par = zeros(size(newborns_par));
            H0_par = zeros(size(newborns_par));
            % adults_par is the structure array for parent Adults
            adults_par(1:sl, 1:test_rep_total2) = comm_struct;
            parfor i = 1:test_rep_total2*sl
                adults_par(i) = simu_one_Mutual(newborns_par(i), comm_struct, const_struct);
                P_par(i) = adults_par(i).P;
                M0_par(i) = adults_par(i).M_t(1);
                H0_par(i) = adults_par(i).H_t(1);
            end
            save([folder_name1 '/ParResults'],'P_par','M0_par','H0_par','spike_all')
            if Pn_sig>pcs
                Pn_par = normrnd(0, Pn_sig, size(P_par));
                save([folder_name1 '/ParResults'], 'Pn_par', '-append')
                P_par = P_par + Pn_par;
            end
            % sort the parent Adults according to their functions
            for i = 1:sl
                [P_par_sorted(i,:), I_par] = sort(P_par(i,:), 'descend');
                adults_par(i,:) = adults_par(i, I_par);
            end
            % only the top test_rep parent Adults reproduce
            adults_par = adults_par(:,1:test_rep_total1);
            rseed = randi(2^32-1,off_rep_max*(sl+1), test_rep_total1, 'uint32');
            % reproduce Adults of the current substitution ratio into Newborns for
            % selection and offspring Newborns for checking heritability
            rep_counter = 0;
            for i = 1 : sel_counter
                [comm_temp, H_isolates_out(i), M_isolates_out(i)] ...
                    = pipette_SpikeMix_Mutual(comm_selected(i), newborn_struct, const_struct, spike_all(1), ...
                    rep_num_M(i), H_isolates_in(comm_selected(i).parentnum), ...
                    M_isolates_in(comm_selected(i).parentnum), spike_clone_num, i);
                rep_num_temp = rep_num_M(i)-off_rep_max;
                if rep_num_temp > 0
                    newborns(rep_counter+1:rep_counter+rep_num_temp) ...
                        = comm_temp(1:rep_num_temp);
                    newborns_off(1 : off_rep_max, i) ...
                        = transpose(comm_temp(rep_num_temp+1:end));
                else
                    newborns(rep_counter+1:rep_counter+rep_num_M(i)) = comm_temp;
                end
                rep_counter = rep_counter+rep_num_temp;
            end
            P_sel_dynamics(n, 1) = mean([comm_selected.P]);
            for i = sel_counter+1:test_rep_total1
                BM = comm_all(I(i)).M_t(end) + comm_all(I(i)).H_t(end);
                dil_factor = floor(BM/BM_target/(1-abs(spike_all(1))));
                off_rep_num = min(off_rep_max, dil_factor);
                if dil_factor > 0
                    [comm_temp, ~, ~] = pipette_SpikeMix_Mutual(comm_all(I(i)), newborn_struct, ...
                        const_struct, spike_all(1), off_rep_num, H_isolates_in(comm_all(I(i)).parentnum),...
                        M_isolates_in(comm_all(I(i)).parentnum), spike_clone_num, i);
                    newborns_off(1 : off_rep_num, i) = transpose(comm_temp);
                    heri_par_counter = heri_par_counter+1;
                    heri_par_idx(heri_par_counter) = i;
                end
            end
            % reproduce each parent Adult into up to test_rep_num offspring Newborns
            for i = 1:sl % i is the index for substitution fraction
                for j = 1:test_rep_total1
                    BM = adults_par(i,j).M_t(end) + adults_par(i,j).H_t(end);
                    dil_factor = floor(BM/BM_target/(1-abs(spike_all(i+1))));
                    off_rep_num = min(off_rep_max, dil_factor);
                    [comm_temp, ~, ~] = pipette_SpikeMix_Mutual(adults_par(i,j), newborn_struct,...
                        const_struct, spike_all(i+1), off_rep_num, H_isolates_in(adults_par(i,j).parentnum),...
                        M_isolates_in(adults_par(i,j).parentnum), spike_clone_num, j);
                    newborns_off(i*off_rep_max+1:i*off_rep_max+off_rep_num, j) ...
                        = transpose(comm_temp);
                end
            end
            H_isolates_in = H_isolates_out;
            M_isolates_in = M_isolates_out;
            for i = 1:off_rep_max * (sl+1) * test_rep_total1
                newborns_off(i).rseed = rseed(i);
            end
            clear newborns_par adults_par
        else
            check_cycle = -1;
            check_cycle_prev = n;
        end
    elseif check_flag == true && abs(n - check_cycle) < 0.1
        % P_off is the array for P(T) of each offspring Adult
        P_off = zeros(size(newborns_off));
        M0_off = zeros(size(newborns_off));
        H0_off = zeros(size(newborns_off));
        parfor i=1 : test_rep_total1*off_rep_max*(sl+1)
            if sum(newborns_off(i).M_L)+sum(newborns_off(i).H_L) < pcs
                P_off(i) = NaN;
                M0_off(i) = 0;
                H0_off(i) = 0;
            else
                comm_temp = simu_one_Mutual(newborns_off(i), comm_struct, const_struct);
                P_off(i) = comm_temp.P;
                if P_off(i) <0
                    error('negative P_off')
                end
                M0_off(i) = comm_temp.M_t(1);
                H0_off(i) = comm_temp.H_t(1);
            end
        end
        save([folder_name1 '/OffResults'],'P_off','M0_off','H0_off')
        clear newborns_off
        load(['C' num2str(check_cycle-1) '/comm_all/P_all'])
        if Pn_sig>pcs
            Pn_off = normrnd(0, Pn_sig, size(P_off));
            save([folder_name1 '/OffResults'], 'Pn_off', '-append')
            P_off = P_off + Pn_off;
            load(['C' num2str(check_cycle-1) '/comm_all/Pn'])
            P_sorted = sort(P_all+Pn, 'descend');
        else
            P_sorted = sort(P_all, 'descend');
        end
        check_cycle_prev = check_cycle;
        [spike_all, heri, lb, ub] = switch_decision(spike_all, P_sorted, P_par_sorted, P_off, heri_par_idx, test_rep_total1, ...
            off_rep_max, HeriSwitch, n_bstrap, q);
        check_flag = false;
        save([folder_name1 '/OffResults'],'heri','lb','ub','heri_par_idx','-append')
        save([folder_name1 '/spike_all'], 'spike_all')
        clear heri_par_idx P_par_sorted
    end
    % Reproducing chosen Adults when heritability is not checked
    if check_flag == false
        if exist('comm_selected', 'var')>0 || max([newborns.rseed])>0
            error('Communities already reproduced')
        end
        rep_num_M = rep_num_cal(comm_all, I, spike_all, N, comm_rep_num, BM_target);
        rep_counter = 0;
        sel_counter = size(rep_num_M, 1);        
        H_isolates_out(1, 1:sel_counter) = H_isolates_struct;
        M_isolates_out(1, 1:sel_counter) = M_isolates_struct;
        comm_selected = comm_all(I(1:sel_counter));
        for i = 1 : sel_counter
            [newborns(rep_counter+1:rep_counter+rep_num_M(i)), H_isolates_out(i), M_isolates_out(i)]...
                = pipette_SpikeMix_Mutual(comm_selected(i), newborn_struct, ...
                const_struct, spike_all(1), rep_num_M(i), H_isolates_in(comm_selected(i).parentnum),...
                M_isolates_in(comm_selected(i).parentnum), spike_clone_num, i);
            rep_counter = rep_counter + rep_num_M(i);
        end
        H_isolates_in = H_isolates_out(1:sel_counter);
        M_isolates_in = M_isolates_out(1:sel_counter);
        P_sel_dynamics(n, 1) = mean([comm_selected.P]);
    end
    % assign random number seed to Newborns of the next cycle
    rseed = randi(2^32-1, N, 1, 'uint32');
    for ri = 1 : N
        newborns(ri).rseed = rseed(ri);
    end
    clear H_isolates_out M_isolates_out
    % save the selected communities and the current state of the random number generator
    save([folder_name1 '/comm_selected'],'comm_selected');
    save([folder_name1 '/H_isolates_in'],'H_isolates_in');
    save([folder_name1 '/M_isolates_in'],'M_isolates_in');
    save([folder_name1 '/distrng'],'distrng');
    clear comm_selected
    n = n+1;
    % save Newborns for the next cycle
    folder_name3=['C' num2str(n) '/comm_all'];
    if ~exist(folder_name3, 'dir')
        mkdir(folder_name3)
    end
    save([folder_name3 '/newborns'],'newborns');
    save('comm_all/P_sel_dynamics', 'P_sel_dynamics')
end
% save('comm_all/newborns', 'newborns')
save('comm_all/check_cycle_m', 'check_cycle_m')
save('comm_all/spike_all', 'spike_all')

