function [spike_all , heri, lb, ub] = switch_decision(spike_all, P_sorted, P_par_sorted, P_off, heri_par_idx, test_rep, off_rep_max, HeriSwitch, n_bstrap, q)
sl = length(spike_all) - 1;
% heri stores the heritability under each substitution fraction
heri = zeros(size(spike_all));
% ub stores the upper confidence interval under each substitution fraction
ub = zeros(size(spike_all));
% lb stores the lower confidence interval under each substitution fraction
lb = zeros(size(spike_all));
% heri stores the p_value of heritability under each substitution fraction
%         p_val = zeros(sl+1, 1);
% heritability is defined as the Spearman correlation coefficient between P(T) of
% parent Adults and average P(T) among offspring Adults from each lineage
heri(1) = slope_func((P_sorted(heri_par_idx))',...
    (nanmedian(P_off(1:off_rep_max, heri_par_idx), 1))');
% the confidence interval is estimated from bootstraping
[lb(1), ub(1)] = bstrap_itvl((P_sorted(heri_par_idx))',...
    (nanmedian(P_off(1:off_rep_max, heri_par_idx), 1))', @slope_func, n_bstrap, q);
for i = 1 : length(spike_all) - 1
    heri(i+1) = slope_func(P_par_sorted(i, 1:test_rep)', ...
        (nanmedian(P_off(i*off_rep_max+1:(i+1)*off_rep_max, 1:test_rep), 1))');
    [lb(i+1), ub(i+1)] = bstrap_itvl(P_par_sorted(i, 1:test_rep)',...
        (nanmedian(P_off(i*off_rep_max+1:(i+1)*off_rep_max, 1:test_rep), 1))', @slope_func, n_bstrap, q);
end

if HeriSwitch == 1
    [heri_sorted, I_heri] = sort(heri(2:end),'descend');
    spike_idx = [];
    for i = 1:sl
        % if the heritability of alternative spiking ratio is larger than
        % the 95% confidence interval, switch.
        if heri_sorted(i) > ub(1)
            
            %             % if the lb of the heritability of alternative spiking ratio is larger than
            %             % the ub of the heritability of the current spiking ratio, switch.
            %             if lb(I_heri(i)+1) > ub(1)
            spike_idx = I_heri(i)+1;
            break
        end
    end
    if ~isempty(spike_idx)
        spike_frac = spike_all(spike_idx);
        spike_all(spike_idx) = [];
        spike_test = spike_all;
        spike_all = [spike_frac spike_test];
    end
elseif HeriSwitch == -1
    [heri_sorted, I_heri] = sort(heri(2:end),'ascend');
    spike_idx = [];
    for i = 1:sl
        if heri_sorted(i) < lb(1)
            spike_idx = I_heri(i)+1;
            break
        end
    end
    if ~isempty(spike_idx)
        spike_frac = spike_all(spike_idx);
        spike_all(spike_idx) = [];
        spike_test = spike_all;
        spike_all = [spike_frac spike_test];
    end
elseif HeriSwitch == 0
    spike_all = spike_all(randperm(sl+1));
else
    error('HeriSwitch value not valid')
end