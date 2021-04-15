function [fp_manu,L_manu,N_manu,n_genos_curr] =...
    mutation(fp_manu,L_manu,N_manu,n_genos_curr,pot_mut_index,...
    p_mut, frac_null, sp0, sn0, g)
% generate vector of mutant cell numbers
if sum(N_manu(pot_mut_index)) > 0
    N_mut = fastbinorv(N_manu(pot_mut_index),p_mut);
    if sum(N_mut) > 0
        L_mut = L_manu(pot_mut_index);
        fp_mut = fp_manu(pot_mut_index);
        % remove mutants from their original genotypes
        N_manu(pot_mut_index) = N_manu(pot_mut_index) - N_mut;
        % remove elements of N_mut, L_mut, and p_mut with 0 cells.
        [fp_mut,L_mut,N_mut] = removeZeros(fp_mut,L_mut, N_mut);
        % generate vector of null mutant cell numbers
        N_null = fastbinorv(N_mut,frac_null);
        N_am = N_mut - N_null;
        % make this more space efficient.
        [fp_null,L_null,N_null] = removeZeros(fp_mut,L_mut,N_null);
        % update fp_manu, L_manu and N_manu, and n_genos_curr
        fp_manu(n_genos_curr + 1 :n_genos_curr + length(N_null)) = fp_null * 0;
        L_manu(n_genos_curr + 1 :n_genos_curr + length(N_null)) = L_null;
        N_manu(n_genos_curr + 1 :n_genos_curr + length(N_null)) = N_null;
        n_genos_curr = n_genos_curr + length(N_null);
        % generate vector of active mutant (am) cell numbers
        [fp_am,L_am,N_am] = removeZeros(fp_mut,L_mut,N_am);
        fp_back = zeros(sum(N_am),1);
        L_back = zeros(sum(N_am),1);
        % reformat vectors to execute non-null mutation.
        x = 1;
        for i = 1 : length(N_am)
            fp_back(x:x+N_am(i)) = fp_am(i);
            L_back(x:x+N_am(i)) = L_am(i);
            x = x + N_am(i) + 1;
        end
        if any(fp_back == 0)
            error('something went wrong here')
        end
        N_back = ones(length(fp_back),1);
        % update fp_manu, L_manu, N_manu, and n_genos_curr
        params = [sp0, sn0, g];
        fp_manu(n_genos_curr + 1 : n_genos_curr + length(fp_back)) = mutrnd_Dunham(params,fp_back);
        L_manu(n_genos_curr + 1 : n_genos_curr + length(fp_back)) = L_back;
        N_manu(n_genos_curr + 1 : n_genos_curr + length(fp_back)) = N_back;
        n_genos_curr = n_genos_curr + length(fp_back);
    end
end
end