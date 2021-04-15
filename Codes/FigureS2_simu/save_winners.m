function save_winners(wellPlate_fp_manu, wellPlate_L_manu, wellPlate_N_manu, wellPlate_L_help, wellPlate_N_help, wellPlate_Bio_M,...
    wellPlate_Bio_H, wellPlate_R, wellPlate_B, wellPlate_P, wellPlate_n_genos, gen, resultsfolder, cw_params, win_inds)

% cw_params = [compressWinner, digits_fp, digits_L];

if cw_params(1)
    for j = 1 : length(win_inds)
        i = win_inds(j);
        [F_temp, L_temp, N_temp, ngenos_temp] = compress(...
            wellPlate_fp_manu{i}, wellPlate_L_manu{i}, wellPlate_N_manu{i}, cw_params(2:3));
        wellPlate_fp_manu{i} = F_temp;
        wellPlate_L_manu{i} = L_temp;
        wellPlate_N_manu{i} = N_temp;
        wellPlate_n_genos{i} = ngenos_temp;
    end
end

wellPlate.fp_manu = wellPlate_fp_manu(win_inds);
wellPlate.L_manu = wellPlate_L_manu(win_inds);
wellPlate.N_manu = wellPlate_N_manu(win_inds);
wellPlate.L_help = wellPlate_L_help(win_inds);
wellPlate.N_help = wellPlate_N_help(win_inds);
wellPlate.Bio_M = wellPlate_Bio_M(win_inds);
wellPlate.Bio_H = wellPlate_Bio_H(win_inds);
wellPlate.R = wellPlate_R(win_inds);
wellPlate.B = wellPlate_B(win_inds);
wellPlate.P = wellPlate_P(win_inds);
wellPlate.n_genos = wellPlate_n_genos(win_inds);
wellPlate.win_inds = win_inds;
wellPlate.gen = gen;
wellPlate.DT = datetime;
save(strcat(resultsfolder,'gen',num2str(gen)),'-struct','wellPlate');

end