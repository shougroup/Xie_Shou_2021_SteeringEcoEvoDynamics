function savedata(wellPlate_fp_manu, wellPlate_L_manu, wellPlate_N_manu, wellPlate_L_help, wellPlate_N_help, wellPlate_Bio_M,...
    wellPlate_Bio_H, wellPlate_R, wellPlate_B, wellPlate_P, wellPlate_n_genos, gen, resultsfolder, cw_params)

% cw_params = [compressWinner, digits_fp, digits_L];

if cw_params(1)
    for i = 1 : length(wellPlate_fp_manu)
        [F_temp, L_temp, N_temp, ngenos_temp] = compress(...
            wellPlate_fp_manu{i}, wellPlate_L_manu{i}, wellPlate_N_manu{i}, cw_params(2:3));
        wellPlate_fp_manu{i} = F_temp;
        wellPlate_L_manu{i} = L_temp;
        wellPlate_N_manu{i} = N_temp;
        wellPlate_n_genos{i} = ngenos_temp;
    end
end

wellPlate.fp_manu = wellPlate_fp_manu;
wellPlate.L_manu = wellPlate_L_manu;
wellPlate.N_manu = wellPlate_N_manu;
wellPlate.L_help = wellPlate_L_help;
wellPlate.N_help = wellPlate_N_help;
wellPlate.Bio_M = wellPlate_Bio_M;
wellPlate.Bio_H = wellPlate_Bio_H;
wellPlate.R = wellPlate_R;
wellPlate.B = wellPlate_B;
wellPlate.P = wellPlate_P;
wellPlate.n_genos = wellPlate_n_genos;
wellPlate.gen = gen;
wellPlate.DT = datetime;
save(strcat(resultsfolder,'gen',num2str(gen)),'-struct','wellPlate');

end