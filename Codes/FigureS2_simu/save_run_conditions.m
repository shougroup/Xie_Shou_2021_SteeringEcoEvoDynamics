function save_run_conditions(cycle_duration, end_cycles, multiplier, num_wells,...
    cw_params, gc, time, ic_fp_manu, ic_L_manu, ic_N_manu, ic_L_help,...
    ic_N_help, ic_n_genos,resultsfolder, seeds, scoremode)

runconditions.cycle_duration = cycle_duration;
runconditions.end_cycles = end_cycles;
runconditions.multiplier = multiplier;
runconditions.num_wells = num_wells;
runconditions.cw_params = cw_params;
runconditions.time = time;
runconditions.ic_fp_manu = ic_fp_manu;
runconditions.ic_L_manu = ic_L_manu;
runconditions.ic_N_manu = ic_N_manu;
runconditions.ic_L_help = ic_L_help;
runconditions.ic_N_help = ic_N_help;
runconditions.ic_n_genos = ic_n_genos;
runconditions.seeds = seeds;
%runconditions.scoremode = scoremode;
save(strcat(resultsfolder,'run_conditions'),'-struct','runconditions');
save(strcat(resultsfolder,'gc','-struct'),'gc');

end
