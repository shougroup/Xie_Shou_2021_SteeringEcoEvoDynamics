function save_newborns(newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
    gen, resultsfolder,win_inds)
nb.newb_fp_manu = newb_fp_manu;
nb.newb_L_manu = newb_L_manu;
nb.newb_N_manu = newb_N_manu;
nb.newb_L_help = newb_L_help;
nb.newb_N_help = newb_N_help;
nb.newb_n_genos = newb_n_genos;
nb.gen = gen + 1;
% nb.win_inds = win_inds;
save(strcat(resultsfolder,'nb',num2str(gen + 1)),'-struct','nb');
end