function [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos] = pipette(...
    newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
    winner_fp_manu,winner_L_manu,winner_N_manu,...
    winner_L_help,winner_N_help, winner_n_genos,target_inds,nD,cw_params,len_fp_manu_init)

%PIPETTE Summary of this function goes here
%   Nov 26, 2017
%   Alex Yuan in the Shou Lab

if length(target_inds) > nD
    error('The length of target_inds must be less than nD')
end

% execute data compression if desired
if cw_params(1)
    [winner_fp_manu, winner_L_manu, winner_N_manu, winner_n_genos] =...
        compress(winner_fp_manu, winner_L_manu, winner_N_manu,cw_params(2:3));
end

% distribute manufacturer cells from winning well
probs = ones(length(winner_N_manu),nD) * 1/nD;
newb_N_mat = mnrnd(winner_N_manu,probs);
for i = 1 : length(target_inds)
    [nb_fp_temp,nb_L_temp,nb_N_temp] = removeZeros(winner_fp_manu, winner_L_manu, newb_N_mat(:,i));
    temp_n_genos = length(nb_fp_temp);
    if temp_n_genos < len_fp_manu_init
        nb_fp_temp = [nb_fp_temp;zeros(len_fp_manu_init-temp_n_genos,1)];
        nb_L_temp = [nb_L_temp;zeros(len_fp_manu_init-temp_n_genos,1)];
        nb_N_temp = [nb_N_temp;zeros(len_fp_manu_init-temp_n_genos,1)];
    end
    newb_fp_manu{target_inds(i)} = nb_fp_temp;
    newb_L_manu{target_inds(i)} = nb_L_temp;
    newb_N_manu{target_inds(i)} = nb_N_temp;
    newb_n_genos{target_inds(i)} = temp_n_genos;
end

% distribute helper cells from winning well
phiM_T = (winner_L_manu' * winner_N_manu) / ...
    (winner_L_manu' * winner_N_manu + winner_L_help * winner_N_help);

N_help_remaining = winner_N_help;
for i = 1 : length(target_inds)
    % calculate bioH as bioH = bioM / phiM - bioM
    bioM_temp = (newb_L_manu{target_inds(i)})' * (newb_N_manu{target_inds(i)});
    bioH_temp = bioM_temp / phiM_T - bioM_temp;
    % calcualate the number of H cells that you allocate to the ith newborn
    N_help_temp = floor(bioH_temp / winner_L_help);
    newb_N_help{target_inds(i)} = N_help_temp;
    newb_L_help{target_inds(i)} = winner_L_help;
    N_help_remaining = N_help_remaining - N_help_temp;
    % check to make sure you don't run out of helper cells.
    if N_help_remaining < 0
        error('you went into helper cell debt!')
    end
end
end

