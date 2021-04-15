function [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos] = pipette(...
    newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
    winner_fp_manu,winner_L_manu,winner_N_manu,...
    winner_L_help,winner_N_help, winner_n_genos,target_inds,nD,cw_params,len_fp_manu_init,N0)

%PIPETTE Summary of this function goes here
%   Detailed explanation goes here

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

newb_NH_vec = mnrnd(winner_N_help,ones(1,nD) * 1/nD);
for i = 1 : length(target_inds)
    newb_N_help{target_inds(i)} = newb_NH_vec(i);
    newb_L_help{target_inds(i)} = winner_L_help;
end
end

