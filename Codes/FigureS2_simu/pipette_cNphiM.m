function [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos]...
    = pipette_const_N(...
    newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
    winner_fp_manu,winner_L_manu,winner_N_manu,...
    winner_L_help,winner_N_help, winner_n_genos,target_inds,nD,cw_params,len_fp_manu_init,N0)

%PIPETTE CONST N and phiM
%   11/27/17
%   Alex Yuan in the Shou Lab

if length(target_inds) > nD
    error('The length of target_inds must be less than nD')
end

% execute data compression if desired
if cw_params(1)
    [winner_fp_manu, winner_L_manu, winner_N_manu, winner_n_genos] =...
        compress(winner_fp_manu, winner_L_manu, winner_N_manu,cw_params(2:3));
end

phiM_T = (winner_L_manu' * winner_N_manu) / ...
    (winner_L_manu' * winner_N_manu + winner_L_help * winner_N_help);
bio_M_target = phiM_T * N0;
bio_H_target = (1 - phiM_T) * N0;

%%

probs_M = ones(length(winner_N_manu),nD) * 1/nD;
% do regular pipetting, using the multinomial distribution
% where helper cells are included in the array
cellmat_M = mnrnd(winner_N_manu,probs_M);
donation_M = zeros(length(winner_N_manu),1);
% take "taxes" or "donations" from all columns of cellmat that have excess
% cells.
for i = 1 : nD
    while winner_L_manu' * cellmat_M(:,i) > bio_M_target
        ne_inds = cellmat_M(:,i)>0;
        temp1 = cellmat_M(ne_inds,i);
        temp2 = donation_M(ne_inds);
        [temp1,temp2] = ...
            draw_one_cell(temp1,temp2);
        cellmat_M(ne_inds,i) = temp1;
        donation_M(ne_inds) = temp2;
    end
end
% give back cells to any wells that are lacking cells, until they achieve
% biomass over N0. Then take back the last cell you gave.

% TODO: deal with unlikely edge case where you don't have enough cells to
% give and then take. I think i took care of this edge case because there
% should always be enough biomass to top every newborn off.
for i = 1 : nD
    % give
    while winner_L_manu' * cellmat_M(:,i) < bio_M_target
        ne_inds = donation_M > 0;
        temp1 = cellmat_M(ne_inds,i);
        temp2 = donation_M(ne_inds);
        [temp2,temp1,indx] = draw_one_cell(temp2,temp1);
        cellmat_M(ne_inds,i) = temp1;
        donation_M(ne_inds) = temp2;
    end
    % take
    if winner_L_manu' * cellmat_M(:,i) > bio_M_target
        temp1(indx) = temp1(indx) - 1;
        temp2(indx) = temp2(indx) + 1;
        cellmat_M(ne_inds,i) = temp1;
        donation_M(ne_inds) = temp2;
    end
end
biomat_M = winner_L_manu' * cellmat_M;

%% distribute cells
for i = 1 : length(target_inds)
    [nb_fp_temp,nb_L_temp,nb_N_temp] = removeZeros(winner_fp_manu, winner_L_manu, cellmat_M(1:end,i));
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

%% TODO: do the same for helper cells.
N_help_remaining = winner_N_help;
num_helpers_per_well = floor(bio_H_target / winner_L_help);
for i = 1 : length(target_inds)
    newb_N_help{target_inds(i)} = num_helpers_per_well;
    N_help_remaining = N_help_remaining - num_helpers_per_well;
    if N_help_remaining < 0
        error('you are in helper cell debt!')
    end
    newb_L_help{target_inds(i)} = winner_L_help;
end
end

