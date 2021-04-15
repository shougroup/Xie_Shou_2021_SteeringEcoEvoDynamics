function [newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos]...
    = pipette_const_N(...
    newb_fp_manu,newb_L_manu,newb_N_manu,newb_L_help,newb_N_help,newb_n_genos,...
    winner_fp_manu,winner_L_manu,winner_N_manu,...
    winner_L_help,winner_N_help, winner_n_genos,target_inds,nD,cw_params,len_fp_manu_init,N0)

%PIPETTE CONST N
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

probs = ones(length(winner_N_manu)+1,nD) * 1/nD;
% do regular pipetting, using the multinomial distribution
% where helper cells are included in the array
cellmat = mnrnd([winner_N_manu;winner_N_help],probs);
donation = zeros(length(winner_N_manu)+1,1);
% take "taxes" or "donations" from all columns of cellmat that have excess
% cells.
for i = 1 : nD
    while [winner_L_manu; winner_L_help]' * cellmat(:,i) > N0
        ne_inds = cellmat(:,i)>0;
        temp1 = cellmat(ne_inds,i);
        temp2 = donation(ne_inds);
        [temp1,temp2] = ...
            draw_one_cell(temp1,temp2);
        cellmat(ne_inds,i) = temp1;
        donation(ne_inds) = temp2;
    end
end
% give back cells to any wells that are lacking cells, until they achieve
% biomass over N0. Then take back the last cell you gave.

% TODO: deal with unlikely edge case where you don't have enough cells to
% give and then take. I think i took care of this edge case because there
% should always be enough biomass to top every newborn off.
for i = 1 : nD
    % give
    while [winner_L_manu; winner_L_help]' * cellmat(:,i) < N0
        ne_inds = donation > 0;
        temp1 = cellmat(ne_inds,i);
        temp2 = donation(ne_inds);
        [temp2,temp1,indx] = draw_one_cell(temp2,temp1);
        cellmat(ne_inds,i) = temp1;
        donation(ne_inds) = temp2;
    end
    % take
    if [winner_L_manu; winner_L_help]' * cellmat(:,i) > N0
        temp1(indx) = temp1(indx) - 1;
        temp2(indx) = temp2(indx) + 1;
        cellmat(ne_inds,i) = temp1;
        donation(ne_inds) = temp2;
    end
end
biomat = [winner_L_manu; winner_L_help]' * cellmat;

%% distribute cells
for i = 1 : length(target_inds)
    [nb_fp_temp,nb_L_temp,nb_N_temp] = removeZeros(winner_fp_manu, winner_L_manu, cellmat(1:end-1,i));
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

for i = 1 : length(target_inds)
    newb_N_help{target_inds(i)} = cellmat(end,i);
    newb_L_help{target_inds(i)} = winner_L_help;
end
end

