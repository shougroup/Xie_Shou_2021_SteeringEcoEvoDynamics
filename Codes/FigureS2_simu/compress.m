function [fp_manu, L_manu, N_manu, n_genos] = compress(fp_manu, L_manu, N_manu, dFP_dL)

% cw_params = [compressWinner, digits_fp, digits_L];

% 'compresses' the three arrrays by binning fp and L.
nbins_fp = 10.^dFP_dL(1);
nbins_L = 10.^dFP_dL(2);
lossy_fp = round(nbins_fp * fp_manu) / nbins_fp;
lossy_L = round(nbins_L * L_manu) / nbins_L;
FLN = sortrows([lossy_fp,lossy_L,N_manu],[1,2]);
for i = 2 : length(FLN)
    if isequal(FLN(i,1:2),FLN(i-1,1:2))
        FLN(i,3) = FLN(i,3) + FLN(i-1,3);
        FLN(i-1,3) = 0;
    end
end
FLN(FLN(:,3) == 0,:) = [];
fp_manu = FLN(:,1);
L_manu = FLN(:,2);
N_manu = FLN(:,3);
n_genos = length(fp_manu);
end