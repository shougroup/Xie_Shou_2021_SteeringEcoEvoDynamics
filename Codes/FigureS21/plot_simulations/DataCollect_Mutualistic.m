function DataCollect_Mutualistic(Fnum,C)
dataname = ['PlotData/Data' num2str(Fnum)];
pcs = 1e-15;
%%
fp_commmean = zeros(C, 1);
K_MR_commmean = zeros(C, 1);
K_HR_commmean = zeros(C, 1);
K_MB_commmean = zeros(C, 1);
g_Mmax_commmean = zeros(C, 1);
g_Hmax_commmean = zeros(C, 1);
B0_commmean = zeros(C, 1);
M0frac_commmean = zeros(C, 1);
MTfrac_commmean = zeros(C, 1);
P_commmean = zeros(C, 1);

for n = 1:C
    load([num2str(Fnum) '/C' num2str(n) '/comm_all/newborns.mat']);
    load([num2str(Fnum) '/C' num2str(n) '/comm_selected.mat']);
    comm_type_num = length(comm_selected);
    rseed_NB = [newborns.rseed];
    pnum_NB = [newborns.parentnum];
    
    fp_mean = zeros(comm_type_num, 1);
    K_MR_mean = zeros(comm_type_num, 1);
    K_HR_mean = zeros(comm_type_num, 1);
    K_MB_mean = zeros(comm_type_num, 1);
    gMmax_mean = zeros(comm_type_num, 1);
    gHmax_mean = zeros(comm_type_num, 1);
    B0_mean=zeros(comm_type_num,1);
    M0 = zeros(comm_type_num, 1);
    H0 = zeros(comm_type_num, 1);
    MT = zeros(comm_type_num, 1);
    HT = zeros(comm_type_num, 1);
    
    for i=1:comm_type_num
        M0(i) = comm_selected(i).M_t(1);
        H0(i) = comm_selected(i).H_t(1);
        MT(i) = comm_selected(i).M_t(end);
        HT(i) = comm_selected(i).H_t(end);
        % find selected Newborns based on rseed and parentnum
        c_num = find((rseed_NB == comm_selected(i).rseed) & (pnum_NB == comm_selected(i).parentnum));
        if length(c_num) ~= 1
            error('error in finding selected Newborns in cycle %d', n)
        end
        if sum(newborns(c_num).M_L) < pcs
            fp_mean(i) = 0;
        else
            fp_mean(i)=sum(newborns(c_num).fp.*newborns(c_num).M_L)/sum(newborns(c_num).M_L);
            K_MR_mean(i) = sum(newborns(c_num).K_MR.*newborns(c_num).M_L)/sum(newborns(c_num).M_L);
            K_MB_mean(i) = sum(newborns(c_num).K_MB.*newborns(c_num).M_L)/sum(newborns(c_num).M_L);
            gMmax_mean(i) = sum(newborns(c_num).gM_max.*newborns(c_num).M_L)/sum(newborns(c_num).M_L);
            B0_mean(i) = sum(newborns(sel_idx(i)).B0.*newborns(sel_idx(i)).H_L)/sum(newborns(sel_idx(i)).H_L);
            K_HR_mean(i) = sum(newborns(c_num).K_HR.*newborns(c_num).H_L)/sum(newborns(c_num).H_L);
            gHmax_mean(i) = sum(newborns(c_num).gH_max.*newborns(c_num).H_L)/sum(newborns(c_num).H_L);
        end
    end
    
    fp_commmean(n) = mean(fp_mean);
    K_MR_commmean(n) = mean(K_MR_mean);
    K_HR_commmean(n) = mean(K_HR_mean);
    K_MB_commmean(n) = mean(K_MB_mean);
    g_Mmax_commmean(n) = mean(gMmax_mean);
    g_Hmax_commmean(n) = mean(gHmax_mean);
    B0_commmean(n) = mean(B0_mean);
    M0frac_commmean(n)=mean(M0./(M0+H0));
    MTfrac_commmean(n)=mean(MT./(MT+HT));
    P_commmean(n)=mean([comm_selected.P]);
end

save(dataname,'fp_commmean','P_commmean','K_MR_commmean','K_HR_commmean',...
    'K_MB_commmean','g_Mmax_commmean','g_Hmax_commmean','M0frac_commmean',...
    'MTfrac_commmean','B0_commmean')