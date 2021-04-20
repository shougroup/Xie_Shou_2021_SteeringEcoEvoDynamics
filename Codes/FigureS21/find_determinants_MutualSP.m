function [para, BM0] = find_determinants_MutualSP(folder_name, cycle, cnum)
load([folder_name  '/C' num2str(cycle) '/comm_selected.mat']);
load([folder_name '/C' num2str(cycle) '/comm_all/newborns.mat']);
NB_rseed = [newborns.rseed];
NB_pnum = [newborns.parentnum];
c_BM = 1/3;
c_RM = 1e-4;
c_RH = 1e-4;
dM = 0.7*5e-3;
dH = 0.3*5e-3;
pcs = 1e-9;
temp_idx = find(NB_rseed == comm_selected(cnum).rseed & NB_pnum == comm_selected(cnum).parentnum);
if length(temp_idx) == 1
    sel_idx = temp_idx;
else
    error('cannot find NB of Adult %d in Cycle %d ', cnum, cycle)
end

M0 = comm_selected(cnum).M_t(1);
H0 = comm_selected(cnum).H_t(1);
if ~isempty(newborns(sel_idx).M_L>pcs)
    fp = sum(newborns(sel_idx).fp.*newborns(sel_idx).M_L)/M0;
    K_MR = sum(newborns(sel_idx).K_MR.*newborns(sel_idx).M_L)/M0;
    K_MB = sum(newborns(sel_idx).K_MB.*newborns(sel_idx).M_L)/M0;
    gMmax = sum(newborns(sel_idx).gM_max.*newborns(sel_idx).M_L)/M0;
    K_HR = sum(newborns(sel_idx).K_HR.*newborns(sel_idx).H_L)/H0;
    B0 = sum(newborns(sel_idx).B0.*newborns(sel_idx).H_L)/H0;
    gHmax = sum(newborns(sel_idx).gH_max.*newborns(sel_idx).H_L)/H0;
end

para = [gMmax; gHmax; c_BM; c_RM; c_RH; fp; K_MR; K_HR; K_MB; dM; dH; B0];
BM0 = [M0 H0];