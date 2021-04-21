function dy=chem_conc_SP(t,y,paras)

B = y(1);
R = y(2);

gM_maxM_L = paras.gM_maxM_L;
gH_maxH_L=paras.gH_maxH_L;
c_BM = paras.c_BM;
c_RM = paras.c_RM;
c_RH = paras.c_RH;
K_MB = paras.K_MB;
K_MR = paras.K_MR;
K_HR =paras.K_HR;

BN = B./K_MB;
RN_M = R./K_MR;
RN_H = R./K_HR;

% the term gM(R, B)/gM_max
M_coef = BN./(RN_M+BN).*RN_M./(RN_M+1) + RN_M./(RN_M+BN).*BN./(BN+1);
% the term gH(R)/gH_max
H_coef=RN_H./(RN_H+1);
dB = -sum(c_BM*gM_maxM_L.*M_coef) + sum(gH_maxH_L.*H_coef);
dR = -sum(c_RM*gM_maxM_L.*M_coef) - sum(c_RH*gH_maxH_L.*H_coef);

dy=[dB; dR];