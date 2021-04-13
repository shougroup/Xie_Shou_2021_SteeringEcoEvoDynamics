function J=chem_conc_jacobian_SP(t, y, paras)

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

% M_coef = BN./(RN_M+BN).*RN_M./(RN_M+1) + RN_M./(RN_M+BN).*BN./(BN+1);
% H_coef=RN_H./(RN_H+1);

dM_coefdB = RN_M.^2./(RN_M+BN).^2./(RN_M+1) ...
    + RN_M.*(1-BN.^2)./(RN_M+BN).^2./(BN+1).^2;
dM_coefdR = BN.^2./(RN_M+BN).^2./(BN+1) ...
    + BN.*(1-RN_M.^2)./(RN_M+BN).^2./(RN_M+1).^2;
dH_coefdR = 1./(RN_H+1).^2./K_HR;
dM_coefdB = dM_coefdB./K_MB;
dM_coefdR = dM_coefdR./K_MR;

% dB = -sum(c_BM*gM_maxM_L.*M_coef) + sum(gH_maxH_L.*H_coef);
% dR = -sum(c_RM*gM_maxM_L.*M_coef) - sum(c_RH*gH_maxH_L.*H_coef);
J11 = -sum(c_BM*gM_maxM_L.*dM_coefdB);
J12 = -sum(c_BM*gM_maxM_L.*dM_coefdR) + sum(gH_maxH_L.*dH_coefdR);
J21 = -sum(c_RM*gM_maxM_L.*dM_coefdB);
J22 = -sum(c_RM*gM_maxM_L.*dM_coefdR) - sum(c_RH*gH_maxH_L.*dH_coefdR);
J=[J11 J12; J21 J22];