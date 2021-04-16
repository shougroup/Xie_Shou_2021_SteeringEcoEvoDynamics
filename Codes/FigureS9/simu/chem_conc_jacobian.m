function J=chem_conc_jacobian(t,y,para)

B=y(1);
R=y(2);

gM_max = para(1);
gH_max = para(2);
c_BM = para(3);
c_RM = para(4);
c_RH = para(5);
K_MR = para(6);
K_HR = para(7);
K_MB = para(8);
M = para(9);
H = para(10);

BN = B/K_MB;
RN = R/K_MR;
KHKM = K_HR/K_MR;

% g_M = BN/(BN+RN)*RN/(RN+1)+RN/(BN+RN)*BN/(BN+1);
% g_H = RN/(RN+KHKM);
% dB = -gM_max*c_BM*g_M*M+gH_max*g_H*H;
% dR = -gM_max*c_RM*g_M*M-gH_max*c_RH*g_H*H;

dg_MdBN = RN^2 / (BN+RN)^2 / (RN+1) + (RN^2-RN*BN^2) / (BN+RN)^2 / (BN+1)^2;
dg_MdRN = BN^2 / (BN+RN)^2 / (BN+1) + (BN^2-BN*RN^2) / (BN+RN)^2 / (RN+1)^2;
dg_HdRN = KHKM / (RN+KHKM)^2;

dBdBN = -gM_max*c_BM*dg_MdBN*M;
dBdRN = -gM_max*c_BM*dg_MdRN*M+gH_max*dg_HdRN*H;
dRdBN = -gM_max*c_RM*dg_MdBN*M;
dRdRN = -gM_max*c_RM*dg_MdRN*M-gH_max*c_RH*dg_HdRN*H; 

J=[dBdBN/K_MB dBdRN/K_MR; dRdBN/K_MB dRdRN/K_MR];