function dy=chem_conc(t, y, para)

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

g_M = BN/(BN+RN)*RN/(RN+1)+RN/(BN+RN)*BN/(BN+1);
g_H = RN/(RN+KHKM);

dB = -gM_max*c_BM*g_M*M+gH_max*g_H*H;
dR = -gM_max*c_RM*g_M*M-gH_max*c_RH*g_H*H;

dy=[dB; dR];