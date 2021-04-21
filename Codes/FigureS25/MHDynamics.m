function dx = MHDynamics(t, x, para)
M = x(1); % biomass of M
H = x(2); % biomass of H
B = x(3); % level of Byproduct
R = x(4); % level of Resource
P = x(5); % level of Product

if M < 1
    M = 0; % threshold for a single M cell
end

if H < 1
    H = 0; % threshold for a single H cell
end

gMmax = para(1);
gHmax = para(2);
c_BM = para(3);
c_RM = para(4);
c_RH = para(5);
fp = para(6);
K_MR = para(7);
K_HR = para(8);
K_MB = para(9);
d1 = para(10);
d2 = para(11);

BN = B/K_MB;
RN = R/K_MR;
KHKM = K_HR/K_MR;
% Y1BN = c_MB*KA;
% Y1sN = c_MR*K_MR;
% Y2sN = c_HR*K_MR;
% K2K1 = K_HR/K_MR;
% r2N = r2/KA;
g_M = BN/(BN+RN)*RN/(RN+1)+RN/(BN+RN)*BN/(BN+1);
g_H = RN/(RN+KHKM);
dM = (1-fp)*gMmax*g_M*M-d1*M;
dH = gHmax*g_H*H-d2*H;
dB = -gMmax*c_BM*g_M*M+gHmax*g_H*H;
dR = -gMmax*c_RM*g_M*M-gHmax*c_RH*g_H*H;
dP = fp*gMmax*g_M*M;


dx = [dM; dH; dB; dR; dP];
