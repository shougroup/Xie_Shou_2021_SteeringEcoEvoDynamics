function [fp_manu,L_manu,N_manu,L_help,N_help,Bio_M,Bio_H,R,B,P,n_genos,gc]...
    = simulateOneWell_NoCost(gc,ic_fp_manu,ic_L_manu,ic_N_manu,ic_L_help,ic_N_help,ic_n_genos)
%UNTITLED2 Summary of this function goes here
% gc is a structure of global constants and is not allowed to change.



%% define anonymous functions
b_H = @(R) gc.b_Hmax .* R ./ (R + gc.K_HR);
R_M = @(R) R ./ gc.K_MR;
B_M = @(B) B ./ gc.K_MB;
b_M = @(R,B) gc.b_Mmax .* (R_M(R) .* B_M(B)) ./ (R_M(R) + B_M(B)) .*...
    (1 ./(R_M(R) + 1) + 1 ./(B_M(B) + 1));
R_prime = @(Bio_H, Bio_M, R, B) -b_H(R) * gc.c_RH * Bio_H - b_M(R,B) * gc.c_RM * Bio_M;
B_prime = @(Bio_H, Bio_M, R, B) b_H(R) * Bio_H - b_M(R,B) * gc.c_BM * Bio_M;
%% set up data structures

R = [gc.R_init,zeros(1,gc.nsteps)];
B = zeros(1,gc.nsteps + 1);
fp_manu = ic_fp_manu;
L_manu = ic_L_manu;
N_manu = ic_N_manu;
L_help = ic_L_help;
N_help = ic_N_help;
n_genos_curr = ic_n_genos;
Bio_H = [N_help * L_help,zeros(1,gc.nsteps)];
Bio_M = [L_manu' * N_manu, zeros(1,gc.nsteps)];
P = zeros(1,gc.nsteps + 1);
%% simulate!
for tstep = 1 : gc.nsteps
    % keep track of L_manu and N_manu that were calculated in previous step
    L_manu_old = L_manu;
    N_manu_old = N_manu;
    % do a numerical integration over the timestep;
    % RBFG is composed of 4 elements:
    % * R = resource
    % * B = byproduct
    % * F = integral of b_H
    % * G = integral of b_M
    RBFG_prime = @(RBFG)...
        [R_prime(Bio_H(tstep), Bio_M(tstep), RBFG(1), RBFG(2));...
        B_prime(Bio_H(tstep), Bio_M(tstep), RBFG(1), RBFG(2));...
        b_H(RBFG(1));...
        b_M(RBFG(1),RBFG(2))];
    RBFG_now = [R(tstep);B(tstep);0;0];
    tspan = [0, gc.dt];
    [t,RBFG] = ode23s(@(t,RBFG) RBFG_prime(RBFG), tspan, RBFG_now);
    % update resource and byproduct
    R(tstep + 1) = RBFG(end,1);
    B(tstep + 1) = RBFG(end,2);
    % update N_help and L_help
    int_H = RBFG(end,3);
    L_help = L_help * exp(int_H);
    help_divides = L_help >=2;
    if help_divides
        N_help = N_help * 2;
        L_help = L_help / 2;
    end
    % update N_manu and L_manu, and also determine lgt2.
    int_M = RBFG(end,4);
    L_manu = L_manu .* exp(int_M);%L_manu = L_manu .* exp((1 - fp_manu) * int_M);
    lgt2 = L_manu > 2;
    pot_mut_index = and(L_manu > 2, fp_manu > 0);
    N_manu(lgt2) = N_manu(lgt2) * 2;
    L_manu(lgt2) = L_manu(lgt2) / 2;
    % update product
    Bio_M(tstep + 1) = L_manu' * N_manu;
%     P(tstep + 1) = P(tstep) + sum(gc.r_P * fp_manu ./ (1 - fp_manu) .* (L_manu .* N_manu - L_manu_old .* N_manu_old));
    P(tstep + 1) = P(tstep) + sum(gc.r_P * fp_manu .* (L_manu .* N_manu - L_manu_old .* N_manu_old));
    % update N_manu and N_help by implementing death
    N_manu = N_manu - fastbinorv(N_manu,gc.d_M * gc.dt);  %round(N_manu * gc.d_M * gc.dt);
    N_help = N_help - fastbinorv(N_help,gc.d_H * gc.dt);  %round(N_help * gc.d_H * gc.dt);
    % mutation...
    [fp_manu,L_manu,N_manu,n_genos_curr] = mutation(fp_manu,L_manu,N_manu,n_genos_curr,pot_mut_index,...
        gc.p_mut, gc.frac_null, gc.sp0, gc.sn0, gc.g);
    % update biomass
    Bio_M(tstep + 1) = L_manu' * N_manu;
    Bio_H(tstep + 1) = L_help * N_help;
end
n_genos = n_genos_curr;
end
