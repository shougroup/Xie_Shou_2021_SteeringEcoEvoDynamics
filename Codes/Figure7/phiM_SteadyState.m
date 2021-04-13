function phiM = phiM_SteadyState(para)
phiM = 0.7;
dphiM = 1;
T0 = 17;
R0 = 1;
options=odeset('RelTol',1e-6,'abstol',1e-10);
while abs(dphiM) > 1e-3
[T, X] = ode15s(@(t,x) MHDynamics(t,x,para), [0 T0], [phiM*100;100*(1-phiM);0;R0;0],options);
dphiM = X(end, 1)/(X(end, 1)+X(end, 2))-phiM;
phiM = X(end, 1)/(X(end, 1)+X(end, 2));
end