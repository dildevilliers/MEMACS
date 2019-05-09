function angEuler = GRASP2Euler(angGRASP)
% function angEuler = GRASP2Euler(angGRASP)
[th,ph,ps] = unpackGRASP(angGRASP);
alpha = ph;
beta = th;
gamma = -ph + ps;
angEuler = [alpha,beta,gamma];
end