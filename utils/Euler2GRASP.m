function angGRASP = Euler2GRASP(angEuler)
% function angGRASP = Euler2GRASP(angEuler)
[alpha,beta,gamma] = unpackEuler(angEuler);
th = beta;
ph = alpha;
ps = alpha + gamma;
angGRASP = [th,ph,ps];
end