function Sr = renormS(S,Z0,Zr)

% function Sr = renormS(S,Z0,Zr)
% Renormalize the standard S-matrix in 'S' from Z0 to Zr impedance level.
% Z0 and Zr may be scalars or vectors of lenght n, where S = [n x n x f].

% Author: Dirk de Villiers
% Date  : 2015/12/18

[n1,n2,nf] = size(S);

if n1 ~= n2, error('S should be square for along the 3rd dimension'); end

nZ0 = length(Z0);
nZr = length(Zr);

if nZ0 == 1,
    Z0 = repmat(Z0,1,n1);
elseif nZ0 ~= n1
    warning('First entry of Z0 used as impedance for all ports');
    Z0 = repmat(Z0(1),1,n1);
end
if nZr == 1,
    Zr = repmat(Zr,1,n1);
elseif nZr ~= n1
    warning('First entry of Zr used as impedance for all ports');
    Zr = repmat(Zr(1),1,n1);
end

rn = (Zr - Z0)./(Zr + Z0);
R = diag(rn);

An = sqrt(Zr./Z0).*(1./(Zr + Z0));
A = diag(An);

Sr = zeros(n1,n2,nf);
for ff = 1:nf
    Sf = S(:,:,ff);
    d = eye(n1)-R*Sf;
    m1 = A\(Sf-R);
    m2 = m1/d;
    Sr(:,:,ff) = m2*A;
end


