function [Xd,M] = rotGenVect(X,U,th)

% [Xd,M] = rotGenVect(X,U,th)
% Returns the rotated about the general vector U vector Xd = [x';y';z'] of 
% the input vector X = [x;y;z] by angle th;
% Vectors can be rows of equal length
% The affine transformation matrix is returned in M

% Make sure U is normalized
assert(numel(U) == 3,'Expected 3 element vector in U')
U = U./norm(U);
ux = U(1);
uy = U(2);
uz = U(3);

cth = cos(th);
sth = sin(th);
cthm1 = (1 - cth);

M = [cth + ux.^2.*cthm1, ux.*uy.*cthm1 - uz.*sth, ux.*uz.*cthm1 + uy.*sth;
     uy.*ux.*cthm1 + uz.*sth, cth + uy.^2.*cthm1, uy.*uz.*cthm1 - ux.*sth;
     uz.*ux.*cthm1 - uy.*sth, uz.*uy.*cthm1 + ux.*sth, cth + uz.^2.*cthm1];
Xd = M*X;




