function [Xd,M] = scaleVect(X,LAM)

% function [Xd,M] = scaleVect(X,LAM)
% Returns the scaled vector Xd = [x';y';z'] of the input
% vector X = [x;y;z] scaled by distance LAM = [lx;ly;lz]
% Vectors can be rows of equal length
% The affine transformation matrix is returned in M

M = eye(4);
M(1,1) = LAM(1);
M(2,2) = LAM(2);
M(3,3) = LAM(3);

Np = length(X(1,:));
Xd = zeros(size(X));
for pp = 1:Np
    Xp = reshape(X(:,pp),3,1);
    [Xf] = M*[Xp;1];
    Xd(:,pp) = Xf(1:3);
end