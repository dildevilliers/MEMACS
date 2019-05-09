function [Xd,M] = trans3D(X,DEL)

% function [Xd,M] = trans3D(X,DEL)
% Returns the translated vector Xd = [x';y';z'] of the input
% vector X = [x;y;z] translated by distance DEL = [dx;dy;dz]
% Vectors can be rows of equal length
% The affine transformation matrix is returned in M

M = eye(4);
M(1:3,4) = DEL;

Np = length(X(1,:));
Xd = zeros(size(X));
for pp = 1:Np
    Xp = reshape(X(:,pp),3,1);
    [Xf] = M*[Xp;1];
    Xd(:,pp) = Xf(1:3);
end