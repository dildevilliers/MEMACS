function [Xd,M] = roty(X,th)

% function [Xd,M] = roty(X,th)
% Returns the rotated about the y-axis vector Xd = [x';y';z'] of the input
% vector X = [x;y;z] by angle th;
% Vectors can be rows of equal length
% The affine transformation matrix is returned in M

M = [cos(th),0,sin(th),0;0,1,0,0;-sin(th),0,cos(th),0;0,0,0,1];

Np = length(X(1,:));
Xd = zeros(size(X));
for pp = 1:Np
    Xp = reshape(X(:,pp),3,1);
    [Xf] = M*[Xp;1];
    Xd(:,pp) = Xf(1:3);
end