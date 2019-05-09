function Xd = rotEuler(X,angEuler)

% function Xd = rotEuler(X,angEuler)
% Returns the Euler angle (angGRASP) rotated vector Xd = [x';y';z'] of the input
% vector X = [x;y;z]

coor_base = CoordinateSystem();     % Work from global coordinate system
% Make a GRASP rotated coordinate system
coor_new = CoordinateSystem();
coor_new = rotEuler(coor_new,angEuler);
X = Pnt3D(X(1,:),X(2,:),X(3,:));
% Xd = changeBase(X,coor_new,coor_base);
Xd = changeBase(X,coor_base,coor_new);
Xd = Xd.pointMatrix;
