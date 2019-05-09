function Xd = rotGRASP(X,angGRASP)

% function Xd = rotGRASP(X,angGRASP)
% Returns the GRASP angle (angGRASP) rotated vector Xd = [x';y';z'] of the input
% vector X = [x;y;z]
% Vectors can be rows of equal length

coor_base = CoordinateSystem();     % Work from global coordinate system
% Make a GRASP rotated coordinate system
coor_new = CoordinateSystem();
coor_new = rotGRASP(coor_new,angGRASP);
X = Pnt3D(X(1,:),X(2,:),X(3,:));
% Xd = changeBase(X,coor_new,coor_base);
Xd = changeBase(X,coor_base,coor_new);
Xd = Xd.pointMatrix;
