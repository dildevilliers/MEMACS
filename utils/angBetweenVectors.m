function A = angBetweenVectors(U,V)

% function A = angBetweenVectors(U,V)
% Returns the angle (in rad) between the 3D vectors U and V ([x;y;z])

A = atan2(norm(cross(U,V)),dot(U,V));
