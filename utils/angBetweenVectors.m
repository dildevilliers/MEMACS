function A = angBetweenVectors(U,V)
% ANGBETWEENVECTORS calculates the angle between the input vectors
% function A = angBetweenVectors(U,V)
% Returns the angle (in rad) between the 3D vectors U and V ([x;y;z])
% U and V must have 3 rows, and will be expanded in columns if one only has
% a single column

assert(size(U,1) == 3, 'U must have 3 rows')
assert(size(V,1) == 3, 'V must have 3 rows')

nU = size(U,2);
nV = size(V,2);

if nU ~= nV
    if nU == 1
        U = repmat(U,1,nV);
    elseif nV == 1
        V = repmat(V,1,nU);
    end
end

cV = cross(U,V);
cVn = (sum(abs(cV).^2,1)).^(1/2);
dV = dot(U,V);
A = atan2(cVn,dV);

