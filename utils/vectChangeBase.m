function Vout = vectChangeBase(Vin,coor_new,coor_base)
% Transforms vector Vin, defined in coor_base, to coor_new

assert(size(Vin,1) == 3,'Vin must have 3 rows')

if nargin < 3 || isempty(coor_base)
    coor_base = CoordinateSystem();
end

p0base = Pnt3D;
p1base = p0base.addVect(Vin);

p0new = p0base.changeBase(coor_new,coor_base);
p1new = p1base.changeBase(coor_new,coor_base);

Vout = pointMatrix(p1new - p0new);
