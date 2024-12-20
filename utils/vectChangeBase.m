function Vout = vectChangeBase(Vin,coor_new,coor_base)
% Transforms vector Vin, defined in coor_base, to coor_new

assert(size(Vin,1) == 3,'Vin must have 3 rows')

if nargin < 3 || isempty(coor_base)
    coor_base = CoordinateSystem();
end

Pbase = [0,0,0;Vin(:).'].';
Pb = Pnt3D(Pbase(1,:),Pbase(2,:),Pbase(3,:));
Pn = pointMatrix(Pb.changeBase(coor_new,coor_base));
Vout = Pn(:,2:end) - Pn(:,1);

