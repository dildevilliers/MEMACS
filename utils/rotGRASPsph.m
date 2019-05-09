function [sphAngOut] = rotGRASPsph(sphAngIn,angGRASP)

% function [sphAngOut] = rotSphCoorsGRASP(sphAngIn,angGRASP)
% rotates the spherical angles, specified in sphAngIn [ph;th] (rad) (can have
% multiple columns), by the GRASP angles in angGRASP [th,ph,ps].

[Ncomp,Nang] = size(sphAngIn);
if Ncomp ~= 2
    if Nang == 2
        sphAngIn = sphAngIn.';
    else
        error('sphAngIn should have 2 rows')
    end
end
[x,y,z] = PhTh2DirCos(sphAngIn(1,:),sphAngIn(2,:));
X = [x;y;z];
Xp = rotGRASP(X,angGRASP);
[php,thp] = DirCos2PhTh(Xp(1,:),Xp(2,:),Xp(3,:));
sphAngOut = [php;thp];

