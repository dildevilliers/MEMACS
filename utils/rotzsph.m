function [sphAngOut] = rotzsph(sphAngIn,th)

% function [sphAngOut] = rotzsph(sphAngIn,angGRASP)
% rotates the spherical angles, specified in sphAngIn [ph;th] (rad) (can have
% multiple columns), around the z-axis by the angle th in (rad).

[Ncomp,Nang] = size(sphAngIn);
if Ncomp ~= 2
    if Nang == 2
        sphAngIn = sphAngIn.';
    else
        error('sphAngIn should have 2 rows')
    end
end
[x,y,z] = PhTh2DirCos(sphAngIn(1,:),sphAngIn(2,:));
Xp = rotz([x;y;z],th);
[php,thp] = DirCos2PhTh(Xp(1,:),Xp(2,:),Xp(3,:));
sphAngOut = [php;thp];

