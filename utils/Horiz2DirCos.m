function [u,v,w] = Horiz2DirCos(az,alt)
% function [u,v,w] = Horiz2DirCos(az,alt)
% Try to maintain the ph angles in the pole

[ph,th] = Horiz2PhTh(az,alt);
[u,v,w] = PhTh2DirCos(ph,th);
end