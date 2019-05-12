function [az,alt] = DirCos2Horiz(u,v,w)
% function [az,alt] = DirCos2Horiz(u,v,w)

[ph,th] = DirCos2PhTh(u,v,w);
[az,alt] = PhTh2Horiz(ph,th);
end