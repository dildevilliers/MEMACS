function [Xg,Yg] = DirCos2TrueView(u,v,w)
% function [Xg,Yg] = DirCos2TrueView(u,v,w)
Ph = atan2(v,u);
Th = acos(w);
Xg = real(Th.*cos(Ph));
Yg = real(Th.*sin(Ph));
end