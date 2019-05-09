function [ep,al] = DirCos2ElAz(u,v,w)
% function [ep,al] = DirCos2ElAz(u,v,w)
al = real(asin(u));
ep = real(atan2(v,w));
end