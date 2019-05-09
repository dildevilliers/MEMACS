function [az,alt] = DirCos2AzAlt(u,v,w)
% function [az,alt] = DirCos2AzAlt(u,v,w)

az = - real(atan2(v,u));
alt = pi/2 - real(acos(w));


end