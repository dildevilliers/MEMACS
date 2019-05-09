function [RA,dec] = DirCos2RAdec(u,v,w)
% function [RA,dec] = DirCos2RAdec(u,v,w)
RA = real(atan2(v,u));
dec = pi/2 - real(acos(w));
end