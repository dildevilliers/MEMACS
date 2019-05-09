function [az,el] = DirCos2AzEl(u,v,w)
% function [az,el] = DirCos2AzEl(u,v,w)
el = real(asin(v));
az = real(atan2(u,w));
end