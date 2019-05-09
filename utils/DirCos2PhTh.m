function [ph,th] = DirCos2PhTh(u,v,w)
% function [ph,th] = DirCos2PhTh(u,v,w)
ph = real(atan2(v,u));
th = real(acos(w));
end