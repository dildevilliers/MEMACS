function [u,v,w] = DirCos2DirCos(u,v,w)
% function [u,v,w] = DirCos2DirCos(u,v,w)
if nargin < 3
    w = sqrt(1 - u.^2 - v.^2);
end
end