function [u,v,w] = PhTh2DirCos(ph,th)
% function [u,v,w] = PhTh2DirCos(ph,th)
% Try to maintain the ph angles in the pole
th(th == 0 & ph ~= 0) = eps;
u = real(sin(th).*cos(ph));
v = real(sin(th).*sin(ph));
w = real(cos(th));
end