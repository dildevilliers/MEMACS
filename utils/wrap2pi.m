function th = wrap2pi(th)
%wrap2pi Wrap angle in radians to [-pi pi]
%
%   th = wrap2pi(th) wraps angles in th, in radians,
%   to the interval [-pi pi] such that pi maps to pi and -pi maps to
%   -pi.  (In general, odd, positive multiples of pi map to pi and odd,
%   negative multiples of pi map to -pi.)

q = (th < -pi) | (pi < th);
th(q) = wrap22pi(th(q) + pi) - pi;
