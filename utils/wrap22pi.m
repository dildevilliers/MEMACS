function th = wrap22pi(th)
%wrap22pi Wrap angle in radians to [0 2*pi]
%
%   th = wrap22pi(th) wraps angles in LAMBDA, in radians,
%   to the interval [0 2*pi] such that zero maps to zero and 2*pi maps
%   to 2*pi. (In general, positive multiples of 2*pi map to 2*pi and
%   negative multiples of 2*pi map to zero.)


positiveInput = (th > 0);
th = mod(th, 2*pi);
th((th == 0) & positiveInput) = 2*pi;
