function x = dec2logical(d,n)
% function x = dec2logical(d,n)
 d = d(:); % Make sure d is a column vector.

[f,e]=log2(max(d)); % How many digits do we need to represent the numbers?
x=rem(floor(d*pow2(1-max(n,e):0)),2);
end