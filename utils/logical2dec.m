function x = logical2dec(s)
% function x = logical2dec(s)
 [m,n] = size(s);

 twos = pow2(n-1:-1:0);
 x = sum(s .* twos(ones(m,1),:),2);
end