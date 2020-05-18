function y = rmlinear(a,X)
% Computes a linear function with independent variables (no interactions)
% of the form y = a0 + a1*x1 + a2*x2 + ... + an*xn.
% Returns:
% y: function value [m x 1]
% Inputs:
% a: coefficients [1 x n+1]
% X: independent variables [m x n]

[m,n] = size(X);

if length(a) ~= (n+1), error('a should have one more column than X'); end

A = repmat(reshape(a,1,n+1),m,1);
y = A(:,1) + sum(A(:,2:end).*X,2);
