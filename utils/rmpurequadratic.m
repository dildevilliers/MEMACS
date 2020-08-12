function y = rmpurequadratic(a,X)
% Computes a quadratic function with independent variables (no interactions)
% of the form y = a0 + a1*x1 + a2*x2 + ... + an*xn + b1*x1^2 + b2*x2^2...
% with b_i = a_(n+i)
% Returns:
% y: function value [m x 1]
% Inputs:
% a: coefficients [1 x 2n+1]
% X: independent variables [m x n]

[m,n] = size(X);

if length(a) ~= (2*n+1), error('a should have 2*n+1 columns if X has n columns'); end

A = repmat(reshape(a,1,2*n+1),m,1);
y = A(:,1) + sum(A(:,2:n+1).*X,2) + sum(A(:,n+2:end).*X.^2,2);
