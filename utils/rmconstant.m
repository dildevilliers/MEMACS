function y = rmconstant(a,X)
% Computes a constant function 
% of the form y = a0 
% Returns:
% y: function value [m x 1]
% Inputs:
% a: coefficients [1]
% X: independent variables [m x n] - not used

if length(a) > 1, error('a should be scalar'); end

[m,n] = size(X);

y = ones(m,1).*a;
