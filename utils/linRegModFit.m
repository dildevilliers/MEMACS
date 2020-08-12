function lm = linRegModFit(X,y,modelspec)

% function lm = linRegModFit(X,y,modelspec)
% Fits a linear regression model, of type modelspec, on the training set
% X,y.
% X is a [m x n] data matrix, with y a [m x 1] vector of responses
% modelspec can be (empty - default 'linear'):
%  'constant' -  Model contains only a constant (intercept) term.
%  'linear' - Model contains an intercept and linear terms for each predictor. (Default)
%  'interactions' - Model contains an intercept, linear terms, and all products of pairs of distinct predictors (no squared terms).
%  'purequadratic' - Model contains an intercept, linear terms, and squared terms.
%  'quadratic' - Model contains an intercept, linear terms, interactions, and squared terms.
% The ouput 'lm' is a linear regression model for use with linRegModCalc.m
% See help lsqcurvefit for more info on the different variables in the
% structure - note X is replaced with a in this case
%
% Dirk de Villiers
% Created: 2015-06-25
% Updated: 2015-07-02
%
% 2015-06-25: Basic shell and functionality
% 2015-07-02: Reformulate to use formal normal equations (see wikipedia page on linear least squares fitting)
% ToDo: 'interactions', 'quadratic'

if nargin < 3, modelspec = 'linear'; end

% Basic error checking
Xsize = size(X);
ysize = size(y);
if Xsize(1) ~= ysize(1), error('X and y should have the same number of rows'); end

m = Xsize(1);
n = Xsize(2);

switch modelspec
    case 'constant'
        Xmat = [ones(m,1)];
        a = (Xmat'*Xmat)\(Xmat'*y);
    case 'linear'
        if m < (n+1), warning('Underdetermined system - consider supplying more samples'); end
        Xmat = [ones(m,1),X];
        a = (Xmat'*Xmat)\(Xmat'*y);
    case 'interactions'
        
    case 'purequadratic'
        if m < (2*n+1), warning('Underdetermined system - consider supplying more samples'); end
        Xmat = [ones(m,1),X,X.^2];
        a = (Xmat'*Xmat)\(Xmat'*y);
    case 'quadratic'
        
    otherwise % Choose default here
        warning('Unknown modelspec: Linear model assumed')
        modelspec = 'linear';
        
        
end

lm.modelspec = modelspec;
lm.a = a;





