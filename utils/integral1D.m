function I = integral1D(X,Y,rule)

% Integrates Y(X) over the full X range using 
% the optional argument rule = ['trap'] or 'simpson'
% X and Y can only be 1D vectors of the same size

if nargin < 3
    % If we have a midpoint - use simpson
    if mod(length(Y),2) == 0
        rule = 'trap';
    else
        rule = 'simp';
    end
end

assert(all(size(X) == size(Y)),'Input vectors must be the same size');
assert(min(size(Y)) == 1,'Inputs X and Y must be 1D vectors');

switch rule
    case 'trap'
        I = trapz(X,Y);
    case 'simp'
%         if mod(length(X),2) == 0
%             warning('Simpsons rule not great for even number of points...')
%         end
        I = simpsons(Y,min(X),max(X));
    otherwise
        error(['Unknown rule: ',rule])
end

