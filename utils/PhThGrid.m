function [ph,th] = PhThGrid(stepDeg,phLimDeg,thLimDeg)

% function [ph,th] = PhThGrid(stepDeg,phLimDeg,thLimDeg)
% returns the FarField class standard ph th grid vectors in radians
% Inputs all in degrees.
% stepDeg - can be 1 or 2 element array [stepPh, stepTh] (default [5,5])
% phLimDeg - limits of ph axis. 2 element array (default [0,360])
% thLimDeg - limits of th axis. 2 element array (default [0,360])


if nargin < 1
    stepDeg = [5,5];
else
    if numel(stepDeg) == 1
        stepDeg = ones(1,2).*stepDeg;
    elseif numel(stepDeg) > 2
        error('Expect stepDeg to have 1 or 2 elements');
    end
end

if nargin < 2
    phLimDeg = [0,360];
else
    if numel(phLimDeg) ~= 2
        error('Expect phLimDeg to have 1 or 2 elements');
    end
end

if nargin < 3
    thLimDeg = [0,180];
else
    if numel(phLimDeg) ~= 2
        error('Expect thLimDeg to have 1 or 2 elements');
    end
end

phVect = deg2rad(phLimDeg(1):stepDeg(1):phLimDeg(2));
thVect = deg2rad(thLimDeg(1):stepDeg(2):thLimDeg(2));

[PH,TH] = meshgrid(phVect,thVect);
ph = PH(:);
th = TH(:);
