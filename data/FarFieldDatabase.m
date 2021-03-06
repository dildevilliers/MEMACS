function FF = FarFieldDatabase(patternName,params)
% FARFIELDDATABASE contains several analytical patterns.

load EMconstants
f0 = 1e9;
if nargin < 1
    patternName = 'DipoleZ';
end
if nargin < 2
    params = [];
end

[ph,th] = PhThGrid;
if isfield(params,'freq')
    freq = params.freq;
else
    freq = f0;
end

switch patternName
    case 'DipoleZ'
        % Half wave z-directed dipole
        Gth = cos(pi./2.*cos(th))./sin(th);
        Gth(isinf(Gth)) = 1;
        Gph = zeros(size(Gth));
        FF = FarField(ph,th,Gth,Gph,freq);
    case 'DipoleKolitsidas'
        % Horizontal dipole from the Kolitsidas MNRAS2018 beam modelling paper
        if isfield(params,'h'), h = params.h; else, h = c0/150e6/2; end
        if isfield(params,'freq'), freq = params.freq; else, freq = 150e6; end
        lam = c0./freq;
        k = 2.*pi/lam;
        kh = k.*h;
        R = (2/3 - sin(2.*kh)./(2.*kh) - cos(2.*kh)./(2.*kh).^2 + sin(2*kh)./(2.*kh).^3);
        D = 4.*sin(kh.*cos(th)).*(1 - sin(th).^2.*sin(ph).^2)./R;
        % Fix back radiation to be smooth
        D(th >= pi/2) = max(max(D)).*lin10(-100);
%         thVect = unique(th);
%         Nback = length(thVect(thVect > pi/2));
%         Dback = D(abs(th - pi/2) < 1e-12);
%         Dback = repmat(Dback,Nback,1);
%         D(th > pi/2) = Dback;
        FF = FarField.farFieldFromPowerPattern(ph,th,D,freq);
    otherwise
        error(['Unknown patternName: ', patternName])
end