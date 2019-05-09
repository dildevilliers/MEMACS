function G = feedPatterns(th0,taper_dB,th,freq,r)

% function G = feedPatterns(th0,taper_dB,th,freq,r)
% computes the cos^n(th/2), cos^n(th), simple tapered as well as gaussian beam feed
% patterns (See GRASP manual for a description of the last two)
% Returns a structure containing the patterns as a function of th as
% Gcosn2, Gcosn, Gsimple, and Ggauss (E-field patterns are returned - not power)
% All the patterns are normalized to 1
% Inputs:
% th0 = subtended half angle in rad
% taper_dB = required edge taper in dB (typically  < 0)
% th = theta vector
% freq = frequency in Hz
% r = nearfield distance (optional)

if nargin < 5
    ff_force = 1;
elseif nargin == 5
    ff_force = 0;
end

%% Calculate the cos^n(th/2) pattern
nc2 = taper_dB./(20*log10(cos(th0/2)));
Gcosn2 = cos(th./2).^nc2;
G.nc2 = nc2;

%% Calculate the cos^n(th) pattern
nc = taper_dB./(20*log10(cos(th0)));
Gcosn = cos(th).^nc;
G.nc = nc;

%% Calculate simple tapered pattern
Gsimple = 10.^(taper_dB./20.*(th./th0).^2);

%% Calculate the gaussian beam pattern
c0 = 299792458;
lambda = c0./freq;
k = 2.*pi./lambda;
b = (20.*log10((1 + cos(th0))./2) - taper_dB)./(20.*k.*(1-cos(th0)).*log10(exp(1)));  % GRASP Tech (2.3-24)
if ff_force
    Ggauss = exp(k.*b.*cos(th)).*(1 + cos(th));
else
    rho = r.*sin(th);
    z = r.*cos(th);
    Ggauss = 1./sqrt(rho.^2 + z.^2 - b.^2 + 1i.*2.*b.*z).*exp(-1i.*k.*sqrt(rho.^2 + z.^2 - b.^2 + 1i.*2.*b.*z)).*(1 + cos(th));
end
G.b = b;

G.Gcosn2 = Gcosn2./max(Gcosn2);
G.Gcosn = Gcosn./max(Gcosn);
G.Gsimple = Gsimple./max(Gsimple); 
G.Ggauss = Ggauss./max(Ggauss); 
