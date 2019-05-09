function Pout = powerPattern(ph,th,Pin,th0,taper_dB,freq)

%Name: powerPattern.m
%Description:
%   Function to create a power pattern vector to be used in
%   farFieldFromPowerPattern function for the @FarField class
%Inputs:
% --ph: column vector [Nang x 1] of ph angles in rad - Nang is no. of unique angles
% --th: column vector [Nang x 1] of th angles in rad - Nang is no. of unique angles
% --Pin: Input power pattern cut [Nth x 1]/string/cut pair [Nth x 2]
% --th0 = subtended half angle in rad
% --taper_dB = required edge taper in dB (typically  < 0)
% --fieldPol: (optional) string denoting how far field should be polarized- options are 'linearX', 'linearY', 'circularLH', 'circularRH' 
% --freq: scalar frequency in Hz
%Outputs:
% --Pout: [(Nth*Nph) x 1] column vector of output power pattern - specified
% over range described by th, ph inputs

%Basic error checking
if ( ~isvector(th) || ~isvector(ph))
    error('Incompatible th and/or ph inputs - must be vectors of same length')
end

% Define a few useful terms
Pin_pattern = Pin;
[Nth,Nph] = deal(length(th),length(ph));
assert(Nth == Nph,'ph and th must be the same length');
Nth_cut = length(unique(th));
Nph_cut = length(unique(ph));

%Generate power pattern if a string was given for Pin
if ischar(Pin)
    G = feedPatterns(th0,taper_dB,th,freq);
    switch Pin
        case 'gauss'
            Pin_pattern = G.Ggauss;
        case 'cosn2'
            Pin_pattern = G.Gcosn2;
        case 'cosn'
            Pin_pattern = G.Gcosn;
        case 'simple'
            Pin_pattern = G.Gsimple;
        otherwise
            error('Power pattern input string not recognised');
    end
end

%% Generate output pattern
[NPin1,NPin2] = size(Pin_pattern);
if (NPin1 == Nth && NPin2 == 1) %case where a full 2D power pattern was given, just return Pin
    Pout = Pin_pattern;
elseif (NPin1 == Nth_cut) %either a 1D axially symmetric pattern or a symmetric pattern with cardinal plane cuts
    switch NPin2
        case 1 %1D axially symmetric pattern
            Pout = repmat(Pin_pattern,Nph_cut,1);
        case 2 %Symmetric pattern defined by (phi = 0) and (phi = pi/2) cuts - final power pattern is sum of both projected rectangularly
            [Pout1,Pout2] = deal(repmat(Pin_pattern(:,1),Nph_cut,1),repmat(Pin_pattern(:,2),Nph_cut,1));
            Pout = Pout1.*abs(cos(ph)) + Pout2.*abs(sin(ph));
        otherwise
            error('Pin must have 1 or 2 columns')
    end
else
    error('Pin input size is incompatible with th and/or ph inputs')
end


end