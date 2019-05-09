function FF = stdFF(th,Gth,freq,polPhi,Nph)

%  function FF = stdFF(th,Gth,freq,polPhi,Nph)
% Function that returns a standard FF structure in FF for a given
% directivity pattern Gth, sampled at th (rad) [0,pi].  FF will have Nph cuts in
% phi (default 181), with polarization angle (polPhi in rad) default pi/2 ('y').
% The field is perfectly linearly polarized in direction phiPol.
% Gth may be frequency dependent, specified then as a matrix of size (Nth x Nf)
% The freq vector in freq may be left out, but then an empty frequency
% matrix is returned in the structure FF.

% Dirk de Villiers
% Create: 2014-09-17
% Last edit: 2014-09-18


if nargin == 2, [freq,Nph,polPhi] = deal([]);
elseif nargin == 3, [polPhi,Nph] = deal([]);
elseif nargin == 4, [Nph] = deal([]);
end

if isempty(Nph), Nph = 181; end
if isempty(polPhi), polPhi = pi/2; end

% if min(th) ~= 0 || abs(max(th) - pi) > eps
%     error('th should be in the range [0,pi]');
% end

Nth = length(th);
Nf = max(length(freq),1);

ph = linspace(0,2*pi,Nph);

% Check if G is frequency dependent
[Ng,Ngf] = size(Gth);
if Ng == 1 % Row vector - no frequency dependence 
    Gth = reshape(Gth,Ngf,1);
    Ngf = 1;    % number of specified directivity frequencies
end
if Ngf ~= Nf
    warning('Only the first column of Gth will be used and repeated at each specified frequency, since different numbers of frequencies are defind in the frequency vector and the Gth matrix');
    Gth = repmat(Gth(:,1),1,Nf);
end

E = sqrt(Gth);
    
% Pre-allocate
[TH,PH,Eth,Eph] = deal(zeros(Nth*Nph,Nf));
TH = repmat(reshape(th,Nth,1),Nph,Nf);
for pp = 1:Nph
    phMAT = ones(Nth,Nf).*ph(pp);
    
    rowPos = (pp-1)*Nth + 1:pp*Nth; 
    PH(rowPos,:) = phMAT;
    
    cosPh = cos(phMAT-polPhi);
    sinPh = sin(phMAT-polPhi);
%     keyboard;
    Eth(rowPos,:) = E.*cosPh;
    Eph(rowPos,:) = -E.*sinPh;
end

FF.th = TH;
FF.ph = PH;
FF.Eth = Eth;
FF.Eph = Eph;
FF.freq = freq;

