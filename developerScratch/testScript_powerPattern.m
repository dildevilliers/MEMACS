Nth_cut = 37;
Nph_cut = 73;
th = linspace(0,pi,Nth_cut);
ph = linspace(0,2*pi,Nph_cut);
th0 = 45;
taper_dB = -10;
freq = 1e9;

[PH,TH] = meshgrid(ph,th);

%% Example 1 - Pin specified by string

P1 = powerPattern(TH(:),PH(:),'gauss',th0,taper_dB,freq);

figure
surf(rad2deg(TH),rad2deg(PH),reshape(P1,Nth_cut,Nph_cut))
title('Test Power Pattern 1')
xlabel('\theta (degrees)')
ylabel('\phi (degrees)')
axis tight

%% Example 2 - Pin specified as vector

P2cut = abs(sinc(3.*th)).'; 
P2 = powerPattern(PH(:),TH(:),P2cut,th0,taper_dB,freq);

figure
surf(rad2deg(TH),rad2deg(PH),reshape(P2,Nth_cut,Nph_cut))
title('Test Power Pattern 2')
xlabel('\theta (degrees)')
ylabel('\phi (degrees)')
axis tight

%% Example 3 - Pin specified as two vectors in cardinal planes

G = feedPatterns(th0,taper_dB,th,freq);
P3cut = G.Ggauss.';

P3 = powerPattern(TH(:),PH(:),[P2cut,P3cut],th0,taper_dB,freq);

figure
surf(rad2deg(TH),rad2deg(PH),reshape(P3,Nth_cut,Nph_cut))
title('Test Power Pattern 3')
xlabel('\theta (degrees)')
ylabel('\phi (degrees)')
axis tight
