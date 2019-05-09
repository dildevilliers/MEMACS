% Test script for the suite of Array classes

close all
clear all

%%
fRF = 1.57542e9;
th_in = deg2rad(0);
ph_in = deg2rad(0);
Ps = -60; % in dBm
phaseSig = deg2rad(0);

Nant = 12;      % Number of elements
d_lam = 0.5;     % Element spacing in wavelengths 
d = d_lam*(physconst('lightspeed')/fRF);    % Element spacing in m
x = (-(Nant-1)/2:(Nant-1)/2).*d;       % Uniform linear array along x-axis
r = pnt3D(x,0,0);
% ampErrors = [1,1,1,1];
% phaseErrors_deg = [0,0,0,0];
ampErrors = ones(1,Nant);
phaseErrors_deg = zeros(1,Nant);
phaseErrors_deg([1,2,3,4,7,9]) = [0,15,90,45,15,50];
ampErrors([1,4,5,6,7]) = [2,0.1,0.7,0.9,1];
channelPhasors = ampErrors.*exp(1i.*deg2rad(phaseErrors_deg));

Pn = -80; % in dBm
LNAGain = 20;
IFGain = 60;
fIF = 4.096e6;
fLO = fRF - fIF; 

fSamp = 4*fIF;          % Sample rate in Hz
Nt = 8000;                % Number of time samples
delT = 1/fSamp;
t0 = 0;
t = t0:delT:(t0+delT*(Nt-1));

Nbits = 1;
A = ArrayADC(Nbits);

S1 = PlaneWaveSignal('compExp',fRF,th_in,ph_in,Ps,phaseSig);
S2 = PlaneWaveSignal('compExp',fRF,th_in + deg2rad(40),ph_in,Ps,phaseSig+pi/2);
S3 = PlaneWaveSignal('noise',fRF,th_in + deg2rad(-20),ph_in,Ps,phaseSig);
% S = [S1,S2,S3];
S = S1;
% E = ArrayElements(r);
% R = ArrayReceiver(Pn,LNAGain,IFGain,fLO);
% portSigMat = E.portSignals(S,t);
% sn = R.sigRec(portSigMat,t);
% % b = A.getBinarySignal(real(sn));
% % x = A.getDigitalSignal(b);
% x = A.ADC(real(sn));
% plot(t,x(1,:))
% 
% plot(t,real(sn(1,:)),'k'), grid on, hold on

arraySys = ArraySystem(r,channelPhasors,Pn,LNAGain,IFGain,fLO,fSamp,Nt,Nbits);
arraySys.plotPortSignal(S);
x = arraySys.getPortSignal(S);

%% Plot signal phases
% figure
% plot(t*1e6,phase(x(1,:))), grid on, hold on
% plot(t*1e6,phase(x(2,:)))
% plot(t*1e6,phase(x(3,:)))
% plot(t*1e6,phase(x(4,:)))
% xlabel('Time (\mu s)')
% ylabel('Phase (deg)')

% % Extract the phase from the received signals
% r = arraySys.elements.antPos.pointMatrix;
% 
% Pmean = rad2deg(mean(angle(x),2));
% figure
% plot(Pmean), grid on
% 
% figure
% plot(real(x(1,:)),imag(x(1,:)),'*'), grid on
y = zeros(Nant-1,length(t));
for ii = 1:Nant-1
    y(ii,:) = x(ii+1,:)./x(ii,:);
end
ym = mean(rad2deg(angle(y)),2);
err = sum(abs(ym));
figure
plot(ym,'o'), grid on


arrayDBE = ArrayDBE(arraySys);
% figure
% arrayDBE.plotPortPSD(x,1,4)
% 
% Test the scanning
Nscan = 501;
th_scan_deg = linspace(-60,60,Nscan);
% ph_scan_deg = ones(size(th_scan_deg))*0;    % Use ph = 0 for array along x-axis
ph_scan_deg = 0;    % Use ph = 0 for array along x-axis
% y = arrayDBE.scanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x);
figure
arrayDBE.plotScanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x);


% angle_deg = -180:2:180;
% 
% z = zeros(1,181);
% for ii = 1:Nant-1
%     for jj = 1:181
%         y(ii,:) = mean(rad2deg(angle(x(ii+1,:)./x(ii,:)) + angle_deg(1,jj))) ;
%         z(1,jj) = abs(mean(y(ii,:)));
%     end
%     [min_z, index] = min(z(1,:))
% %     figure
% %     a = mean(rad2deg(angle(x(ii+1,:))))
% %     plot(mean(rad2deg(angle(x(ii+1,:)))))
% %     hold on
%     x(ii+1,:) = x(ii+1,:).*exp(-1i*deg2rad(min_z));
% %     b = mean(rad2deg(angle(x(ii+1,:))))
% %     plot(mean(rad2deg(angle(x(ii+1,:)))))
% end



xMod = x;
delP = zeros(1,Nant);
for aa = 2:Nant
    delP(aa) = mean((angle(xMod(aa,:)./xMod(aa-1,:))));
    xMod(aa,:) = xMod(aa,:).*exp(-1i.*delP(aa));
end

calVect = zeros(1,Nant);
for aa = 2:Nant
    calVect(aa) = mean((angle(x(aa,:)./x(1,:))));
end


calVector = ones(1,Nant);
for aa = 2:Nant
    calVector(aa) = mean((x(aa,:)./x(1,:)));
end

y = zeros(Nant-1,length(t));
for ii = 1:Nant-1
    y(ii,:) = x(ii+1,:)./x(ii,:);
end

ym = mean(rad2deg(angle(y)),2);
err = sum(abs(ym));
figure
plot(ym,'o'), grid on

fileName = 'Channel signals';
fidDPA = fopen(fileName,'wt');
fprintf([repmat('%f\t', 1, size(x, 2)) '\n'], x')
fclose(fidDPA);
            
%%
% 
arrayDBE = ArrayDBE(arraySys);
% figure
% arrayDBE.plotPortPSD(x,1,4)
% 
% Test the scanning
Nscan = 501;
th_scan_deg = linspace(-60,60,Nscan);
% ph_scan_deg = ones(size(th_scan_deg))*0;    % Use ph = 0 for array along x-axis
ph_scan_deg = 0;    % Use ph = 0 for array along x-axis
% y = arrayDBE.scanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x);
figure
arrayDBE.plotScanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),xMod);

figure
arrayDBE.plotScanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x,delP);

figure
arrayDBE.plotScanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x,calVect);

figure
arrayDBE.plotScanBeam(fRF,deg2rad(th_scan_deg),deg2rad(ph_scan_deg),x,calVector);