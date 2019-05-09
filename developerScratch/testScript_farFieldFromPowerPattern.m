clear all
close all

Nth_cut = 37;
Nph_cut = 73;
th = linspace(0,pi,37);
ph = linspace(0,2*pi,73);
th0 = 45;
taper_dB = -10;
freq = 1e9;

[PH,TH] = meshgrid(ph,th);

%Create Gaussian power pattern to build FarField object with
P1 = powerPattern(PH(:),TH(:),'gauss',th0,taper_dB,freq);

%Build example FarField objects from P1 power pattern (different
%polarizations)

[plotlinX plotlinY plotcircLH plotcircRH] = deal([1 1 1 1]);

if plotlinX
    FF1 = FarField.farFieldFromPowerPattern(PH(:),TH(:),P1,freq,'fieldPol','linearX');
    figure
    %E1 is -inf dB (perfectly X-linear field in Ludwig3 coords) and will break plot
    subplot(3,1,1)
    FF1.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    FF1 = coor2spherical(FF1);
    subplot(3,1,2)
    FF1.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    subplot(3,1,3)
    FF1.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    % FF1.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30)
end

if plotlinY
    FF2 = FarField.farFieldFromPowerPattern(PH(:),TH(:),P1,freq,'fieldPol','linearY');
    figure
    %E1 is -inf dB (perfectly Y-linear field in Ludwig3 coords) and will break plot
    % FF2.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    subplot(3,1,1)
    FF2.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    FF2 = coor2spherical(FF2);
    subplot(3,1,2)
    FF2.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    subplot(3,1,3)
    FF2.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
end

if plotcircLH
    FF3 = FarField.farFieldFromPowerPattern(PH(:),TH(:),P1,freq,'fieldPol','circularLH');
    figure
    %E2 is -inf dB (perfectly LH circular field in spherical pol basis) and will break plot
    subplot(3,1,1)
    FF3.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    FF3 = pol2linear(FF3);
    subplot(3,1,2)
    FF3.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    subplot(3,1,3)
    FF3.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
end

if plotcircRH
    FF4 = FarField.farFieldFromPowerPattern(PH(:),TH(:),P1,freq,'fieldPol','circularRH');
    figure
    %E1 is -inf dB (perfectly RH circular field in spherical pol basis) and will break plot
    subplot(3,1,1)
    FF3.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    FF3 = pol2linear(FF3);
    subplot(3,1,2)
    FF3.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
    subplot(3,1,3)
    FF3.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
end
