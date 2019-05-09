% Test the symmetry implementation of the FarField class
close all
clear all

%% Get something to work with
pathName = 'Farfield Source [1]';
FF = FarField.readCSTffs(pathName);
FF = FF.coor2Ludwig3(false);
FF = FF.setXrange('sym');
FF = FF.currentForm2Base;

% Plot it
scaleMag = 'lin';
outputType = 'real';
output = 'E2';
plotGridHandle = @grid2PhTh;
FFplot = plotGridHandle(FF);
FFplot.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)

%% Get a version symmetric about x
% phXind = FF.phBase >= pi;
% phXind = FF.phBase >= -eps;
% FFx = FarField(FF.phBase(phXind),FF.thBase(phXind),FF.E1(phXind,:),FF.E2(phXind,:),FF.E3(phXind,:),FF.freq,FF.Prad./2,FF.radEff,FF.coorSys,FF.polType,'PhTh',FF.freqUnit);
% FFx = FFx.setSymmetryXZ('electric');
% FFx = FFx.setXrange('sym');
% FFx = plotGridHandle(FFx);
% figure
% FFx.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% % Fix to full range and plot
% FFxFull = FFx.mirrorSymmetricPattern;
% FFxFull = plotGridHandle(FFxFull);
% figure
% FFxFull.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% 
% FFdelta = FFxFull - FF;
% figure
% FFdelta.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType','mag','scaleMag',scaleMag)
% [eEx] = FFdelta.norm

%% Get a version symmetric about y
phYind = FF.phBase >= -pi/2-eps & FF.phBase <= pi/2+eps;
FFy = FarField(FF.phBase(phYind),FF.thBase(phYind),FF.E1(phYind,:),FF.E2(phYind,:),FF.E3(phYind,:),FF.freq,FF.Prad./2,FF.radEff,FF.coorType,FF.polType,'PhTh',FF.freqUnit);
FFy = FFy.setSymmetryYZ('magnetic');
FFy = FFy.setXrange('sym');
FFy = plotGridHandle(FFy);
figure
FFy.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% Fix to full range and plot
FFyFull = FFy.mirrorSymmetricPattern;
FFyFull = plotGridHandle(FFyFull);
figure
FFyFull.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)

FFdelta = FFyFull - FF;
figure
FFdelta.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType','mag','scaleMag',scaleMag)
[eEy] = FFdelta.norm

%% Get a version with 2 symmetry planes
% phXYind = FF.phBase >= -eps & FF.phBase <= pi/2+eps;
% FFxy = FarField(FF.phBase(phXYind),FF.thBase(phXYind),FF.E1(phXYind,:),FF.E2(phXYind,:),FF.E3(phXYind,:),FF.freq,FF.Prad./2,FF.radEff,FF.coorSys,FF.polType,'PhTh',FF.freqUnit);
% FFxy = FFxy.setSymmetryXZ('electric');
% FFxy = FFxy.setSymmetryYZ('magnetic');
% FFxy = FFxy.setXrange('sym');
% FFxy = plotGridHandle(FFxy);
% figure
% FFxy.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% % Fix to full range and plot
% FFxyFull = FFxy.mirrorSymmetricPattern;
% FFxyFull = plotGridHandle(FFxyFull);
% figure
% FFxyFull.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
% 
