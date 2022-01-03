% Script to test FarField rotations
clear all
close all

p = mfilename('fullpath');
[filepath] = fileparts(p);
dataPath = [filepath,'\..\data\SimPatterns\rotateFF\'];

plotDim = 3;    % Select 1, 2 or 3 for 2D or 3D plots
grid2Dplot = 'PhTh'; % Can be DirCos, TrueView, PhTh, etc
coorPlot = 'Ludwig3';
output = 'Directivity';
outputType = 'mag';
scaleMag = 'dB';
fi = 6;
onlyPower = true;

% GRASP rotation
GRASPangDeg = [65,-110,-15];
cGlob = CoordinateSystem;
cRot = cGlob.rotGRASP(deg2rad(GRASPangDeg));

% Read a test field
pathName = [dataPath,'FForigin'];
FF = FarField.readGRASPgrd(pathName);

handle2Dgrid = str2func(['grid2',grid2Dplot]);
handleCoor = str2func(['coor2',coorPlot]);
% Plot the original version
if plotDim == 3
    FF.plot('plotType','3D','output',output)
elseif plotDim == 2
    FF = handle2Dgrid(FF);
    FF = handleCoor(FF,false);
    FF.plot('plotType','2D','output',output,'outputType',outputType,'step',1,'showGrid',true,'dynamicRange_dB',40,'scaleMag',scaleMag,'freqIndex',fi), hold on
elseif plotDim == 1
    FF = handleCoor(FF,false);
    FF.plotPrincipleCuts('output',output);
end

% Rotate the field
% FFr = FF.rotate(rotHandle,rotAng);
FFr = FF.rotate(cRot);

if plotDim == 3
    figure
    FFr.plot('plotType','3D','output',output)
elseif plotDim == 2
    figure
    FFr = handle2Dgrid(FFr);
    FFr = handleCoor(FFr,false);
    FFr.plot('plotType','2D','output',output,'outputType',outputType,'step',1,'showGrid',true,'dynamicRange_dB',40,'scaleMag',scaleMag,'freqIndex',fi), hold on
elseif plotDim == 1
    FFr = handleCoor(FFr,false);
    FFr.plotPrincipleCuts('output',output);
end

% Read the rotated validation file
FFrotVal = FarField.readGRASPgrd([dataPath,'FFrot']);

if plotDim == 3
    figure
    FFrotVal.plot('plotType','3D','output',output)
elseif plotDim == 2
    figure
    FFrotVal = handle2Dgrid(FFrotVal);
    FFrotVal = handleCoor(FFrotVal,false);
    FFrotVal.plot('plotType','2D','output',output,'outputType',outputType,'step',1,'showGrid',true,'dynamicRange_dB',40,'scaleMag',scaleMag,'freqIndex',fi), hold on
elseif plotDim == 1
    FFrotVal = handleCoor(FFrotVal,false);
    FFrotVal.plotPrincipleCuts('output',output);
end

%% Plot the difference field
FFdel = FFrotVal - FFr;
FFdel = handle2Dgrid(FFdel);
FFdel = handleCoor(FFdel,false);
figure
FFdel.plot('plotType','2D','output',output,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40,'scaleMag',scaleMag,'freqIndex',fi), hold on

