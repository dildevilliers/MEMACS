% Script to test FarField rotations
clear all
close all

plotDim = 3;    % Select 1, 2 or 3 for 2D or 3D plots
grid2Dplot = 'PhTh'; % Can be DirCos, TrueView, PhTh, etc
coorPlot = 'spherical';
output = 'Directivity';
onlyPower = true;

% CST like axis rotations
rotX = deg2rad(45);
rotY = deg2rad(45);
rotZ = deg2rad(45);
cGlob = coordinateSystem;
cRot = cGlob.rotX(rotX);
cRot = cRot.rotY(rotY);
cRot = cRot.rotZ(rotZ);
rotAng = cRot.getGRASPangles;

rotHandle = @rotGRASP;
% rotAng = deg2rad([45,45,45]);
% rotAng = deg2rad([0,0,90]);

% Read a test field
pathName = 'CircWG_origin';
FF = FarField.readCSTffs(pathName);

handle2Dgrid = str2func(['grid2',grid2Dplot]);
handleCoor = str2func(['coor2',coorPlot]);
% Plot the original version
if plotDim == 3
    FF.plot('plotType','3D','output',output)
elseif plotDim == 2
    FF = handle2Dgrid(FF);
    FF = handleCoor(FF,false);
    FF.plot('plotType','2D','output',output,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40), hold on
elseif plotDim == 1
    FF = handleCoor(FF,false);
    FF.plotPrincipleCuts('output',output);
end

% Rotate the field
FFr = FF.rotate(rotHandle,rotAng,onlyPower);
if plotDim == 3
    figure
    FFr.plot('plotType','3D','output',output)
elseif plotDim == 2
    figure
    FFr = handle2Dgrid(FFr);
    FFr = handleCoor(FFr,false);
    FFr.plot('plotType','2D','output',output,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40), hold on
elseif plotDim == 1
    FFr = handleCoor(FFr,false);
    FFr.plotPrincipleCuts('output',output);
end

% Read the rotated validation file
FFrotVal = FarField.readCSTffs('CircWG_rot');
if plotDim == 3
    figure
    FFrotVal.plot('plotType','3D','output',output)
elseif plotDim == 2
    figure
    FFrotVal = handle2Dgrid(FFrotVal);
    FFrotVal = handleCoor(FFrotVal,false);
    FFrotVal.plot('plotType','2D','output',output,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40), hold on
elseif plotDim == 1
    FFrotVal = handleCoor(FFrotVal,false);
    FFrotVal.plotPrincipleCuts('output',output);
end

%% Plot the difference field
FFdel = FFrotVal - FFr;
FFdel = handle2Dgrid(FFdel);
FFdel = handleCoor(FFdel,false);
figure
FFdel.plot('plotType','2D','output',output,'outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40), hold on

