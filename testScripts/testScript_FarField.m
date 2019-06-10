% Test script for the FarField class
% Created: 2019-05-08, Dirk de Villiers
% Updated: 2019-05-24, Dirk de Villiers

close all
clearvars

disp('-------------------------------------------------------------------')
disp('...Testing FarField...');

p = mfilename('fullpath');
[filepath] = fileparts(p);
dataPath = [filepath,'\..\data\'];

%% Constructors - just run them through
try
    FF = FarField;
    disp('Pass: constructor')
    constructorPass = true;
    clear FF
catch constructor_errInfo
    disp('FAIL: constructors')
    constructorPass = false;
end

try
    FF = FarField.readCSTffs([dataPath,'CircWG_origin']);
    disp('Pass: readCSTffs')
    readCSTffsPass = true;
    clear FF
catch readCSTffs_errInfo
    disp('FAIL: readCSTffs')
    readCSTffsPass = false;
end

% % ReadNFSscan
% output = 'E1';
% outputType = 'imag';
% FFth180ph360 = FarField.readNFSscan([dataPath,'NFSscan_th180ph360']);
% figure, FFth180ph360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
% FFth360ph180 = FarField.readNFSscan([dataPath,'NFSscan_th360ph180']);
% figure, FFth360ph180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
% FFaz360el180 = FarField.readNFSscan([dataPath,'NFSscan_az360el180']);
% figure, FFaz360el180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
% FFaz180el360 = FarField.readNFSscan([dataPath,'NFSscan_az180el360']);
% figure, FFaz180el360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
% 
% FFphth = FFth360ph180.coor2Ludwig3(1);
% figure, FFphth.plot('plotType','2D','step',1,'showGrid',false,'output',output,'outputType',outputType,'scaleMag','lin')
% FFazel = FFaz180el360.coor2Ludwig3(1);
% figure, FFazel.plot('plotType','2D','step',1,'showGrid',false,'output',output,'outputType',outputType,'scaleMag','lin')
% 
% stepDeg = 1;
% xylimsDeg = [-180,180;0,180]; 
% FFphth = FFphth.currentForm2Base(stepDeg,xylimsDeg);
% FFazel = FFazel.currentForm2Base(stepDeg,xylimsDeg);
% FFdelta = FFphth - FFazel;
% figure, FFdelta.plot('plotType','2D','step',1,'showGrid',false,'output',output,'outputType','mag','scaleMag','dB')

%% Grid range
% getRange
FF = FarField.readCSTffs([dataPath,'CircWG_origin']);
yRangeNew = [0,pi/2];
xRangeNew = [0,deg2rad(355)];
FFn = FF.getRange(xRangeNew,yRangeNew);
if all(FFn.xRange == xRangeNew) && all(FFn.yRange == yRangeNew)
    disp('Pass: getRange')
    getRangePass = true;
else
    disp('FAIL: getRange')
    getRangePass = false;
end


FF = FarField.readCSTffs([dataPath,'CircWG_origin']);
figure, FF.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

FFsym180 = FF.setRangeSph('sym','180');
figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

FFpos180 = FFsym180.setRangeSph('pos');
figure, FFpos180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

FFsym360 = FFpos180.setRangeSph('sym','360');
figure, FFsym360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

FFsym180 = FFsym360.setRangeSph;
figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

FFpos360 = FFpos180.setRangeSph('pos',360);
figure, FFpos360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

FFsym180 = FFpos360.setRangeSph;
figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

% Test everything there and back and mixed...
FF1 = FF.setRangeSph('sym','180');
FF0 = FF1.setRangeSph('pos','180');
FFd0 = FF0 - FF;
err1 = FFd0.norm

FF2 = FF.setRangeSph('pos','180');
FFb = FF2.setRangeSph;
FFdb = FFb - FF1;
errb = norm(FFdb)


% figure, FF.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

% FFx = FF.setXrange('sym');
% figure, FFx.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FFx = FFx.setXrange('pos');
% figure, FFx.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

% FFy1 = FF.setYrangeSph(360);
% figure, FFy1.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FFy2 = FFy1.setXrange('sym');
% figure, FFy1.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

% FF = FF.setXrange('pos');
% figure, FF.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FFy360 = FF.setYrangeSph(360);
% figure, FFy360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% 
% FFy180 = FFy360.setYrangeSph(180);
% figure, FFy180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')


%% Final test
FarFieldPass = all([constructorPass,readCSTffsPass,...
    getRangePass,...
    ]);
if FarFieldPass
    disp('Pass: FarField');
else
    disp('FAIL: FarField');
end

disp('-------------------------------------------------------------------')


