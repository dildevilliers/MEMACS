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

% setRange - plot a variety and calculate the differences
% First PhTh
FF = FarField.readCSTffs([dataPath,'CircWG_origin']);
figure, FF.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% Get a tolerance from the field values
tol = 1e-6*min(sqrt(lin10(FF.Directivity_dBi)));

% p180 -> s180 -> p180
FFsym180 = FF.setRangeSph('sym','180');
figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FFpos180 = FFsym180.setRangeSph('pos');
% figure, FFpos180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
FFs180 = FFsym180.setRangeSph('pos') - FF;
err_s180 = norm(norm(FFs180));

% p180 -> s360 -> p180
FFsym360 = FF.setRangeSph('sym','360');
figure, FFsym360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FFsym180 = FFsym360.setRangeSph;
% figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
FFs360 = FFsym360.setRangeSph('pos') - FF;
err_s360 = norm(norm(FFs360));

% p180 -> p360 -> p180
FFpos360 = FF.setRangeSph('pos',360);
figure, FFpos360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FFsym180 = FFpos360.setRangeSph;
% figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
FFp360 = FFpos360.setRangeSph('pos') - FF;
err_p360 = norm(norm(FFp360));

% Now the azel/elaz
FF1 = FF.coor2Ludwig2AE;
FF1 = FF1.currentForm2Base;
figure, FF1.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')

% s180 -> p180 -> s180
FF1pos180 = FF1.setRangeSph('pos','180');
figure, FF1pos180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FF1sym180 = FF1pos180.setRangeSph('sym','180');
% figure, FF1sym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
FF1s180 = FF1pos180.setRangeSph - FF1;
err1_s180 = norm(norm(FF1s180));

% p180 -> s360 -> s180
FF1sym360 = FF1.setRangeSph('sym','360');
figure, FF1sym360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FF1sym180 = FF1sym360.setRangeSph;
% figure, FF1sym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
FF1s360 = FF1sym360.setRangeSph - FF1;
err1_s360 = norm(norm(FF1s360));

% s180 -> p360 -> s180
FF1pos360 = FF1.setRangeSph('pos','360');
figure, FF1pos360.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
% FF1sym180 = FF1pos360.setRangeSph;
% figure, FF1sym180.plot('plotType','2D','step',1,'showGrid',1,'output','E1','outputType','real','scaleMag','lin')
FF1p360 = FF1pos360.setRangeSph - FF1;
err1_p360 = norm(norm(FF1p360));

setTestVect = [err_s180,err_s360,err_p360,err1_s180,err1_s360,err1_p360];
if all(setTestVect < tol)
    disp('Pass: setRange')
    setRangePass = true;
else
    disp('FAIL: setRange')
    setRangePass = false;
end



%% Final test
FarFieldPass = all([constructorPass,readCSTffsPass,...
    getRangePass,setRangePass,...
    ]);
if FarFieldPass
    disp('Pass: FarField');
else
    disp('FAIL: FarField');
end

disp('-------------------------------------------------------------------')


