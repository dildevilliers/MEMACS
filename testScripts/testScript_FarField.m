% Test script for the FarField class
% Created: 2019-05-08, Dirk de Villiers
% Updated: 2019-05-24, Dirk de Villiers

close all
clearvars

disp('-------------------------------------------------------------------')
disp('...Testing FarField...');

p = mfilename('fullpath');
[filepath] = fileparts(p);
dataPath = [filepath,'\..\data\SimPatterns\'];

%% Constructors - just run them through
tol = 1e-10;
% Hand constructors
try
    FF = FarField;
    errConstrEmpty = abs(FF.Prad - FF.pradInt);
catch constructorEmpty_errInfo
    errConstrEmpty = false;
end

try
    FF1 = FarField(FF.x,FF.y,FF.E1,FF.E2);
    errConstrBasic = abs(FF.Prad - FF.pradInt);
catch constructorBasic_errInfo
    errConstrBasic = true;
end

try
    FF2 = FarField(FF.x,FF.y,FF.E1,FF.E2,1);
    errConstrPower = abs(FF.Prad - FF.pradInt);
catch constructorPower_errInfo
    errConstrPower = true;
end

constructTestVect = [errConstrEmpty,errConstrBasic,errConstrPower];
if all(constructTestVect < tol)
    disp('Pass: constructor')
    constructorPass = true;
    clear FF FF1 FF2
else
    disp('FAIL: constructor')
    constructorPass = false;
end

% GRASP constructors
try
    FF1 = FarField.readGRASPgrd([dataPath,'FF_phth_spherical_pos180']);
    FF2 = FarField.readGRASPgrd([dataPath,'FF_uv_spherical']);
    FF3 = FarField.readGRASPgrd([dataPath,'FF_azel_spherical_pos180']);
    FF4 = FarField.readGRASPgrd([dataPath,'FF_elaz_spherical_pos180']);
    FF5 = FarField.readGRASPgrd([dataPath,'FF_trueview_spherical']);
    FF6 = FarField.readGRASPgrd([dataPath,'FF_phth_ludwig3_pos180']);
    FF7 = FarField.readGRASPgrd([dataPath,'FF_phth_circular_pos180']);
    disp('Pass: readGRASPgrd')
    readGRASPgrdPass = true;
    clear FF1 FF2 FF3 FF4 FF5 FF6 FF7
catch readCSTffs_errInfo
    disp('FAIL: readGRASPgrd')
    readGRASPgrdPass = false;
end

try
    Nf = 11;
    Ncut = 37;
    f = linspace(1,1.5,Nf);
    FF1 = FarField.readGRASPcut([dataPath,'FFcut_spherical'],Nf,Ncut,'freq',f,'freqUnit','GHz');
    FF2 = FarField.readGRASPcut([dataPath,'FFcut_ludwig3'],Nf,Ncut,'freq',f,'freqUnit','GHz');
    FF3 = FarField.readGRASPcut([dataPath,'FFcut_circular'],Nf,Ncut,'freq',f,'freqUnit','GHz');
    disp('Pass: readGRASPcut')
    readGRASPcutPass = true;
    clear FF1 FF2 FF3
catch readCSTffs_errInfo
    disp('FAIL: readGRASPcut')
    readGRASPcutPass = false;
end

% CST constructor
try
    FF = FarField.readCSTffs([dataPath,'CircWG_origin']);
    disp('Pass: readCSTffs')
    readCSTffsPass = true;
    clear FF
catch readCSTffs_errInfo
    disp('FAIL: readCSTffs')
    readCSTffsPass = false;
end

%% Field normalization
% pradInt - use the known power from GRASP as test
tol = 1e-3;
FF1 = FarField.readGRASPgrd([dataPath,'FF_phth_spherical_pos180']);
errIntPhTh = abs(FF1.Prad - FF1.pradInt)./FF1.Prad;
FF2 = FarField.readGRASPgrd([dataPath,'FF_azel_spherical_sym180']);
errIntAzEl = abs(FF2.Prad - FF2.pradInt)./FF2.Prad;

pradIntTestVect = [errIntPhTh,errIntAzEl];
if all(pradIntTestVect < tol)
    disp('Pass: pradInt')
    pradIntPass = true;
    clear FF1 FF2 FF3
else
    disp('FAIL: pradInt')
    pradIntPass = false;
end

% setPower
FF1 = FarField.readGRASPgrd([dataPath,'FF_phth_spherical_pos180']);
FF2 = FF1.setPower(1);
errSetPower = abs(FF2.pradInt - 1);

if all(errSetPower < tol)
    disp('Pass: setPower')
    setPowerPass = true;
    clear FF1 FF2
else
    disp('FAIL: setPower')
    setPowerPass = false;
end

%% Grid transformations
showGridPlots = true;

gridLocal = {'PhTh','AzEl','ElAz'};
gridProj = {'DirCos','TrueView','ArcSin'}; 
gridAstro = {'Horiz','RADec','GalLongLat'};

% Test PhTh input to all local and projections
FF1 = FarField.readGRASPgrd([dataPath,'FF_phth_spherical_pos180']);
gridPhThTest = [gridLocal(2:end),gridProj];
for ii = 1:length(gridPhThTest)
    FF2 = FF1.changeGrid(gridPhThTest{ii});
    if showGridPlots && 0
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2PhTh;
    gridTransPhThTestVect(ii) = isGridEqual(FF3,FF1);
end

% Test AzEl input to all local and projections
FF1 = FarField.readGRASPgrd([dataPath,'FF_azel_spherical_sym180']);
gridAzElTest = [gridLocal([1,3]),gridProj];
for ii = 1:length(gridAzElTest)
    FF2 = FF1.changeGrid(gridAzElTest{ii});
    if showGridPlots && 0
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2AzEl;
    gridTransAzElTestVect(ii) = isGridEqual(FF3,FF1);
end

% Test ElAz input to all local and projections
FF1 = FarField.readGRASPgrd([dataPath,'FF_elaz_spherical_sym180']);
gridElAzTest = [gridLocal([1,2]),gridProj];
for ii = 1:length(gridElAzTest)
    FF2 = FF1.changeGrid(gridElAzTest{ii});
    if showGridPlots
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2ElAz;
    gridTransElAzTestVect(ii) = isGridEqual(FF3,FF1);
end

gridTransTestVect = [gridTransPhThTestVect,gridTransAzElTestVect,gridTransElAzTestVect]

% % Loop through the combos, and build up the TestVect
% % Also read in other base grids from GRASP and loop those
% FF2 = FF1.grid2DirCos;
% FF3 = FF2.grid2PhTh;
% 
% gridTransTestVect = isGridEqual(FF3,FF1)

keyboard

% % ReadNFSscan
% output = 'E1';
% outputType = 'phase';
% FFth180ph360 = FarField.readNFSscan([dataPath,'NFSscan_th180ph360']);
% figure, FFth180ph360.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin'), axis normal
% FFth360ph180 = FarField.readNFSscan([dataPath,'NFSscan_th360ph180']);
% figure, FFth360ph180.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin'), axis normal
% FFaz360el180 = FarField.readNFSscan([dataPath,'NFSscan_az360el180']);
% figure, FFaz360el180.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin'), axis normal
% FFaz180el360 = FarField.readNFSscan([dataPath,'NFSscan_az180el360']);
% figure, FFaz180el360.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin'), axis normal
% 
% FFphth = FFth360ph180.coor2Ludwig3(1);
% figure, FFphth.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin'), axis normal
% FFazel = FFaz180el360.coor2Ludwig3(1);
% figure, FFazel.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin'), axis normal
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
outputType = 'real';
output = 'E2';
% setRange - plot a variety and calculate the differences
% First PhTh
FF = FarField.readCSTffs([dataPath,'CircWG_origin']);
figure, FF.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
% Get a tolerance from the field values
tol = 1e-6*min(sqrt(lin10(FF.Directivity_dBi)));

% p180 -> s180 -> p180
FFsym180 = FF.setRangeSph('sym','180');
figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
% FFpos180 = FFsym180.setRangeSph('pos');
% figure, FFpos180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
FFs180 = FFsym180.setRangeSph('pos') - FF;
err_s180 = norm(norm(FFs180));

% p180 -> s360 -> p180
FFsym360 = FF.setRangeSph('sym','360');
figure, FFsym360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
FFsym180 = FFsym360.setRangeSph;
figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
FFs360 = FFsym360.setRangeSph('pos') - FF;
err_s360 = norm(norm(FFs360));

% p180 -> p360 -> p180
FFpos360 = FF.setRangeSph('pos',360);
figure, FFpos360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
% FFsym180 = FFpos360.setRangeSph;
% figure, FFsym180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
FFp360 = FFpos360.setRangeSph('pos') - FF;
err_p360 = norm(norm(FFp360));

% Now the azel/elaz
FF1 = FF.coor2Ludwig2AE;
FF1 = FF1.currentForm2Base;
figure, FF1.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')

% s180 -> p180 -> s180
FF1pos180 = FF1.setRangeSph('pos','180');
figure, FF1pos180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
% FF1sym180 = FF1pos180.setRangeSph('sym','180');
% figure, FF1sym180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
FF1s180 = FF1pos180.setRangeSph - FF1;
err1_s180 = norm(norm(FF1s180));

% p180 -> s360 -> s180
FF1sym360 = FF1.setRangeSph('sym','360');
figure, FF1sym360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
% FF1sym180 = FF1sym360.setRangeSph;
% figure, FF1sym180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
FF1s360 = FF1sym360.setRangeSph - FF1;
err1_s360 = norm(norm(FF1s360));

% s180 -> p360 -> s180
FF1pos360 = FF1.setRangeSph('pos','360');
figure, FF1pos360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
% FF1sym180 = FF1pos360.setRangeSph;
% figure, FF1sym180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag','lin')
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
FarFieldPass = all([constructorPass,readGRASPgrdPass,readGRASPcutPass,readCSTffsPass,...
    pradIntPass,setPowerPass,...
    getRangePass,setRangePass,...
    ]);
if FarFieldPass
    disp('Pass: FarField');
else
    disp('FAIL: FarField');
end

disp('-------------------------------------------------------------------')


