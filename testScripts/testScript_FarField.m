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
dataPathGRASPgrd = [dataPath,'GRASPgrd\'];
dataPathGRASPcut = [dataPath,'GRASPcut\'];
dataPathCSTffs = [dataPath,'CSTffs\'];
dataPathCSTtxt = [dataPath,'CSTtxt\'];


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
    FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_phth_spherical_pos180']);
    FF2 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_uv_spherical']);
    FF3 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_azel_spherical_pos180']);
    FF4 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_elaz_spherical_pos180']);
    FF5 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_trueview_spherical']);
    FF6 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_phth_ludwig3_pos180']);
    FF7 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_phth_circular_pos180']);
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
    FF1 = FarField.readGRASPcut([dataPathGRASPcut,'FFcut_spherical'],Nf,Ncut,'freq',f,'freqUnit','GHz');
    FF2 = FarField.readGRASPcut([dataPathGRASPcut,'FFcut_ludwig3'],Nf,Ncut,'freq',f,'freqUnit','GHz');
    FF3 = FarField.readGRASPcut([dataPathGRASPcut,'FFcut_circular'],Nf,Ncut,'freq',f,'freqUnit','GHz');
    disp('Pass: readGRASPcut')
    readGRASPcutPass = true;
    clear FF1 FF2 FF3
catch readCSTffs_errInfo
    disp('FAIL: readGRASPcut')
    readGRASPcutPass = false;
end

% CST constructors
try
    FF = FarField.readCSTffs([dataPathCSTffs,'FF']);
    disp('Pass: readCSTffs')
    readCSTffsPass = true;
    clear FF
catch readCSTffs_errInfo
    disp('FAIL: readCSTffs')
    readCSTffsPass = false;
end

readCSTtxtPass = true;
coorVect = {'spherical','azel','elaz','L3'};
polVect = {'lin','circ'};
xRangeVect = {'pos','sym'};
yRangeVect = {'180','360'};
for cc = 1:length(coorVect)
    for pp = 1:length(polVect)
        for xx = 1:length(xRangeVect)
            for yy = 1:length(yRangeVect)
                fileName = ['FF_',coorVect{cc},'_',polVect{pp},'_',xRangeVect{xx},yRangeVect{yy}];
                try
                    FF = FarField.readCSTtxt([dataPathCSTtxt,fileName]);
                    if 0
                        figure, FF.plot('plotType','2D','showGrid',1)
                    end
                    readCSTtxtPass = readCSTtxtPass && 1;
                catch readCSTtxt_errInfo
                    readCSTtxtPass = readCSTtxtPass && 0;
                    keyboard
                end
            end
        end
    end
end
if readCSTtxtPass
    disp('Pass: readCSTtxt')
    clear FF 
else
    disp('FAIL: readCSTtxt')
end


%% Field normalization
% pradInt - use the known power from GRASP as test
tol = 1e-3;
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_phth_spherical_pos180']);
errIntPhTh = abs(FF1.Prad - FF1.pradInt)./FF1.Prad;
FF2 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_azel_spherical_sym180']);
errIntAzEl = abs(FF2.Prad - FF2.pradInt)./FF2.Prad;
FF3 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_elaz_spherical_sym180']);
errIntElAz = abs(FF3.Prad - FF3.pradInt)./FF3.Prad;
% Also integrate the CST fields - should all give the same answer
PradCST = zeros(1,length(coorVect)*length(polVect)*length(xRangeVect)*length(yRangeVect));
ii = 1;
for cc = 1:length(coorVect)
    for pp = 1:length(polVect)
        for xx = 1:length(xRangeVect)
            for yy = 1:length(yRangeVect)
                fileName = ['FF_',coorVect{cc},'_',polVect{pp},'_',xRangeVect{xx},yRangeVect{yy}];
                FF = FarField.readCSTtxt([dataPathCSTtxt,fileName]);
                PradCST(ii) = FF.pradInt;
                ii = ii+1;
            end
        end
    end
end
errIntCST = diff(round(PradCST*1e4));
pradIntTestVect = [errIntPhTh,errIntAzEl,errIntCST];
if all(pradIntTestVect < tol)
    disp('Pass: pradInt')
    pradIntPass = true;
    clear FF1 FF2 FF3
else
    disp('FAIL: pradInt')
    pradIntPass = false;
end

% setPower
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_phth_spherical_pos180']);
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
showGridPlots = false;

gridLocal = {'PhTh','AzEl','ElAz'};
gridProj = {'DirCos','TrueView','ArcSin'}; 
gridAstro = {'Horiz','RAdec','GalLongLat'};
gridTest = [gridLocal,gridProj,gridAstro];


% Test PhTh input 
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_phth_spherical_sym180']);
for ii = 1:length(gridTest)
    FF2 = FF1.changeGrid(gridTest{ii});
    if showGridPlots && 1
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2PhTh;
    gridTransPhThTestVect(ii) = isGridEqual(FF3,FF1);
end
clear FF1 FF2 FF3

% Test AzEl input
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_azel_spherical_sym180']);
for ii = 1:length(gridTest)
    FF2 = FF1.changeGrid(gridTest{ii});
    if showGridPlots  && 0
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2AzEl;
    gridTransAzElTestVect(ii) = isGridEqual(FF3,FF1);
end
clear FF1 FF2 FF3

% Test ElAz input
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_elaz_spherical_sym180']);
for ii = 1:length(gridTest)
    FF2 = FF1.changeGrid(gridTest{ii});
    if showGridPlots && 0
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2ElAz;
    gridTransElAzTestVect(ii) = isGridEqual(FF3,FF1);
end
clear FF1 FF2 FF3

% Test DirCos input
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_uv_spherical']);
for ii = 1:length(gridTest)
    FF2 = FF1.changeGrid(gridTest{ii});
    if showGridPlots && 0
        figure, FF2.plot('plotType','2D','showGrid',1,'output','E1','outputType','real')
    end
    FF3 = FF2.grid2DirCos;
    gridTransDirCosTestVect(ii) = isGridEqual(FF3,FF1);
end
clear FF1 FF2 FF3

% Test TrueView input
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_trueview_spherical']);
for ii = 1:length(gridTest)
    FF2 = FF1.changeGrid(gridTest{ii});
    if showGridPlots && 1
%         figure, FF2.plot('plotType','2D','showGrid',1,'output','E1','outputType','real')
        figure, FF2.plot('plotType','2D','showGrid',1)
    end
    FF3 = FF2.grid2TrueView;
    gridTransTrueViewTestVect(ii) = isGridEqual(FF3,FF1);
end
clear FF1 FF2 FF3

gridTransTestVect = [gridTransPhThTestVect,gridTransAzElTestVect,gridTransElAzTestVect,gridTransDirCosTestVect,gridTransTrueViewTestVect];

if all(gridTransTestVect)
    disp('Pass: grid2*')
    gridTransPass = true;
else
    disp('FAIL: grid2*')
    gridTransPass = false;
end

%% Coordinate transformations
showCoorPlots = false;
output = 'E1';
outputType = 'real';
scaleMag = 'lin';
setStdGrd = true;
coorVect = {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3','power'};
% coorVect = {'spherical'};
tol = 1e-10;

% Test spherical input
FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_spherical_lin_sym180']);
for ii = 1:length(coorVect)
    FF2 = FF1.changeCoor(coorVect{ii},setStdGrd);
    if showCoorPlots && 1
        figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
    end
    % Don't reset for power coordinate - just run through
    if ~strcmp(coorVect(ii),'power')
        FF3 = FF2.coor2spherical;
        FFd = FF3 - FF1;
        normE = FFd.norm;
        coorTransSphTestVect(ii) = max(normE) < tol;
    end
end
clear FF1 FF2 FF3

% Test the GRASP spherical input coor on different grids
inputGridNameLocal = {'phth','azel','elaz'};
inputGridNameProj = {'uv','trueview'};
inputGridName = [inputGridNameLocal,inputGridNameProj];
jj = 1;
for gg = 1:length(inputGridName)
    if any(strcmp(inputGridName(gg),inputGridNameLocal))
        FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_',inputGridName{gg},'_spherical_sym180']);
    else
        FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FF_',inputGridName{gg},'_spherical']);
    end
    for ii = 1:length(coorVect)
        FF2 = FF1.changeCoor(coorVect{ii},setStdGrd);
        if showCoorPlots && 1
            figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
        end
        % Don't reset for power coordinate - just run through
        if ~strcmp(coorVect(ii),'power')
            FF3 = FF2.coor2spherical;
            FFd = FF3 - FF1;
            normE = FFd.norm;
            coorTransSphGRASPTestVect(jj) = max(normE) < tol;
            jj = jj+1;
        end
    end
    clear FF1 FF2 FF3
end

% Test AzEl input
FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_azel_lin_sym180']);
for ii = 1:length(coorVect)
    FF2 = FF1.changeCoor(coorVect{ii},setStdGrd);
    if showCoorPlots && 1
        figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
    end
    % Don't reset for power coordinate - just run through
    if ~strcmp(coorVect(ii),'power')
        FF3 = FF2.coor2Ludwig2AE;
        FFd = FF3 - FF1;
        normE = FFd.norm;
        coorTransAzElTestVect(ii) = max(normE) < tol;
    end
end
clear FF1 FF2 FF3

% Test ElAz input
FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_elaz_lin_sym180']);
for ii = 1:length(coorVect)
    FF2 = FF1.changeCoor(coorVect{ii},setStdGrd);
    if showCoorPlots && 1
        figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
    end
    % Don't reset for power coordinate - just run through
    if ~strcmp(coorVect(ii),'power')
        FF3 = FF2.coor2Ludwig2EA;
        FFd = FF3 - FF1;
        normE = FFd.norm;
        coorTransElAzTestVect(ii) = max(normE) < tol;
    end
end
clear FF1 FF2 FF3

% Test Ludwig3 input
FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_L3_lin_sym180']);
for ii = 1:length(coorVect)
    FF2 = FF1.changeCoor(coorVect{ii},setStdGrd);
    if showCoorPlots && 1
        figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
    end
    % Don't reset for power coordinate - just run through
    if ~strcmp(coorVect(ii),'power')
        FF3 = FF2.coor2Ludwig3;
        FFd = FF3 - FF1;
        normE = FFd.norm;
        coorTransL3TestVect(ii) = max(normE) < tol;
    end
end
clear FF1 FF2 FF3

coorTransTestVect = [coorTransSphTestVect,coorTransSphGRASPTestVect,coorTransAzElTestVect,coorTransElAzTestVect,coorTransL3TestVect];

if all(coorTransTestVect)
    disp('Pass: coor2*')
    coorTransPass = true;
else
    disp('FAIL: coor2*')
    coorTransPass = false;
end

%% Polarisation transformations
% Run through a big set from the CST simulations.
showPolPlots = true;
output = 'E1';
outputType = 'phase';
scaleMag = 'lin';
gridVect = {'spherical','azel','elaz','L3'};
polVectFileName = {'lin','circ'};
polVect = {'linear','circular','slant'};
tol = 1e-10;
cc = 1;
for gg = 1:1%length(gridVect)
    for pp = 1:length(polVectFileName)
        FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_',gridVect{gg},'_',polVectFileName{pp},'_sym180']);
        for ii = 1:length(polVect)
            FF2 = FF1.changePol(polVect{ii});
            if showPolPlots && 1
                figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
            end
            FF3 = FF2.changePol(polVect{pp});   % Change back to the type that we read - names are different
            FFd = FF3 - FF1;
            normE = FFd.norm;
            polTransTestVect(cc) = max(normE) < tol;
            cc = cc + 1;
        end
    end
end
clear FF1 FF2 FF3

if all(polTransTestVect)
    disp('Pass: pol2*')
    polTransPass = true;
else
    disp('FAIL: pol2*')
    polTransPass = false;
end

%%
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
FarFieldPass = all([constructorPass,readGRASPgrdPass,readGRASPcutPass,readCSTffsPass,readCSTtxtPass,...
    pradIntPass,setPowerPass,...
    gridTransPass,coorTransPass,polTransPass...
    getRangePass,setRangePass,...
    ]);
if FarFieldPass
    disp('Pass: FarField');
else
    disp('FAIL: FarField');
end

disp('-------------------------------------------------------------------')


