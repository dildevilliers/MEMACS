% Test script for the FarField class
% Created: 2019-05-08, Dirk de Villiers
% Updated: 2019-05-24, Dirk de Villiers

close all
clearvars
clear all

disp('-------------------------------------------------------------------')
disp('...Testing FarField...');

p = mfilename('fullpath');
[filepath] = fileparts(p);
dataPath = [filepath,'\..\data\SimPatterns\'];
dataPathGRASPgrd = [dataPath,'GRASPgrd\'];
dataPathGRASPcut = [dataPath,'GRASPcut\'];
dataPathCSTffs = [dataPath,'CSTffs\'];
dataPathCSTtxt = [dataPath,'CSTtxt\'];
dataPathNFSscan = [dataPath,'NFSscan\'];
dataPathFITS = [dataPath,'FITS\'];


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

try
    inStruct = struct('coorType','Ludwig1','gridType','AzEl','polType','circular',...
                'symmetryXZ','electric','symmetryYZ','none','symmetryXY','none','symmetryBOR','none',...
                'r',2,'slant',deg2rad(30),'freqUnit','MHz',...
                'orientation',[0,0,0],'earthLocation',[0 0 0],'time',FF2.time);
    FF3 = FarField(FF.x,FF.y,FF.E1,FF.E2,FF.freq,1,1,inStruct);
    errConstrStruct = ~isequal(FF3.auxParamStruct,inStruct);
catch constructorStruct_errInfo
    errConstrStruct = true;
end

constructTestVect = [errConstrEmpty,errConstrBasic,errConstrPower,errConstrStruct];
if all(constructTestVect < tol)
    disp('Pass: constructor')
    constructorPass = true;
    clear FF FF1 FF2 FF3
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

% NFS constructor
readNFSscanPass = true;
NFSnameVect = {'az180el360','az360el180','th180ph360','th360ph180'};
for nn = 1:length(NFSnameVect)
    fileName = ['NFSscan_',NFSnameVect{nn}];
    try
        FF = FarField.readNFSscan([dataPathNFSscan,fileName]);
        if 0
            figure, FF.plot('plotType','2D','showGrid',1,'output','E1','outputType','real')
        end
        readNFSscanPass = readNFSscanPass && 1;
    catch readNFSscan_errInfo
        readNFSscanPass = readNFSscanPass && 0;
        keyboard
    end
end
if readNFSscanPass
    disp('Pass: readNFSscan')
    clear FF 
else
    disp('FAIL: readNFSscan')
end

% FITS constructor
readFITSpass = true;
% Seperate real and imag and components - Rick Perley Holography version...
fitsPath = [dataPathFITS,'\ant5'];
inputStructL.pathName1 = {[fitsPath,'LLreal_small.fits'],[fitsPath,'LRreal_small.fits']};
inputStructL.pathName2 = {[fitsPath,'LLimag_small.fits'],[fitsPath,'LRimag_small.fits']};
inputStructL.type1 = 'real';
inputStructL.type2 = 'imag';
[inputStructL.scale1,inputStructL.scale2] = deal('lin');
inputStructL.scaleFuncGrid = @sind;
% Seperate Real and Imaginary - Khan Asad EM models version
fitsPath = [dataPathFITS,'\MeerKAT_EM_'];
inputStructEM.pathName1 = [fitsPath,'real_small.fits'];
inputStructEM.pathName2 = [fitsPath,'imag_small.fits'];
inputStructEM.type1 = 'real';
inputStructEM.type2 = 'imag';
[inputStructEM.scale1,inputStructEM.scale2] = deal('lin');
inputStructEM.scaleFuncGrid = @sind;
% All in one file - Khan Asad MeerKAT Holography version
fitsPath = [dataPathFITS,'\MeerKAT_Holo_'];
inputStructH.pathName1 = [fitsPath,'small.fits'];
inputStructH.scaleFuncGrid = @sind;
try
    FF = FarField.readFITS(inputStructL,'DirCos','Ludwig3','circular');
    FF1 = FarField.readFITS(inputStructEM,'DirCos','Ludwig3','linear','xRange',[-3,3],'yRange',[-3,3],'fRange',[950e6,1670e6]);
    FF2 = FarField.readFITS(inputStructH,'DirCos','Ludwig3','linear','xRange',[-3,3],'yRange',[-3,3],'fRange',[899e6,899e6]);
    if 0
        figure, FF.plot('plotType','2D','showGrid',1)
        figure, plotJones(FF1(1),FF1(2))
        figure, plotJones(FF2(1),FF2(2))
    end
    readFITSpass = readFITSpass && 1;
catch readFITS_sepComp_errInfo
    readFITSpass = readFITSpass && 0;
    keyboard
end
if readFITSpass
    disp('Pass: readFITS')
    clear FF FF1 FF2
else
    disp('FAIL: readFITS')
end

% Struct constructor
FF = FarField;
FFpattern = FF.getFarFieldStruct;
FF1 = FarField.fromStruct(FFpattern);
FarFieldfromStructPass = isequal(FF,FF1);
if FarFieldfromStructPass
    disp('Pass: fromStruct')
    clear FF FF1
else
    disp('FAIL: fromStruct')
end

%% Performance metrics
% Beamwidth
FF1 = FarField.readGRASPgrd([dataPathGRASPgrd,'FFmainBeam']);
% FF1.plot('plotType','2D','showGrid',1,'norm',true,'freqIndex',3), axis normal
[bw,bwStruct] = FF1.getBeamwidth([-10,-inf]);
bwPass = any(~isnan(bwStruct(:)));
if bwPass
    disp('Pass: getBeamwidth')
else
    disp('FAIL: getBeamwidth')
end

% SLL
[sll1,sll2,sllStruct] = FF1.getSLL;
sllPass = any(~isnan([sll1(:);sll2(:)]));
if sllPass
    disp('Pass: getSLL')
else
    disp('FAIL: getSLL')
end

clear FF1

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
    if showGridPlots && 0
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
    if showGridPlots && 0
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
clear FF1 FF2 FF3 FFd

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
    clear FF1 FF2 FF3 FFd
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
clear FF1 FF2 FF3 FFd

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
clear FF1 FF2 FF3 FFd

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
clear FF1 FF2 FF3 FFd

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
showPolPlots = false;
output = 'E1';
outputType = 'phase';
scaleMag = 'lin';
gridVect = {'spherical','azel','elaz','L3'};
polVectFileName = {'lin','circ'};
polVect = {'linear','circular','slant'};
tol = 1e-10;
cc = 1;
for gg = 1:length(gridVect)
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
clear FF1 FF2 FF3 FFd

if all(polTransTestVect)
    disp('Pass: pol2*')
    polTransPass = true;
else
    disp('FAIL: pol2*')
    polTransPass = false;
end

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

% setRangeSph
showGridShiftPlots = true;
output = 'E1';
outputType = 'phase';
scaleMag = 'lin';

tol = 1e-6;
xGrids = {'pos','sym'};
yGrids = {'180','360'};
% PhTh grid
cc = 1;
for ii = 1:2
    xGridIn = xGrids{ii};
    for jj = 1:2
        yGridIn = yGrids{jj};
        FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_spherical_lin_',xGridIn,yGridIn]);
        for mm = 1:2
            xGridTrans = xGrids{mm};
            for nn = 1:2
                yGridTrans = yGrids{nn};
                FF2 = FF1.setRangeSph(xGridTrans,yGridTrans);
                FF3 = FF2.setRangeSph(xGridIn,yGridIn);
                FFd = FF3 - FF1;
                if showGridShiftPlots
%                     figure, FF1.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
                    figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
                    axis normal
%                     figure, FF3.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
%                     figure, FFd.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
                end
                normE = FFd.norm;
                gridChangePhThTestVect(cc) = max(normE) < tol;
                cc = cc + 1;
            end
        end
    end
end
clear FF1 FF2 FF3 FFd

% AzEl grid
cc = 1;
for ii = 1:2
    xGridIn = xGrids{ii};
    for jj = 1:2
        yGridIn = yGrids{jj};
        FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_azel_lin_',xGridIn,yGridIn]);
        for mm = 1:2
            xGridTrans = xGrids{mm};
            for nn = 1:2
                yGridTrans = yGrids{nn};
                FF2 = FF1.setRangeSph(xGridTrans,yGridTrans);
                FF3 = FF2.setRangeSph(xGridIn,yGridIn);
                FFd = FF3 - FF1;
                if showGridShiftPlots
%                     figure, FF1.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
                    figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
                    axis normal
%                     figure, FF3.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
%                     figure, FFd.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
                end
                normE = FFd.norm;
                gridChangeAzElTestVect(cc) = max(normE) < tol;
                cc = cc + 1;
            end
        end
    end
end
clear FF1 FF2 FF3 FFd

% ElAz grid
cc = 1;
for ii = 1:2
    xGridIn = xGrids{ii};
    for jj = 1:2
        yGridIn = yGrids{jj};
        FF1 = FarField.readCSTtxt([dataPathCSTtxt,'FF_elaz_lin_',xGridIn,yGridIn]);
        for mm = 1:2
            xGridTrans = xGrids{mm};
            for nn = 1:2
                yGridTrans = yGrids{nn};
                FF2 = FF1.setRangeSph(xGridTrans,yGridTrans);
                FF3 = FF2.setRangeSph(xGridIn,yGridIn);
                FFd = FF3 - FF1;
                if showGridShiftPlots
%                     figure, FF1.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
                    figure, FF2.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
                    axis normal
%                     figure, FF3.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
%                     figure, FFd.plot('plotType','2D','showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
%                     axis normal
                end
                normE = FFd.norm;
                gridChangeelAzTestVect(cc) = max(normE) < tol;
                cc = cc + 1;
            end
        end
    end
end
clear FF1 FF2 FF3 FFd

if all([gridChangePhThTestVect,gridChangeAzElTestVect,gridChangeelAzTestVect])
    disp('Pass: setRange')
    setRangePass = true;
else
    disp('FAIL: setRange')
    setRangePass = false;
end




%% Static methods
FF = FarField.SimpleTaper(55, -12, -12, 1e9);



%% Final test
FarFieldPass = all([constructorPass,readGRASPgrdPass,readGRASPcutPass,readCSTffsPass,readCSTtxtPass,readNFSscanPass,readFITSpass,FarFieldfromStructPass,...
    bwPass,sllPass,...
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


