%Name: testScript_FarField.m
%Description:
%   Script to test the Farfield object. A patch antenna Farfield is read in from a FEKO ffe file
%   and the use cases of its various methods are shown.
clear all
close all
FF = FarField.readFEKOffe([pwd, '\patch']);

gridType = {'PhTh' 'DirCos' 'AzEl' 'ElAz' 'TrueView' 'ArcSin'}; 
coorType = {'spherical' 'Ludwig1' 'Ludwig2AE' 'Ludwig2EA' 'Ludwig3'};
polType = {'linear' 'circular' 'slant'};
xRangeType = 'pos';

% for ii = 1:6
%     handleGridType = str2func(['grid2',gridType{ii}]);
%     FF = handleGridType(FF);
%     figure
%     plotGrid(FF);
% end
a = 0;
for ii = 1:6
    k = ii;
    if ii > 5
        k = 1;
    end
    for jj = 1:3
        a = a+1;
        if  a == 7
%             keyboard;
        end
        handleGridType = str2func(['grid2',gridType{ii}]);
        handleCoorType = str2func(['coor2',coorType{k}]);
        handlePolType = str2func(['pol2',polType{jj}]);
        FF = handleGridType(FF);
        FF = handleCoorType(FF,0);
        FF = handlePolType(FF);
        FF = FF.setXrange(xRangeType);
        figure
        FF.plot('plotType','2D','step',1,'showGrid',true,'output','E1','outputType','mag');
        figure
        FF.plot('plotType','2D','step',1,'showGrid',true,'output','E2','outputType','mag');
%         figure
%         FF.plot('plotType','2D','step',1,'showGrid',true,'output','CO_XP','outputType','mag');
%         figure
%         FF.plot('plotType','2D','step',1,'showGrid',true,'output','XP_CO','outputType','mag');
    end  
end

figure
plotType = 'polar'; % '3D', '2D', 'polar', 'cartesian'
FF.plot('plotType', plotType);

figure
FF.plotPrincipleCuts;

%% Maths
FF = FarField.readFEKOffe([pwd, '\patch']);

gridType = 'AzEl';
coorType = 'Ludwig3';
polType = 'circular';
xRangeType = 'pos';

handleGridType = str2func(['grid2',gridType]);
handleCoorType = str2func(['coor2',coorType]);
handlePolType = str2func(['pol2',polType]);
FF = handleGridType(FF);
FF = handleCoorType(FF,0);
FF = handlePolType(FF);
FF = FF.setXrange(xRangeType);

%MAKE ANOTHER FARFIELD THATS THE SAME. WE ARE DOING MATHS ON BOTH.
format compact
FF2 = FarField.readFEKOffe([pwd, '\patch']);
FF2 = handleGridType(FF2);
FF2 = handleCoorType(FF2,0);
FF2 = handlePolType(FF2);
FF2 = FF2.setXrange(xRangeType);

FF_plus = plus(FF,FF2);
FF_plus = handleGridType(FF_plus);
FF_plus = handleCoorType(FF_plus,0);
FF_plus = handlePolType(FF_plus);
FF_plus = FF_plus.setXrange(xRangeType);

%FF_plus.E1(1,1)
%FF.E1(1,1)
%FF2.E1(1,1)
if any(FF_plus.E1 ~= (FF.E1 + FF2.E1))
    warning('Plus function failed.')
end
FF_minus = minus(FF,FF2);
if any(FF_minus.E1 ~= (FF.E1 - FF2.E1))
    warning('Minus function failed.')
end
FF_times = times(FF,FF2);
if any(FF_times.E1 ~= (FF.E1.*FF2.E1))
    warning('Times function failed.')
end
FF_conj = conj(FF);
if any(FF_conj.E1 ~= conj(FF.E1))
    warning('Conj function failed.')
end
FF_abs = abs(FF);
if any(FF_abs.E1 ~= abs(FF.E1))
    warning('Abs function failed.')
end
FF_scale = scale(FF,2);
if any(FF_scale.E1 ~= FF.E1.*2)
    warning('Scale function failed.')
end
[normE] = norm(FF);
FF_T = convPower(FF,FF2);

%% Grid getters and setters.
%THE PLAN HERE IS 6 STEPS LONG:
%1) SET TO A SPECIFIC GRID.
%2) GET THAT SPECIFIC GRID.
%3) SET TO A DIFFERENT GRID.
%4) SET TO THE FIRST GRID AGAIN.
%5) GET THAT GRID AGAIN.
%6) COMPARE. IF THEY MATCH -> SUCCESS!

FF = grid2DirCos(FF);
[u, v, w] = getDirCos(FF);
FF = grid2PhTh(FF);
FF = grid2DirCos(FF);
[u2, v2, w2] = getDirCos(FF);
if  (u(1,1) ~= u2(1,1)) || (v(1,1) ~= v2(1,1)) || (w(1,1) ~= w2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Grid transformation failure. DirCos or PhTh')
end

FF = grid2AzEl(FF);
[az, el] = getAzEl(FF);
FF = grid2TrueView(FF);
FF = grid2AzEl(FF);
[az2, el2] = getAzEl(FF);
if  (az(1,1) ~= az2(1,1)) || (el(1,1) ~= el2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Grid transformation failure. AzEl or Trueview')
end

FF = grid2ElAz(FF);
[ep, al] = getElAz(FF);
FF = grid2ArcSin(FF);
FF = grid2ElAz(FF);
[ep2, al2] = getElAz(FF);
if  (ep(1,1) ~= ep2(1,1)) || (al(1,1) ~= al2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Grid transformation failure. ElAz or Arcsin')
end

FF = grid2PhTh(FF);
[ph, th] = getPhTh(FF);
FF = grid2DirCos(FF);
FF = grid2PhTh(FF);
[ph2, th2] = getPhTh(FF);
if  (ph(1,1) ~= ph2(1,1)) || (th(1,1) ~= th2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Grid transformation failure. PhTh or DirCos')
end

FF = grid2TrueView(FF);
[Xg, Yg] = getTrueView(FF);
FF = grid2PhTh(FF);
FF = grid2TrueView(FF);
[Xg2, Yg2] = getTrueView(FF);
if  (Xg(1,1) ~= Xg2(1,1)) || (Yg(1,1) ~= Yg2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Grid transformation failure. TrueView or PhTh')
end

FF = grid2ArcSin(FF);
[asinu, asinv] = getArcSin(FF);
FF = grid2PhTh(FF);
FF = grid2ArcSin(FF);
[asinu2, asinv2] = getArcSin(FF);
if  (asinu(1,1) ~= asinu2(1,1)) || (asinv(1,1) ~= asinv2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Grid transformation failure. ArcSin or PhTh')
end

%% Coordinate system getters and setters.
%THE PLAN HERE IS 6 STEPS LONG:
%1) SET TO A SPECIFIC GRID.
%2) GET THAT SPECIFIC GRID.
%3) SET TO A DIFFERENT GRID.
%4) SET TO THE FIRST GRID AGAIN.
%5) GET THAT GRID AGAIN.
%6) COMPARE. IF THEY MATCH -> SUCCESS!

FF = coor2spherical(FF);
[Eth, Eph, Er] = getEspherical(FF);
FF = coor2Ludwig1(FF);
FF = coor2spherical(FF);
[Eth2, Eph2, Er2] = getEspherical(FF);
if  (Eth(1,1) ~= Eth2(1,1)) || (Eph(1,1) ~= Eph2(1,1)) || (Er(1,1) ~= Er2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Coor transformation failure. Spherical or Ludwig1')
end

FF = coor2Ludwig2AE(FF);
[Eaz, Eel, E3] = getELudwig2AE(FF);
FF = coor2Ludwig2EA(FF);
FF = coor2Ludwig2AE(FF);
[Eaz2, Eel2, E32] = getELudwig2AE(FF);
if  (Eaz(1,1) ~= Eaz2(1,1)) || (Eel(1,1) ~= Eel2(1,1)) || (E3(1,1) ~= E32(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Coor transformation failure. Ludwig2EA or Ludwig2AE')
end

FF = coor2Ludwig3(FF);
[Eh, Ev, E3] = getELudwig3(FF);
FF = coor2Ludwig2EA(FF);
FF = coor2Ludwig3(FF);
[Eh2, Ev2, E32] = getELudwig3(FF);
if  (Eh(1,1) ~= Eh2(1,1)) || (Ev(1,1) ~= Ev2(1,1)) || (E3(1,1) ~= E32(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Coor transformation failure. Ludwig2EA or Ludwig3')
end

FF = coor2Ludwig2EA(FF);
[Eal, Eep, E3] = getELudwig2EA(FF);
FF = coor2Ludwig3(FF);
FF = coor2Ludwig2EA(FF);
[Eal2, Eep2, E32] = getELudwig2EA(FF);
if  (Eal(1,1) ~= Eal2(1,1)) || (Eep(1,1) ~= Eep2(1,1)) || (E3(1,1) ~= E32(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Coor transformation failure. Ludwig3 or Ludwig2EA')
end

FF = coor2Ludwig1(FF);
[Ex, Ey, Ez] = getELudwig1(FF);
FF = coor2spherical(FF);
FF = coor2Ludwig1(FF);
[Ex2, Ey2, Ez2] = getELudwig1(FF);
if  (Ex(1,1) ~= Ex2(1,1)) || (Ey(1,1) ~= Ey2(1,1)) || (Ez(1,1) ~= Ez2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Coor transformation failure. Ludwig1 or Spherical')
end

%% Polarization type getters and setters.
%THE PLAN HERE IS 6 STEPS LONG:
%1) SET TO A SPECIFIC GRID.
%2) GET THAT SPECIFIC GRID.
%3) SET TO A DIFFERENT GRID.
%4) SET TO THE FIRST GRID AGAIN.
%5) GET THAT GRID AGAIN.
%6) COMPARE. IF THEY MATCH -> SUCCESS!

FF = pol2linear(FF);
[E1lin, E2lin, E3lin] = getElin(FF);
FF = pol2circular(FF);
FF = pol2linear(FF);
[E1lin2, E2lin2, E3lin2] = getElin(FF);
if  (E1lin(1,1) ~= E1lin2(1,1)) || (E2lin(1,1) ~= E2lin2(1,1)) || (E3lin(1,1) ~= E3lin2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Pol transformation failure. Linear or Circular')
end

FF = pol2slant(FF);
[Exp,Eco,E3slant] = getEslant(FF);
FF = pol2circular(FF);
FF = pol2slant(FF);
[Exp2,Eco2] = getEslant(FF);
if  (Exp(1,1) ~= Exp2(1,1)) || (Eco(1,1) ~= Eco2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Pol transformation failure. Slant or Circular')
end

FF = pol2circular(FF);
[Elh,Erh,E3circ] = getEcircular(FF);
FF = pol2slant(FF);
FF = pol2circular(FF);
[Elh2,Erh2] = getEcircular(FF);
if  (Elh(1,1) ~= Elh2(1,1)) || (Erh(1,1) ~= Erh2(1,1)) %DOESNT WORK IF YOU DONT SPECIFY ELEMENT.
    warning('Pol transformation failure. Circular or Slant')
end
%===================================================================================================================================

