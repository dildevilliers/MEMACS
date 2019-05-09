% Test script for the CoordinateSystem class
% Created: 2019-05-09, Dirk de Villiers
% Updated: 2019-05-09, Dirk de Villiers

close all
clearvars

disp('-------------------------------------------------------------------')
disp('...Testing CoordinateSystem...');

%% Constructor - just run them through
try
    C0 = CoordinateSystem;
    C = CoordinateSystem(Pnt3D(1,2,3),[1;0;0],[0;-2;0],[]);
    disp('Pass: constructors')
    constructorPass = true;
catch contructor_errInfo
    disp('FAIL: constructors')
    constructorPass = false;
end

%% Setters
Cb = CoordinateSystem(Pnt3D(1,2,3),[1;0;0],[0;-2;0],[]);
C = CoordinateSystem(Pnt3D(1,-1,0),[1;0;0],[0;1;0],Cb);
C = C.set2Base;
if (isequal(Cb.origin,C.origin) && all(Cb.x_axis == C.x_axis) && all(Cb.y_axis == C.y_axis))
    disp('Pass: set2Base')
    set2BasePass = true;
else
    disp('FAIL: set2Base')
    set2BasePass = false;
end

%% Testers
C0 = CoordinateSystem;
C1 = CoordinateSystem(Pnt3D(1,0,0));
B = isequal(C0,C1,1e-12);
Be = isequal(C0,C0);
if ~B && Be
    disp('Pass: isequal')
    isequalPass = true;
else
    disp('FAIL: isequal')
    isequalPass = false;
end 

%% Translate
C0 = CoordinateSystem;
delta = [1,2,-0.5];
Cd = C0.translate(delta);
if all(Cd.origin.pointMatrix == delta.')
    disp('Pass: translate')
    translatePass = true;
else
    disp('FAIL: translate')
    translatePass = false;
end

%% Rotate
C0 = CoordinateSystem;
Cr = C0.rotX(deg2rad(45));
if all(abs(Cr.z_axis - [0;-sqrt(2);sqrt(2)]./2) < 1e-12)
    disp('Pass: rotX')
    rotXPass = true;
else
    disp('FAIL: rotX')
    rotXPass = false;
end

C0 = CoordinateSystem;
Cr = C0.rotY(deg2rad(45));
if all(abs(Cr.z_axis - [sqrt(2);0;sqrt(2)]./2) < 1e-12)
    disp('Pass: rotY')
    rotYPass = true;
else
    disp('FAIL: rotY')
    rotYPass = false;
end

C0 = CoordinateSystem;
Cr = C0.rotZ(deg2rad(45));
if all(abs(Cr.x_axis - [sqrt(2);sqrt(2);0]./2) < 1e-12)
    disp('Pass: rotZ')
    rotZPass = true;
else
    disp('FAIL: rotZ')
    rotZPass = false;
end

C0 = CoordinateSystem;
Cr = C0.rotGRASP(deg2rad([45,45,0]));
if all(abs(Cr.z_axis - [1;1;sqrt(2)]./2) < 1e-12)
    disp('Pass: rotGRASP')
    rotGRASPPass = true;
else
    disp('FAIL: rotGRASP')
    rotGRASPPass = false;
end

C0 = CoordinateSystem;
Cr = C0.rotEuler(deg2rad([45,45,0]));
if all(abs(Cr.z_axis - [1;1;sqrt(2)]./2) < 1e-12) && all(abs(Cr.x_axis - [1;1;-sqrt(2)]./2) < 1e-12)
    disp('Pass: rotEuler')
    rotEulerPass = true;
else
    disp('FAIL: rotEuler')
    rotEulerPass = false;
end

%% Change of basis
C0 = CoordinateSystem;
Cr = C0.rotEuler(deg2rad([25,-35,157]));
Q = dirCosine(Cr,C0);
if all(abs(Cr.x_axis - Q(:,1)) < 1e-12) && all(abs(Cr.y_axis - Q(:,2)) < 1e-12) && all(abs(Cr.z_axis - Q(:,3)) < 1e-12)
    disp('Pass: dirCosine')
    dirCosinePass = true;
else
    disp('FAIL: dirCosine')
    dirCosinePass = false;
end

C0 = CoordinateSystem;
angRot = [30,60,45];
Cr = C0.rotGRASP(deg2rad(angRot));
angGRASP = rad2deg(getGRASPangBetweenCoors(Cr,C0));
if all(abs(angRot - angGRASP) < 1e-12) 
    disp('Pass: getGRASPangBetweenCoors')
    getGRASPangBetweenCoorsPass = true;
else
    disp('FAIL: getGRASPangBetweenCoors')
    getGRASPangBetweenCoorsPass = false;
end

C0 = CoordinateSystem;
angRot = [30,60,45];
Cr = C0.rotEuler(deg2rad(angRot));
angEuler = rad2deg(getEulerangBetweenCoors(Cr,C0));
if all(abs(angRot - angEuler) < 1e-12) 
    disp('Pass: getEulerangBetweenCoors')
    getEulerangBetweenCoorsPass = true;
else
    disp('FAIL: getEulerangBetweenCoors')
    getEulerangBetweenCoorsPass = false;
end

C0 = CoordinateSystem;
C1 = C0.rotGRASP(deg2rad([20,-35,139]));
C0.base = C1;
C2 = C0.getInGlobal; 
if isequal(C1,C2,1e-12)
    disp('Pass: getInGlobal')
    getInGlobalPass = true;
else
    disp('FAIL: getInGlobal')
    getInGlobalPass = false;
end

C0 = CoordinateSystem;
C1 = C0.rotGRASP(deg2rad([20,-35,139]));
C0.base = C1;
C2 = C0.redefineToOtherBase(CoordinateSystem); 
if all(abs(C2.x_axis - C1.x_axis) < 1e-12) && all(abs(C2.y_axis - C1.y_axis) < 1e-12) && all(abs(C2.z_axis - C1.z_axis) < 1e-12)
    disp('Pass: redefineToOtherBase')
    redefineToOtherBasePass = true;
else
    disp('FAIL: redefineToOtherBase')
    redefineToOtherBasePass = false;
end

%% Plotters - just run them through
try
    C0 = CoordinateSystem;
    C1 = C0.rotGRASP(deg2rad([20,-35,139]));
    C0.base = C1;
    C0.plot(0.5), hold on
    C0.plotLocal(2)
    disp('Pass: plotters')
    plottersPass = true;
catch contructor_errInfo
    disp('FAIL: plotters')
    plottersPass = false;
end

%% Final test
CoordinateSystemPass = all([constructorPass,set2BasePass,translatePass,...
    isequalPass,rotXPass,rotYPass,rotZPass,rotGRASPPass,rotEulerPass,...
    dirCosinePass,getGRASPangBetweenCoorsPass,getEulerangBetweenCoorsPass,...
    getInGlobalPass,redefineToOtherBasePass,plottersPass]);

if CoordinateSystemPass
    disp('Pass: CoordinateSystem');
else
    disp('FAIL: CoordinateSystem');
end

disp('-------------------------------------------------------------------')


