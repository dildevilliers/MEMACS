% Test script for the Pnt3D class
% Created: 2019-05-07, Dirk de Villiers
% Updated: 2019-05-07, Dirk de Villiers

close all
clearvars

disp('-------------------------------------------------------------------')
disp('...Testing Pnt3D...');

%% Constructors - just run them through
try
    N = 91;
    PH = linspace(0,2*pi,N);
    TH = linspace(0,pi/2,N);
    R = linspace(1,2,N);
    pSph = Pnt3D.sph(PH,TH,R);
    
    PH = linspace(0,2*pi,N);
    RHO = linspace(0,4,N);
    Z = linspace(1,2,N);
    p = Pnt3D.sph(PH,RHO,Z);
    
    p = Pnt3D;
    
    disp('Pass: constructors')
    constructorPass = true;
catch contructor_errInfo
    disp('FAIL: constructors')
    constructorPass = false;
end

%% Setters - just run them through
p = Pnt3D;
try
    p = p.setX([1,2]);
    p = p.setY([2,1]);
    p = p.setZ([-0.5,3]);
    
    disp('Pass: Setters')
    setterPass = true;
catch setter_errInfo
    disp('FAIL: Setters')
    setterPass = false;
end

%% Object operations
[x,y,z] = deal(1:5);
p0 = Pnt3D(x,y,z);

% getNpts
iN = [1,5];
p1 = p0.getNpts(iN);
if (all(p1.x == x(iN)) && all(p1.y == y(iN)) && all(p1.z == z(iN)))
    disp('Pass: getNpts')
    getNptsPass = true;
else
    disp('FAIL: getNpts')
    getNptsPass = false;
end

% pointMatrix
if (all(p0.x == x(:).') && all(p0.y == y(:).') && all(p0.z == z(:).'))
    disp('Pass: pointMatrix')
    pointMatrixPass = true;
else
    disp('FAIL: pointMatrix')
    pointMatrixPass = false;
end


%% Overloaded methods
p0 = Pnt3D;
[x,y,z] = deal(1:5);
p1 = Pnt3D(x,y,z);

% plus
p2 = p0+p1;
if (all(p2.x == x) && all(p2.y == y) && all(p2.z == z))
    disp('Pass: plus')
    plusPass = true;
else
    disp('FAIL: plus')
    plusPass = false;
end

% minus
p3 = p2-p1;
if (all(p3.x == 0) && all(p3.y == 0) && all(p3.z == 0))
    disp('Pass: minus')
    minusPass = true;
else
    disp('FAIL: minus')
    minusPass = false;
end

% size
if all(size(p3) == size(x))
    disp('Pass: size')
    sizePass = true;
else
    disp('FAIL: size')
    sizePass = false;
end 

% isequal
p4 = Pnt3D(zeros(1,5),zeros(1,5),zeros(1,5));
if isequal(p0,Pnt3D) && isequal(p4,p3) && ~isequal(p3,p2) && ~isequal(p4,p0)
    disp('Pass: isequal')
    isequalPass = true;
else
    disp('FAIL: isequal')
    isequalPass = false;
end 

%% Object element operations
[x,y,z] = deal(1:5);
p1 = Pnt3D(x,y,z);

% scale
scaleVal = [3,3,0.5].';
pS = p1.scale(scaleVal);
if all(pS.x == p1.x*scaleVal(1)) && all(pS.y == p1.y*scaleVal(2)) && all(pS.z == p1.z*scaleVal(3))
    disp('Pass: scale')
    scalePass = true;
else
    disp('FAIL: scale')
    scalePass = false;
end

% distanceCart
p0 = Pnt3D(1,1,1);
xyz = 1:5;
[x,y,z] = deal(xyz);
p1 = Pnt3D(x,y,z);
pd = distanceCart(p0,p1);
d = sqrt(sum(([1,1,1].' - [xyz;xyz;xyz]).^2,1));
if all(abs(pd - d) < 1e-15)
    disp('Pass: distanceCart')
    distanceCartPass = true;
else
    disp('FAIL: distanceCart')
    distanceCartPass = false;
end

% addVect
[x,y,z] = deal(xyz);
p0 = Pnt3D(x,y,z);
DEL = [-1,-1,-1].';
p1 = p0.addVect(DEL);
if all(p1.x == p0.x+DEL(1)) && all(p1.y == p0.y+DEL(2)) && all(p1.z == p0.z+DEL(3))
    disp('Pass: addVect')
    addVectPass = true;
else
    disp('FAIL: addVect')
    addVectPass = false;
end

% changeBase
p0 = Pnt3D(1,1,1);
c0 = CoordinateSystem;
c1 = c0.translate([0,0,1]);
p1 = p0.changeBase(c0,c1);
[x1,y1,z1] = deal(p0.x,p0.y,p0.z+1);
p0 = Pnt3D(0,0,1);
c2 = c0.rotGRASP(deg2rad([45,45,0]));
p2 = p0.changeBase(c0,c2);
[x2,y2,z2] = deal(0.5,0.5,cosd(45));
if all(p1.pointMatrix == [x1,y1,z1].') && all(abs(p2.pointMatrix - [x2,y2,z2].') < 1e-15)
    disp('Pass: changeBase')
    changeBasePass = true;
else
    disp('FAIL: changeBase')
    changeBasePass = false;
end

%% Plotting - just run through
try
    % plotVect
    [x,y] = deal(1:5);
    z = 0;
    p = Pnt3D(x,y,z);
    V = [0.2,2,0.7].';
    figure
    p.plotVect(V,'lineStyle','-','lineWidth',2), hold on
    p.plot('marker','*')
    
    p1 = Pnt3D;
    V = [1:5;zeros(1,5);ones(1,5)];
    figure
    p1.plotVect(V,'lineStyle','-','lineWidth',2), hold on
    p1.plot('marker','*')
    
    % plotLines
    p1 = Pnt3D;
    x = -4:2:4;
    [y,z] = deal(1:5);
    y = y.^2 - 2;
    p2 = Pnt3D(x,y,z);
    figure
    plotLines(p1,p2,'lineStyle','--','lineWidth',1.5,'lineColor',[0.5,0.5,0.5]), hold on
    p1.plot('marker','*')
    p2.plot('marker','o')
    
    % plot
    p = Pnt3D;
    figure
    p.plot
    
    [x,y] = deal(1:5);
    z = 0;
    p = Pnt3D(x,y,z);
    figure
    p.plot('marker','o','lineStyle','-','lineWidth',2)
    
    disp('Pass: plotters')
    plotterPass = true;
catch plotter_errInfo
    disp('FAIL: plotters')
    plotterPass = false;
end


%% Final test
Pnt3DPass = all([constructorPass,setterPass,...
    getNptsPass,pointMatrixPass,plusPass,minusPass,sizePass,...
    isequalPass,scalePass,distanceCartPass,addVectPass,changeBasePass,plotterPass]);
if Pnt3DPass
    disp('Pass: Pnt3D');
else
    disp('FAIL: Pnt3D');
end

disp('-------------------------------------------------------------------')
