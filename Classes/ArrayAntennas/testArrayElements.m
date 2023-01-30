close all
clear all

freq = 1.57542e9;
c = physconst('LightSpeed'); %propagation velocity
lambda = c/freq;
r = 0.5*lambda; %distance between elements

%Number of antenna elements
Nant = 12;

%number of sensors elements dense array

xPos = (-(Nant-1)/2:(Nant-1)/2).*r; %Position of antennas
antPos = Pnt3D(xPos,0,0);

array = ArrayElements(antPos);

th = deg2rad(90);
ph = deg2rad(linspace(-180,180,720));
th = repmat(th,1,length(ph));
array.plotAF(freq,th,ph)