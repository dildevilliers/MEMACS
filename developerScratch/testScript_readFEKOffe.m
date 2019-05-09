%Name: testScript_readFEKOffe.m
%Description:
%   Script to test Farfield object creation from FEKO .ffe file. Creates
%   the object with a .ffe generated from a patch antenna simulated in
%   FEKO.

%create farfield object
FF = FarField.readFEKOffe([pwd, '\patch']);

%2D plot of mag(E1), 1 degree step, 40dB dynamic range and visible grid 
figure
FF.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);
%2D plot of mag(E2), 1 degree step, 40dB dynamic range and visible grid 
figure
FF.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30);

%1D Cartesian plot of directivities
F = openfig('patch_FEKO_directivities.fig');
figure(F)
hold on
FF.plot('plotType','cartesian','output','Directivity','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30,'cutConstant','x','cutValue',0);
FF.plot('plotType','cartesian','output','Directivity','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30,'cutConstant','x','cutValue',deg2rad(40),'Color','r');
FF.plot('plotType','cartesian','output','Directivity','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',30,'cutConstant','x','cutValue',deg2rad(90),'Color','b');
legend('\theta = 0^{\circ}','\theta = 45^{\circ}','\theta = 90^{\circ}')
