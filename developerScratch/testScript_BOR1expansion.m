% Test the BOR1 expansion functionality of the FarField class
close all
clear all

pathName = 'Farfield Source [1]';
FF = FarField.readCSTffs(pathName);
FFBOR1 = FF.getBORpattern(1);
FFexpand = FFBOR1.expandBORpattern;
FFstruct = FF.getFarFieldStruct;
BOR1struct = BOR1(FFstruct);

% Plot it...
figure
FFBOR1.plot('plotType','cartesian','step',1,'cutValue',0), hold on
plot(rad2deg(unique(FF.th)),dB10(BOR1struct.C1(:,1).^2./sqrt(2*FF.eta0)),'r--')
FFBOR1.plot('plotType','cartesian','step',1,'cutValue',pi/2), hold on
plot(rad2deg(unique(FF.th)),dB10(BOR1struct.A1(:,1).^2./sqrt(2*FF.eta0)),'b--')
xlabel('\theta^\circ')

scaleMag = 'lin';
outputType = 'real';
output = 'E2';
plotGridHandle = @grid2PhTh;
plotCoorHandle = @coor2Ludwig3;
FFplot = plotCoorHandle(FF);
FFplot = plotGridHandle(FFplot);
FFexpandPlot = plotCoorHandle(FFexpand);
FFexpandPlot = plotGridHandle(FFexpandPlot);


figure
FFplot.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
figure
FFexpandPlot.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)


FFdelta = FFexpand - FF;
errorE = FFdelta.norm
figure
FFdelta.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType,'scaleMag',scaleMag)
