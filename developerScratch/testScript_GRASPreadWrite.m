% Script to test the GRASP file reading and writing functionality

close all
clear all

% Make a simple pattern
FF = FarField;
FF = FF.coor2spherical;
% Write it to GRASP cut file
path = '..\ExampleSims\GRASP\GRASPcutTest\working\';
Name = 'Feed';
pathName = [path,Name];
FF.writeGRASPcut(pathName);

phValDeg = 10;
FF.plot('plotType','cartesian','step',1,'output','E1','cutValue',deg2rad(phValDeg)), hold on
FF.plot('plotType','cartesian','step',1,'output','E2','cutValue',deg2rad(phValDeg))

% Read it back in
Name = 'FFout';
pathName = [path,Name];
Nf = 1;
Ncuts = 73;
FFin = FarField.readGRASPcut(pathName,Nf,Ncuts);
FFin = FFin.setFreq(FF.freq,FF.freqUnit);
FFin.plot('plotType','cartesian','step',1,'output','E1','cutValue',deg2rad(phValDeg),'lineStyle','--','Color','r')
FFin.plot('plotType','cartesian','step',1,'output','E2','cutValue',deg2rad(phValDeg),'lineStyle','--','Color','r')

figure
FFin.plot('plotType','cartesian','step',1,'output','Directivity','cutValue',deg2rad(phValDeg),'lineStyle','--','Color','r')
% Compare