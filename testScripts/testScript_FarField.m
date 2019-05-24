% Test script for the FarField class
% Created: 2019-05-08, Dirk de Villiers
% Updated: 2019-05-24, Dirk de Villiers

close all
clearvars

disp('-------------------------------------------------------------------')
disp('...Testing FarField...');

p = mfilename('fullpath');
[filepath] = fileparts(p);

%% Constructors 

% ReadNFSscan
dataPath = [filepath,'\..\data\'];
output = 'E1';
outputType = 'imag';
FFth180ph360 = FarField.readNFSscan([dataPath,'NFSscan_th180ph360']);
figure, FFth180ph360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
FFth360ph180 = FarField.readNFSscan([dataPath,'NFSscan_th360ph180']);
figure, FFth360ph180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
FFaz360el180 = FarField.readNFSscan([dataPath,'NFSscan_az360el180']);
figure, FFaz360el180.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)
FFaz180el360 = FarField.readNFSscan([dataPath,'NFSscan_az180el360']);
figure, FFaz180el360.plot('plotType','2D','step',1,'showGrid',1,'output',output,'outputType',outputType)

%% Final test


disp('-------------------------------------------------------------------')


