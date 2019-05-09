% Name: Test readCSTffs
%
% Description: 
%   Script to test the FarField objective creation from a CST .ffs file
%   Tests;  FF obj creation
%           Plot: - 1D - Directivity
%                 - 1D - CO/XP
%                 - 2D - mag of E1
%                 - 2D - mag of E2


% Create FarField obj
pathName = 'CircWG_origin';

FF = FarField.readCSTffs(pathName);

% 1D Plots - plotPrincipleCuts
plotPrincipleCuts(FF,'dynamicRange_dB',40,'plotType','cartesian','output','Directivity');

plotPrincipleCuts(FF,'dynamicRange_dB',40,'plotType','cartesian','output','CO_XP');
 
% 2D Plots - plot
figure('Name','2D - E1')
FF.plot('plotType','2D','output','E1','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);

figure('Name','2D - E2')
FF.plot('plotType','2D','output','E2','outputType','mag','step',1,'showGrid',true,'dynamicRange_dB',40);

