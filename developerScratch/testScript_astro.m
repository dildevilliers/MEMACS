scaleMag = 'lin';

th = linspace(0,pi,37);
ph = linspace(-pi,pi,73);
[PH, TH] = meshgrid(ph,th);
D = dipoleD(TH(:),PH(:),1e9,(3e8/1e9)/8);
FF = FarField.farFieldFromPowerPattern(PH(:),TH(:),D(:),1e9,'linearY');

% FF = FarField();

% FF = FF.setOrientation([pi/6 pi/4]);
% FF = FF.setOrientation([0 pi/2]);
% 
% FF = FF.grid2AzAlt;
% figure
% subplot(3,1,1)
% FF.plot('plotType','2D','outputType','mag','scaleMag',scaleMag,'step',1,'showGrid',true,'dynamicRange_dB',30)
% title('Az-Alt Directivity')
% % FF = FF.grid2Mollweide;
% % FF = FF.grid2AzAlt;
% % figure
% % FF.plot('plotType','2D','outputType','mag','scaleMag',scaleMag,'step',1,'showGrid',true,'dynamicRange_dB',30)
% % % FFa = FF - FFa;
% % % FFa.plot('plotType','2D','outputType','mag','scaleMag','dB','step',1,'showGrid',true,'dynamicRange_dB',30)
% % 
% % % FF = FF.grid2PhTh;
% % % figure
% % % plot(FF)
% % 
% FF1 = FF.grid2RAdec;
% % FF1 = FF1.grid2Mollweide;
% subplot(3,1,2)
% FF1.plot('plotType','2D','outputType','mag','scaleMag',scaleMag,'step',1,'showGrid',true,'dynamicRange_dB',30)
% title('RA-dec Directivity')
% 
% FF2 = FF1.grid2GalLongLat;
% subplot(3,1,3)
% FF2.plot('plotType','2D','outputType','mag','scaleMag',scaleMag,'step',1,'showGrid',true,'dynamicRange_dB',30)
% title('Gal. Lat-Long Directivity')

%% test PhTh-to-AzAlt-to-PhTh grid conversion accuracy
FFa = FF.grid2AzAlt;
FFa = FF.grid2PhTh;
FFa = FF - FFa;
figure
FFa.plot('plotType','2D','outputType','mag','scaleMag','lin','step',1,'showGrid',true,'dynamicRange_dB',30)

%% test AzAlt-to-RAdec-to-AzAlt grid conversion accuracy
FF = FF.grid2AzAlt;
FFb = FF.grid2RAdec;
FFb = FF.grid2AzAlt;
FFb = FF - FFb;
figure
FFb.plot('plotType','2D','outputType','mag','scaleMag','lin','step',1,'showGrid',true,'dynamicRange_dB',30)

%% test AzAlt-to-GalLongLat-to-AzAlt grid conversion accuracy
FF = FF.grid2AzAlt;
FFc = FF.grid2GalLongLat;
FFc = FF.grid2AzAlt;
FFc = FF - FFc;
figure
FFc.plot('plotType','2D','outputType','mag','scaleMag','lin','step',1,'showGrid',true,'dynamicRange_dB',30)

%% test RAdec-to-GalLongLat-to-RAdec grid conversion accuracy
FF = FF.grid2RAdec;
FFd = FF.grid2GalLongLat;
FFd = FF.grid2RAdec;
FFd = FF - FFd;
figure
FFd.plot('plotType','2D','outputType','mag','scaleMag','lin','step',1,'showGrid',true,'dynamicRange_dB',30)