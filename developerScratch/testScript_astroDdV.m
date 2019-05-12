close all
clearvars

load('c:\Work\GitHub\MEMACS\data\HaslamDNDSdata','longU','latU','skymapU')

freq = 408e6;
FF = FarField.farFieldFromPowerPattern(longU(:),latU(:),skymapU(:),freq,'gridType','GalLongLat');

FF = FF.setTime(datetime('now') + hours(10));
FF = FF.grid2PhTh;
FF.plot('plotType','2D','step',1,'showGrid',false)

