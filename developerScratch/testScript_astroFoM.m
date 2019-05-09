% load sky map in galactic coordinates

fits = fitsread('pixel_coords_map_nested_galactic_res9.fits','binarytable');
long = wrapToPi(deg2rad(fits{1}));
lat = wrapToPi(deg2rad(fits{2}));
fits = fitsread('lambda_haslam408_dsds.fits','binarytable');
skymap = fits{1}./1e3;%fits{1};
scatint = scatteredInterpolant([long,lat],skymap);

longgrid = linspace(min(long),max(long),3073);
latgrid = linspace(min(lat),max(lat),1025);
[LATGRID,LONGGRID] = ndgrid(latgrid,longgrid);
skymap_grid = scatint(LONGGRID(:),LATGRID(:));

FFsky = FarField.farFieldFromPowerPattern(LONGGRID(:),LATGRID(:),skymap_grid(:),408e6,'linearY','Ludwig3','GalLongLat');

%Load or generate antenna farfield in PhTh coordinates
th = linspace(0,pi,37);
ph = linspace(-pi,pi,73);
[PH, TH] = meshgrid(ph,th);
D = dipoleD(TH(:),PH(:),1e9,(3e8/1e9)/8);
FFd = FarField.farFieldFromPowerPattern(PH(:),TH(:),D(:),408e6,'linearY','Ludwig3','PhTh');

% set time, location and orientation of antenna
% FFd = FFd.setOrientation([0 pi/2]);
% FFd = FFd.setLocation( );
% FFd = FFd.setTime( );

% convert antenna farfield to galactic coordinates
FFdt = fastPhTh2AstroGrid(FFd,'GalLongLat',FFsky.x,FFsky.y);

%Get Tsky by equ. 12 of Kolitsidas et al, at desired frequency

freq = linspace(50e6,250e6,21);
time = 0:3:23;
beta = 2.5;
Tcmb = 0; %TEMPORARY
fsky = 408e6;
tic;
for tt = 1:length(time)
    for ff = 1:length(freq)
        f = freq(ff);
        Tsky_curr = (skymap_grid - Tcmb).*((f/fsky).^(-beta)) + Tcmb;
        FFsky_curr = FarField.farFieldFromPowerPattern(LONGGRID(:),LATGRID(:),Tsky_curr(:),freq(ff),'linearY','Ludwig3','GalLongLat');
        D = dipoleD(TH(:),PH(:),freq(ff),(3e8/max(freq))/3);
        FFd = FarField.farFieldFromPowerPattern(PH(:),TH(:),D(:),freq(ff),'linearY','Ludwig3','PhTh');
        FFd = FFd.setTime(datetime(2018,7,22,time(tt),0,0));
        FFdt = fastPhTh2AstroGrid(FFd,'GalLongLat',FFsky_curr.x,FFsky_curr.y);
        Tantf(ff,tt) = convPower(FFsky_curr,FFdt)/(4*pi);
        disp(['hour: ',num2str(time(tt)),', frequency: ',num2str(freq(ff)./1e6)])
    end
end
toc
% Perform convolution integral to get antenna sky temperature (should be function of frequency and time!)

% Now we can get waterfall diagram

% Calculate sky figure-of-merit

%Fit residual over frequency

% FFd = FarField.farFieldFromPowerPattern(PH(:),TH(:),D(:),freq(ff),'linearY','Ludwig3','PhTh');
% FFdt = fastPhTh2AstroGrid(FFd,'GalLongLat',FFsky_curr.x,FFsky_curr.y);
% FFdt.plot('plotType','2D','outputType','mag','scaleMag','lin','step',1,'showGrid',false,'dynamicRange_dB',30)