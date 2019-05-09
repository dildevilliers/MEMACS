classdef FarField
    %FARFIELD  Class of radiated FarField patterns (2 field orthogonal field components) 
    % Based largely on information in CST help file: Farfield Calculation.
    % Overview, and the AMTA paper in this folder 'COORDINATE SYSTEM
    % PLOTTING FOR ANTENNA MEASUREMENTS', GF Masters and SF Gregson.
    % Provides a wide variety of methods for manupulating farfield patterns
    % including plotting, comparison, transformations between grids and coordinate
    % system and polarisation types, etc.
    %
    % Constructor methods
    % - FarField
    % - readGRASPgrd
    % - readFEKOffe
    % - readCSTffs
    % - readGRASPcut
    % - farFieldFromPowerPattern
    
    properties
        r(1,1) double {mustBeReal, mustBeFinite}        % Radius where E-field is evaluated in (m)
        radEff(1,:) double {mustBeReal, mustBeFinite}   % Radiation efficiency per unit
        slant(1,1) double {mustBeReal, mustBeFinite}    % slant angle (for polType=slant) in radians
    end
    
    properties (SetAccess = private)
        x           % First (azimuth) grid parameter
        y           % Second (polar) grid parameter
        E1          % First E-field component
        E2          % Second E-field component
        freq        % Frequency
        Prad        % Radiated power per frequency
        coorType    % Coordinate system type {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3'}
        polType     % polarization type {'linear','circular','slant'}
        gridType    % Grid type {'PhTh','DirCos','AzEl','ElAz','AzAlt','TrueView','ArcSin','Mollweide','RAdec','GalLongLat'}
        freqUnit    % Frequency Unit {'Hz','kHz','MHz','GHz','THz'}
        symmetryXZ  % XZ plane symmetry type {'none','electric','magnetic'}
        symmetryYZ  % YZ plane symmetry type {'none','electric','magnetic'}
        symmetryXY  % XY plane symmetry type {'none','electric','magnetic'}
        symmetryBOR % BOR symmetry type {'none','BOR0','BOR1'}
    end
    
    properties (Dependent = true)
        ph          % Spherical coordinate phi angle of grid
        th          % Spherical coordinate theta angle of grid
        xname       % Name of the x-grid variable {'\phi','u','az','\epsilon','Xg=asin(u)','Xg','North-az','RA','long'}
        yname       % Name of the y-grid variable {'\theta','v','el','\alpha','Yg=asin(v)','Yg','Alt','dec','lat'}
        E1name      % Name of the E1-field component {'Eth','Ex','Eaz','Eal','Eh','Exp','Elh'}
        E2name      % Name of the E2-field component {'Eph','Ey','Eel','Eep','Ev','Eco','Erh'}
        Nf          % Number of frequencies
        Nx          % Number of unique x points
        Ny          % Number of unique y points
        Nang        % Number of point combinations (total number of directions)
        freqHz      % Frequency in Hz
        Directivity_dBi % Directivity in dBi
        Gain_dB         % Gain in dB 
        radEff_dB       % Radiation efficiency in dB    
        xRangeType      % Type of x-range: 'sym' or 'pos'
        yRangeType      % Type of y-range: 180 or 360
    end
    
    properties (Dependent = true, Hidden = true)
        symXZ % XZ plane symmetry
        symYZ % YZ plane symmetry
        symXY % XY plane symmetry
    end
    
    properties (SetAccess = private, Hidden = true)
        E3 = []         % Radial E-field component
        orientation     % Antenna orientation - altitude/azimuth in radians relative to zenith/North
        earthLocation   % Antenna location on the Earth longitude and latitude in radians, and height above sea level in meters - defaults to the roof of the Stellenbosch University E&E Engineering Dept. :)
        time            % Time as a datetime object (used, for instance, in astronomical observation)
        
        % Keep the input data here to not lose some info when going through
        % transformations and back...
        xBase
        yBase
        phBase
        thBase
        NxBase
        NyBase
        gridTypeBase
        E1Base
        E2Base
        E3Base
        coorTypeBase
        polTypeBase
    end
    
    properties (Constant = true, Hidden = true)
        c0 = physconst('Lightspeed');
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
        nSigDig = 8;
        projectionGrids = {'TrueView','Arcsin','Mollweide'};
        astroGrids = {'AzAlt','RAdec','GalLongLat'};
    end
    
    methods
        function obj = FarField(varargin)
            % FARFIELD class constructor method.
            % obj = FarField(x,y,E1,E2,freq,Prad,radEff,varargin) - can be
            % empty, which constructs a z-directed incremental dipole farfield pattern.
            % 
            % Inputs
            % - x: column vector of ph angles in rad, [Nang x 1]
            % - y: column vector of th angles in rad, [Nang x 1]
            % - E1: First E-field pattern component, [Nang x Nf]
            % - E2: Second E-field pattern component, [Nang x Nf]
            % - freq: Frequencies where the fields are defined in Hz, [1 x Nf]
            % - Prad: Radiated power at each frequency in W, [1 x Nf]
            % - radEff: Radiation efficiency at each frequency, [1 x Nf]
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   keywords and values are from the sets 
            %   -- E3:          Radial component of E-field ([]), [Nang x 1]
            %   -- coorType:    {('spherical')|'Ludwig1'|'Ludwig2AE'|'Ludwig2EA'|'Ludwig3'}
            %   -- polType:     {('linear')|'circular'|'slant'}
            %   -- gridType:    {('PhTh')|'DirCos'|'AzEl'|'ElAz'|'AzAlt'|'TrueView'|'ArcSin'|'Mollweide'|'RAdec'|'GalLongLat'}
            %   -- freqUnit:    {('Hz'),'kHz','MHz','GHz','THz'}
            %   -- symmetryXZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryYZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryXY:  {('none')|'electric'|'magnetic'}
            %   -- symBOR:      {('none')|'BOR0'|'BOR1'}
            %   -- r:           Radius where the E-Field is evaluated in m, (1)
            %   -- slant:       slant angle in rad for polType=slant, (pi/4)
            %   -- orientation: Antenna orientation, altitude/azimuth in radians relative to zenith/North ([0,pi/2]), [1 x 2]
            %   -- earthLocation: Antenna location, East-longitude/North-latitude in rad
            %                     and height above sea level in m, ([deg2rad(18.86) deg2rad(-33.93) 300])
            %   -- time:        Time, datetime object data type (datetime(2018,7,22,0,0,0))
            %
            % Outputs
            % - obj:    Farfield object
            %
            % Dependencies
            % -
            %
            % Created: 2019-05-09, Dirk de Villiers
            % Updated: 2019-05-09, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 0
            %   File : 
            %
            % Example
            %   F = FarField;
            %   F.plot
            
            % Set up defaults: z-directed incremental dipole
            [ph0,th0] = PhThGrid;
            f0 = 1e9;
            Eth0 = sqrt(3*obj.eta0).*sin(th0);
            Eph0 = zeros(size(Eth0));
            Prad0 = [];
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'FarField';
            
            % Optional inputs
            typeValidation_grid = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','ncols',1},'FarField');
            parseobj.addOptional('x',ph0,typeValidation_grid);
            parseobj.addOptional('y',th0,typeValidation_grid);
            
            typeValidation_fields = @(x) validateattributes(x,{'numeric'},{'finite','nonnan'},'FarField');
            parseobj.addOptional('E1',Eth0,typeValidation_fields);
            parseobj.addOptional('E2',Eph0,typeValidation_fields);
            
            typeValidation_freq = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','increasing','nrows',1},'FarField');
            parseobj.addOptional('freq',f0,typeValidation_freq);
            
            typeValidation_power = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nrows',1},'FarField');
            parseobj.addOptional('Prad',Prad0,typeValidation_power);
            parseobj.addOptional('radEff',1,typeValidation_power);
            
            % Name-value pairs
            parseobj.addParameter('E3',[],typeValidation_fields);
            
            expected_coorType = {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3'};
            parseobj.addParameter('coorType','spherical', @(x) any(validatestring(x,expected_coorType)));
            
            expected_polType = {'linear','circular','slant'};
            parseobj.addParameter('polType','linear', @(x) any(validatestring(x,expected_polType)));
            
            expected_gridType = {'PhTh','DirCos','AzEl','ElAz','AzAlt','TrueView','ArcSin','Mollweide','RAdec','GalLongLat'};
            parseobj.addParameter('gridType','PhTh', @(x) any(validatestring(x,expected_gridType)));
            
            expected_freqUnit = {'Hz','kHz','MHz','GHz','THz'};
            parseobj.addParameter('freqUnit','Hz', @(x) any(validatestring(x,expected_freqUnit)));
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR0','BOR1'};
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'FarField');
            parseobj.addParameter('r',1,typeValidation_scalar);
            parseobj.addParameter('slant',pi/4,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'FarField');
            parseobj.addParameter('orientation',[0,pi/2],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'FarField');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(varargin{:})
            
            obj.x = parseobj.Results.x;
            obj.y = parseobj.Results.y;
            obj.E1 = parseobj.Results.E1;
            obj.E2 = parseobj.Results.E2;
            obj.freq = parseobj.Results.freq;
            obj.Prad = parseobj.Results.Prad;
            obj.radEff = parseobj.Results.radEff;
            obj.E3 = parseobj.Results.E3;
            obj.coorType = parseobj.Results.coorType;
            obj.polType = parseobj.Results.polType;
            obj.gridType = parseobj.Results.gridType;
            obj.freqUnit = parseobj.Results.freqUnit;
            obj.symmetryXZ = parseobj.Results.symmetryXZ;
            obj.symmetryYZ = parseobj.Results.symmetryYZ;
            obj.symmetryXY = parseobj.Results.symmetryXY;
            obj.symmetryBOR = parseobj.Results.symmetryBOR;
            obj.r = parseobj.Results.r;
            obj.slant = parseobj.Results.slant;
            obj.orientation = parseobj.Results.orientation;
            obj.earthLocation = parseobj.Results.earthLocation;
            obj.time = parseobj.Results.time;
            
            % Set the base for power integration and symmetry checks
            obj = setBase(obj);
                
            % Check input sizes
            Nang = size(obj.x,1);
            Nf = numel(obj.freq);
            assert(size(obj.y,1)==Nang && size(obj.E1,1)==Nang && size(obj.E2,1)==Nang,'x,y,E1 and E2 must have the same number of rows')
            % Integrate power if none provided, but a full sphere of fields is 
            if isempty(obj.Prad)
                if obj.isGrid4pi
                    obj.Prad = obj.pradInt;
                else
                    obj.Prad = 4*pi;
                end
            end   
            if isscalar(obj.Prad) == 1, obj.Prad = repmat(obj.Prad,1,Nf); end
            if isscalar(obj.radEff) == 1, obj.radEff = repmat(obj.radEff,1,Nf); end
            assert(size(obj.E1,2)==Nf && size(obj.E2,2)==Nf && size(obj.Prad,2)==Nf && size(obj.radEff,2)==Nf,'E1, E2, freq, Prad and radEff must have the same number of columns')
            
            tol = 10^(-obj.nSigDig);
            % Check input symmetry validity
            if sum(abs([obj.symXZ,obj.symYZ,obj.symXY])) > 0 && ~strcmp(obj.symmetryBOR,'none')
                error('Cannot specify a plane- and BOR-symmetry for the same field')
            end
            if ~strcmp(obj.symmetryXZ,'none')
                obj1 = obj.grid2TrueView;
                assert(all(sign(obj1.y+tol) > 0) || all(sign(obj1.y-tol) < 0),'Invalid range for XZ symmetry')
            end
            if ~strcmp(obj.symmetryYZ,'none')
                obj1 = obj.grid2TrueView;
                assert(all(sign(obj1.x+tol) > 0) || all(sign(obj1.x-tol) < 0),'Invalid range for YZ symmetry')
            end
            if ~strcmp(obj.symmetryXY,'none')
                error('function: setSymmetryXY not implemented yet - please redefine with the full field and symmetryXY = none');
            end
            if strcmp(obj.symmetryBOR,'BOR0')
                assert(obj.Nx == 1,'Invalid range for BOR0 symmetry (one x-cut maximum)')
            elseif strcmp(obj.symmetryBOR,'BOR1')
                assert(all(abs(unique(obj.ph) - [0;pi/2]) < tol),'Invalid range for BOR1 symmetry (E-plane and H-plane required)')
            end
            
            obj = setBase(obj);
        end
        
        %% Dependency-based Setters
        function ph = get.ph(obj)
            [ph,~] = getPhThCurrent(obj);
        end
        
        function th = get.th(obj)
            [~,th] = getPhThCurrent(obj);
        end
        
        function Nf = get.Nf(obj)
            Nf = numel(obj.freq);
        end
        
        function Nx = get.Nx(obj)
            Nx = length(unique(obj.x));
        end
        
        function Ny = get.Ny(obj)
            Ny = length(unique(obj.y));
        end
        
        function Nang = get.Nang(obj)
            Nang = size(obj.x,1);
        end
        
        function freqHz = get.freqHz(obj)
            switch obj.freqUnit
                case 'Hz'
                    freqMult = 1;
                case 'kHz'
                    freqMult = 1e3;
                case 'MHz'
                    freqMult = 1e6;
                case 'GHz'
                    freqMult = 1e9;
                case 'THz'
                    freqMult = 1e12;
            end
            freqHz = obj.freq*freqMult;
        end
        
        function Directivity_dBi = get.Directivity_dBi(obj)
            Directivity_dBi = dB10(max(obj.getDirectivity()));
        end
        
        function Gain_dB = get.Gain_dB(obj)
            Gain_dB = dB10(max(obj.getGain()));
        end
        
        function radEff_dB = get.radEff_dB(obj)
            radEff_dB = dB10(obj.radEff);
        end
        
        function xname = get.xname(obj)
            [xname,~] = setXYnames(obj);
        end
        
        function yname = get.yname(obj)
            [~,yname] = setXYnames(obj);
        end
        
        function E1name = get.E1name(obj)
            [E1name,~] = setEnames(obj);
        end
        
        function E2name = get.E2name(obj)
            [~,E2name] = setEnames(obj);
        end
        
        function xRangeType = get.xRangeType(obj)
            [xRangeType,~] = setRangeTypes(obj);
        end
        
        function yRangeType = get.yRangeType(obj)
            [~,yRangeType] = setRangeTypes(obj);
        end
        
        function symXZ = get.symXZ(obj)
            switch obj.symmetryXZ
                case 'none'
                    symXZ = 0;
                case 'electric'
                    symXZ = -1;
                case 'magnetic'
                    symXZ = 1;
            end
        end
        
        function symYZ = get.symYZ(obj)
            switch obj.symmetryYZ
                case 'none'
                    symYZ = 0;
                case 'electric'
                    symYZ = -1;
                case 'magnetic'
                    symYZ = 1;
            end
        end
        
        function symXY = get.symXY(obj)
            switch obj.symmetryXY
                case 'none'
                    symXY = 0;
                case 'electric'
                    symXY = -1;
                case 'magnetic'
                    symXY = 1;
            end
        end
        
        function obj = setFreq(obj,freq,freqUnit)
            if nargin > 1
                assert(numel(freq) == size(obj.E1,2),'Error, freq must be the same length as the number of columns in E1')
                obj.freq = freq;
            end
            if nargin > 2
                obj.freqUnit = freqUnit;
            end
        end
        
        %% Pattern getters
        function FFpattern = getFarFieldStruct(obj)
            % GETFARFIELDSTRUCT Returns the legacy FarField struct data format.
            % FFpattern = getFarFieldStruct(obj) Converts a FarField
            % object into a struct data type. The struct format is compatible
            % with old script-based code.
            % 
            % Inputs
            % - obj:    FarField object
            % 
            % Outputs
            % - FFpattern: FarField struct
            %
            % Dependencies
            % -
            %
            % Created: 2019-05-09, Dirk de Villiers
            % Updated: 2019-05-09, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Fahmi Mokhupuki
            %  Level : 1
            %   File : 
            %
            % Example
            %   F = FarField;
            %   FFpattern = getFarFieldStruct(F);
            
            % This returns the legacy structure format for testing with all
            % the tons of old code
            
            obj = obj.coor2spherical(true);
            FFpattern.th = repmat(obj.y,1,obj.Nf);
            FFpattern.ph = repmat(obj.x,1,obj.Nf);
            FFpattern.Eth = obj.E1;
            FFpattern.Eph = obj.E2;
            FFpattern.freq = obj.freqHz;
            FFpattern.Nth = obj.Ny;
            FFpattern.Nph = obj.Nx;
            FFpattern.Nf = obj.Nf;
            FFpattern.Prad = obj.Prad;
        end
        
        function [E1field, E2field, E3field] = getEfield(obj)
            % GETEFIELD Returns E-field components from the radiated pattern.
            % [E1field, E2field, E3field] = getEfield(obj) Returns E-field
            % components from a FarField object in matrices of size [Nang x Nf]
            % E-field = E*exp(-jkr)/r
            %
            % Inputs
            % - obj:    FarField object
            %
            % Outputs
            % - E1field:    coorType dependent E-field component 
            % - E2field:    coorType dependent E-field component
            % - E3field:    coorType dependent E-field component ([])
            %
            % Dependencies
            % -
            %
            % Created: 2019-05-07, Dirk de Villiers
            % Updated: 2019-05-07, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Fahmi Mokhupuki
            %  Level : 1
            %   File : 
            %
            % Example
            %   F = FarField;
            %   [E1,E2,E3] = getEfield(F)
            
            % function [E1field, E2field, E3field] = getEfield(obj)
            % Returns the Efield matrices of size [Nang x Nf]
            % Efield = E*exp(-jkr)/r
            k = 2.*pi.*obj.freqHz./obj.c0;
            FFfact = exp(-1i.*k.*obj.r)./obj.r;
            E1field = bsxfun(@times,obj.E1,FFfact);
            E2field = bsxfun(@times,obj.E2,FFfact);
            if ~isempty(obj.E3)
                E3field = bsxfun(@times,obj.E3,FFfact);
            else
                E3field = [];
            end
        end
        
        function [W] = getW(obj)
            % GETW Returns the radiation density in W/m2.
            
            % function [W] = getW(obj)
            % returns the radiation density in W/m2 [Nang x Nf]
            [E1f, E2f] = getEfield(obj);    % Can use any orthogonal pair
            W = 1./(2.*obj.eta0).*(abs(E1f).^2 + abs(E2f).^2);
        end
        
        function [U] = getU(obj)
            % GETU Returns the radiation intensity in W/unit solid angle.
            
            % function [U] = getU(obj)
            % returns the radiation intensity in W/unit solid angle [Nang x Nf]
            U = obj.r^2.*getW(obj);
        end
        
        function [D] = getDirectivity(obj)
            % GETDIRECTIVITY Returns the directivity (linear).
            
            % function [D] = getDirectivity(obj)
            % returns the directivity (linear) in D [Nang x Nf]
            D = 4.*pi.*bsxfun(@times,getU(obj),1./obj.Prad);
        end
        
        function [G] = getGain(obj)
            % GETGAIN Returns the gain in a (linear).
            
            % function [G] = getGain(obj)
            % returns the gain (linear) in G [Nang x Nf]
            G = bsxfun(@times,getDirectivity(obj),obj.radEff);
        end
        
        function [AR] = getAxialRatio(obj)
            % GETAXIALRATIO Returns the Axial Ratio (linear).
            
            % function [AR] = getAxialRatio(obj)
            % returns the Axial Ratio (linear) [Nang x Nf]
            AR = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        function [ARinv] = getAxialRatioInv(obj)
            % GETAXIALRATIOINV Returns the inverted Axial Ration (linear).
            
            % function [ARinv] = getAxialRatioInv(obj)
            % returns the inverted Axial Ratio in ARinv [Nang x Nf]
            ARinv = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        function [Xpol] = getCO_XP(obj)
            % GETCO_XP Returns the CO/XP ratio (linear).
            
            % function [Xpol] = getCO_XP(obj)
            % returns the CO/XP ratio (linear) [Nang x Nf]
            Xpol = (abs(obj.E2)./abs(obj.E1)).^2;
        end
        
        function [Xpol] = getXP_CO(obj)
            % GETXP_CO Returns the XP/CO ratio (linear)
            
            % function [Xpol] = getXP_CO(obj)
            % returns the XP/CO ratio (linear) [Nang x Nf]
            Xpol = (abs(obj.E1)./abs(obj.E2)).^2;
        end
        
        %% Grid transformation setters
        function obj = changeGrid(obj,gridTypeString)
            % CHANGEGRID Change the current FarField object grid.  
            
            mustBeMember(gridTypeString, {'PhTh','DirCos','AzEl','ElAz','TrueView','ArcSin','AzAlt','RAdec','GalLongLat'});
            handleGridType = str2func(['grid2',gridTypeString]);
            obj = handleGridType(obj);  
        end
        
        function obj = grid2PhTh(obj)
            % GRID2PHTH Change the current FarField object grid to a PhTh grid.
            
            formerGridType = obj.gridType;
            formerNx = obj.Nx;
            formerNy = obj.Ny;
            if any(strcmp(obj.gridType,obj.astroGrids))
                obj = obj.grid2AzAlt;
            end
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'PhTh')
                [obj.x,obj.y] = getPhTh(obj);
                obj.gridType = 'PhTh';
                if strcmp(formerGridType,'AzAlt')
                    % Get the grid step sizes from the original
                    obj = obj.rotate(@rotEulersph,[0,0,obj.orientation(1)]);
                    obj = obj.rotate(@rotEulersph,[0,-pi/2+obj.orientation(2),0]);
                    xmin = min(obj.x);
                    xmax = max(obj.x);
                    ymin = min(obj.y);
                    ymax = max(obj.y);
                    stepx = (max(obj.x) - min(obj.x))./(formerNx-1);
                    stepy = (max(obj.y) - min(obj.y))./(formerNy-1);
                    stepDeg = rad2deg([stepx,stepy]);
                    % Set the baseGrid of the rotated object.  This is required
                    % since all transformations operate from the base grid
                    obj = obj.sortGrid;
                    obj = obj.setBase;
                    obj = obj.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
                end
            end
        end
        
        function obj = grid2DirCos(obj)
            % GRID2DIRCOS Change the current FarField object grid to a DirCos grid.
            
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'DirCos')
                [obj.x,obj.y] = getDirCos(obj);
                obj.gridType = 'DirCos';
            end
        end
        
        function obj = grid2AzEl(obj)
            % GRID2AZEL Change the current FarField object grid to a AzEl grid.
            
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'AzEl')
                [obj.x,obj.y] = getAzEl(obj);
                obj.gridType = 'AzEl';
            end
        end
        
        function obj = grid2ElAz(obj)
            % GRID2ELAZ Change the current FarField object grid to a ElAz grid.
            
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'ElAz')
                [obj.x,obj.y] = getElAz(obj);
                obj.gridType = 'ElAz';
            end
        end
        
        function obj = grid2TrueView(obj)
            % GRID2TRUEVIEW Change the current FarField object grid to a TrueView grid.
            
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'TrueView')
                [obj.x,obj.y] = getTrueView(obj);
                obj.gridType = 'TrueView';
            end
        end
        
        function obj = grid2ArcSin(obj)
            % GRID2ARCSIN Change the current FarField object grid to an ArcSin grid.
            
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'ArcSin')
                [obj.x,obj.y] = getArcSin(obj);
                obj.gridType = 'ArcSin';
            end
        end
        
        function obj = grid2AzAlt(obj)
            % GRID2AZALT Change the current FarField object grid to an AzAlt grid.
            
            formerNx = obj.Nx;
            formerNy = obj.Ny;
            if ~any(strcmp([obj.projectionGrids,obj.astroGrids],obj.gridType))
                currGridType = obj.gridType;
                if ~strcmp(currGridType,'PhTh')
                    obj = obj.grid2PhTh;
                end
                obj = obj.rotate(@rotGRASPsph,[wrap2pi(pi/2-obj.orientation(2)),-obj.orientation(1),0]);
                eval(['obj = grid2',currGridType,'(obj);']);
            end
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'AzAlt')
                [obj.x,obj.y] = getAzAlt(obj);
                obj.gridType = 'AzAlt';
            end
            % Get the grid step sizes from the original
            xmin = min(obj.x);
            xmax = max(obj.x);
            ymin = min(obj.y);
            ymax = max(obj.y);
            stepx = (max(obj.x) - min(obj.x))./(formerNx-1);
            stepy = (max(obj.y) - min(obj.y))./(formerNy-1);
            stepDeg = rad2deg([stepx,stepy]);
            % Set the baseGrid of the rotated object.  This is required
            % since all transformations operate from the base grid
            obj = obj.sortGrid;
            obj = obj.setBase;
            obj = obj.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
        end
        
        function obj = grid2RAdec(obj)
            % GRID2RAdec Change the current FarField object grid to a RAdec grid.
            
            formerNx = obj.Nx;
            formerNy = obj.Ny;
            assert(any(strcmp(obj.gridTypeBase,[obj.projectionGrids,obj.astroGrids])),'Current grid must be a projection or be in an astronomical reference frame')
            %must also check RA/dec lengths here...
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'RAdec')
                [obj.x,obj.y] = getRAdec(obj);
                obj.gridType = 'RAdec';
            end
            % Get the grid step sizes from the original
            xmin = min(obj.x);
            xmax = max(obj.x);
            ymin = min(obj.y);
            ymax = max(obj.y);
            stepx = (max(obj.x) - min(obj.x))./(formerNx-1);
            stepy = (max(obj.y) - min(obj.y))./(formerNy-1);
            stepDeg = rad2deg([stepx,stepy]);
            % Set the baseGrid of the rotated object.  This is required
            % since all transformations operate from the base grid
            obj = obj.sortGrid;
            obj = obj.setBase;
            obj = obj.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
        end
        
        function obj = grid2GalLongLat(obj)
            % GRID2GALLONGLAT Change the current FarField object grid to a GalLongLat grid.
            
            formerNx = obj.Nx;
            formerNy = obj.Ny;
            assert(any(strcmp(obj.gridTypeBase,[obj.projectionGrids,obj.astroGrids])),'Current grid must be a projection or be in an astronomical reference frame');
            %must also check RA/dec lengths here...
            obj = obj.grid2Base;
            if ~strcmp(obj.gridType,'GalLongLat')
                [obj.x,obj.y] = getGalLongLat(obj);
                obj.gridType = 'GalLongLat';
            end
            % Get the grid step sizes from the original
            xmin = min(obj.x);
            xmax = max(obj.x);
            ymin = min(obj.y);
            ymax = max(obj.y);
            stepx = (max(obj.x) - min(obj.x))./(formerNx-1);
            stepy = (max(obj.y) - min(obj.y))./(formerNy-1);
            stepDeg = rad2deg([stepx,stepy]);
            % Set the baseGrid of the rotated object.  This is required
            % since all transformations operate from the base grid
            obj = obj.sortGrid;
            obj = obj.setBase;
            obj = obj.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
        end
        
        %% Grid range shifters
        function obj = sortGrid(obj)
            % SORTGRID Sort grid with corresponding E-field values 
            % in ascending order, according to x then y.
            
            obj = roundGrid(obj);
            [~,iSort] = sortrows([obj.x,obj.y],[1 2]);
            obj.x = obj.x(iSort);
            obj.y = obj.y(iSort);
            obj.E1 = obj.E1(iSort,:);
            obj.E2 = obj.E2(iSort,:);
            if ~isempty(obj.E3), obj.E3 = obj.E3(iSort,:); end
        end
        
        function obj = roundGrid(obj,nSigDig)
            % ROUNDGRID Round grid entries to some significant number of digits.
            
            % Round to some significant digits for sorting (some issues can
            % arise in deg2rad and rad2deg
            if nargin < 2
                nSigDig = obj.nSigDig;
            end
            xRound = round(obj.x*10^nSigDig)/10^nSigDig;
            yRound = round(obj.y*10^nSigDig)/10^nSigDig;
            obj.x = xRound;
            obj.y = yRound;
        end
        
        function obj = copyAndInsertXcut(obj1,xvalCopy,xvalPaste,tol)
            % COPYANDINSERTXCUT Copy a FarField cut into another position.
            
            % Use this to copy an X cut into another position.  Typically
            % handy when some transformation does not include the closing
            % cut - that is the 0 and 360 or -180 and 180 cuts.  Can in
            % principle be used to do random stuff - so careful.
            
            if nargin < 4
                tol = mean(diff(unique(obj1.x)));
            end
            % Make a whole new object to initialise the base
            % correctly
            inInd = find(abs(obj1.x - xvalCopy) < tol);
            xNew = [obj1.x;xvalPaste.*ones(size(inInd))];
            yNew = [obj1.y;obj1.y(inInd)];
            E1New = [obj1.E1;obj1.E1(inInd,:)];
            E2New = [obj1.E2;obj1.E2(inInd,:)];
            if ~isempty(obj1.E3)
                E3New = [obj1.E3;obj1.E3(inInd,:)];
            else
                E3New = [];
            end
            obj = FarField(xNew,yNew,E1New,E2New,E3New,obj1.freq,...
                obj1.Prad.*2,obj1.radEff,obj1.coorType,obj1.polType,obj1.gridType,obj1.freqUnit,obj1.slant);
            obj = obj.sortGrid;
            obj = FarField(obj.x,obj.y,obj.E1,obj.E2,obj.E3,obj.freq,...
                obj.Prad.*2,obj.radEff,obj.coorType,obj.polType,obj.gridType,obj.freqUnit,obj.slant);
        end
        
        function obj = setXrange(obj,type)
            % SETXRANGE Set the x-range (ph, az or ep) for the angular gridTypes
            % in the FarField object.
            
            % Attempts to set the x-range (ph, az, or ep) for the angular
            % gridTypes in the FarField object
            % type = 'pos':
            % [-180,180] is transformed to [0,360]
            % with the redundant -180 cut replaced by a redundant 360 cut.
            % type = 'sym':
            % [0,360] is transformed to [-180,180]
            % with the redundant 360 cut replaced by a redundant -180 cut.
            % The resulting object has the same number
            % of field points as the input object.
            
            mustBeMember(type,{'pos','sym'})
            if strcmp(type,'pos')
                t = 'p';
                if strcmp(obj.xRangeType,'pos'), return; end
            elseif strcmp(type,'sym')
                t = 's';
                if strcmp(obj.xRangeType,'sym'), return; end
            end
            % tol = 1e-10;
            tol = 10^(-obj.nSigDig);
            if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                if t == 'p'
                    %         iout = find(obj.x == -pi);   % Redundant
                    %         iin = find(obj.x == 0);      % Will become redundant after inserting
                    iout = find(abs(obj.x + pi) < tol);   % Redundant
                    iin = find(abs(obj.x - 0) < tol);      % Will become redundant after inserting
                    if strcmp(obj.xRangeType,'pos'), return; end
                elseif t == 's'
                    %         iout = find(obj.x == 2*pi);   % Redundant
                    %         iin = find(obj.x == pi);      % Will become redundant after inserting
                    iout = find(abs(obj.x - 2*pi) < tol);   % Redundant
                    iin = find(abs(obj.x - pi) < tol);      % Will become redundant after inserting
                    if strcmp(obj.xRangeType,'sym'), return; end
                end
                Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add
                % Truncate the indexes to the shortest length
                iout = iout(1:Nredun);
                iin = iin(1:Nredun);
                yAdd = obj.y(iin);
                E1Add = obj.E1(iin,:);
                E2Add = obj.E2(iin,:);
                redunFound = Nredun > 0;
                % First remove the ph=-180 from the matrix...
                if redunFound
                    obj.x(iout) = [];
                    obj.y(iout) = [];
                    obj.E1(iout,:) = [];
                    obj.E2(iout,:) = [];
                end
                % Now shift the x values
                if t == 'p'
                    outOfRangeInd = find(obj.x < 0);
                    obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + 2*pi;
                elseif t == 's'
                    outOfRangeInd = find(obj.x > pi);
                    obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - 2*pi;
                end
                % if the -180 was removed, add a 360 cut
                if redunFound
                    if t == 'p'
                        xAdd = ones(numel(iin),1).*2.*pi;
                    elseif t == 's'
                        xAdd = ones(numel(iin),1).*-pi;
                    end
                    obj.x = [obj.x;xAdd(1:Nredun)];
                    obj.y = [obj.y;yAdd(1:Nredun)];
                    obj.E1 = [obj.E1;E1Add(1:Nredun,:)];
                    obj.E2 = [obj.E2;E2Add(1:Nredun,:)];
                end
                % Sort
                obj = obj.sortGrid;
            else
                warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
            end
            
        end
        
        function obj = setYrange(obj,type)
            % SETXRANGE Set the y-range (th, el, or al) for the angular gridTypes
            % in the FarField object.
            
            % Attempts to set the y-range (th, el, or al) for the angular
            % gridTypes in the FarField object
            % type = 180:
            % [0,180] gridType = 'PhTh'
            % [-90,90] gridType = 'AzEl' | 'ElAz'
            % type = 360:
            % [0,360] for xRangeType = 'pos' - xRange = [0,180]
            % [-180,180] for xRangeType = 'sym' - xRange = [-90,90]
            % Redundant fields are replaced by valid ones as far as possible
            % The resulting object has the same number
            % of field points as the input object.
            
            assert(type == 180 || type == 360,'Unknown type: Should be 180 or 360');
            
            % [iout,iin,xAdd,yAdd,E1Add,E2Add,E3Add] = deal([]);
            % Nredun = 0;
            
            % eps = 1e-10;
            eps = 1e-1*min([diff(unique(obj.x));diff(unique(obj.y))]);
            
            % For all cases some redundant points will be removed, and some new ones
            % inserted into the grid. It depends on how the space is cut and rotated
            % where to put in and take out points...  Done now on a case-by-case basis
            % because I can't think of a general algorithm yet
            if strcmp(obj.gridType,'PhTh')
                if type == 360
                    if strcmp(obj.xRangeType,'pos')
                        iout = find(obj.x  > (pi+eps) & abs(obj.y - pi) < eps);   % Redundant
                        iin = find(abs(obj.x - 0) < eps);                   % Will become redundant after inserting
                        xAdd = obj.x(iin);
                        yAdd = 2*pi - obj.y(iin);
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd = find(obj.x > (pi+eps));
                        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - pi;
                        obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
                    elseif strcmp(obj.xRangeType,'sym')
                        iout = [find((obj.x >= (-pi/2-eps) & obj.x <= (pi/2+eps)) & abs(obj.y - 0) < eps); find(abs(obj.x + pi) < eps)];   % Redundant
                        iin = find(abs(obj.x + pi/2) < eps | abs(obj.x - pi/2) < eps);
                        xAdd = obj.x(iin) - sign(obj.x(iin)).*pi;
                        yAdd = -obj.y(iin);
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd = find(obj.x > (pi/2+eps) | obj.x < (-pi/2-eps));
                        signPos = sign(obj.x(outOfRangeInd));
                        signPos(signPos == 0) = 1;
                        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
                        obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
                    end
                elseif type == 180
                    if strcmp(obj.xRangeType,'pos')
                        iout = find(obj.y <= (pi+eps) & abs(obj.x - pi) < eps);   % Redundant
                        iin = find(obj.x <= (pi+eps) & abs(obj.y - pi) < eps);    % Will become redundant after inserting
                        xAdd = obj.x(iin) + pi;
                        yAdd = 2*pi - obj.y(iin);
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd = find(obj.y > (pi+eps));
                        obj.y(outOfRangeInd) = 2*pi - obj.y(outOfRangeInd);
                        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) + pi;
                    elseif strcmp(obj.xRangeType,'sym')
                        %             iout = unique(find((abs(obj.x - pi/2) < eps | abs(obj.x + pi/2) < eps) & (obj.y <= eps)));   % Redundant
                        iout = find((abs(obj.x - pi/2) < eps | abs(obj.x + pi/2) < eps) & (obj.y <= eps));   % Redundant
                        iin1 = find(abs(obj.y - 0) < eps);
                        iin2 = find(abs(obj.x - 0) < eps & obj.y <= eps);
                        iin = [iin1;iin2];
                        signPosIn = sign(obj.x(iin1));
                        signPosIn(signPosIn == 0) = 1;
                        xAdd1 = obj.x(iin1) - signPosIn.*pi;
                        xAdd2 = obj.x(iin2) + pi;
                        xAdd = [xAdd1;xAdd2];
                        yAdd = -obj.y(iin);
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd = find(obj.y < -eps);
                        obj.y(outOfRangeInd) = -obj.y(outOfRangeInd);
                        signPos = sign(obj.x(outOfRangeInd));
                        signPos(signPos == 0) = 1;
                        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
                    end
                end
            elseif strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                if type == 360
                    if strcmp(obj.xRangeType,'pos')
                        iout = find(obj.x <= pi+eps & (abs(obj.y - pi/2) < eps | abs(obj.y + pi/2) < eps));   % Redundant
                        iin1 = find(abs(obj.x - pi) < eps);                   % Will become redundant after inserting
                        iin2 = find(abs(obj.y - 0) < eps);
                        iin = [iin1;iin2];
                        xAdd = [pi - obj.x(iin1);obj.x(iin2)];
                        yAdd = [pi + obj.y(iin1);2*pi + obj.y(iin2)];
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift - first x > 180
                        outOfRangeInd1 = find(obj.x > (pi+eps));
                        obj.x(outOfRangeInd1) = obj.x(outOfRangeInd1) - pi;
                        obj.y(outOfRangeInd1) = pi - obj.y(outOfRangeInd1);
                        % Now for negative y
                        outOfRangeInd2 = find(obj.y < -eps);
                        obj.x(outOfRangeInd2) = obj.x(outOfRangeInd2);
                        obj.y(outOfRangeInd2) = 2*pi + obj.y(outOfRangeInd2);
                    elseif strcmp(obj.xRangeType,'sym')
                        iout = find(abs(obj.x - pi) < eps | (obj.x <= (pi/2+eps) & obj.x >= (-pi/2-eps) & abs(obj.y - pi/2) < eps) | (obj.x <= (pi/2+eps) & obj.x >= (-pi/2-eps) & abs(obj.y + pi/2) < eps));   % Redundant
                        iin1 = find(abs(obj.x + pi/2) < eps);
                        iin2 = find(abs(obj.x - pi/2) < eps);
                        iin3 = find(abs(obj.y - 0) < eps & (obj.x >= (pi/2-eps) | obj.x <= (-pi/2+eps)));
                        xAdd1 = obj.x(iin1);
                        xAdd2 = obj.x(iin2);
                        xAdd3 = obj.x(iin3) - sign(obj.x(iin3)).*pi;    % Never zero so don't worry about making positive
                        xAdd = [xAdd1;xAdd2;xAdd3];
                        iin12 = [iin1;iin2];
                        signPos = sign(obj.y(iin12));
                        signPos(signPos == 0) = 1;
                        yAdd12 = [signPos.*pi - obj.y(iin12)];  % Include the -180 y-cut
                        yAdd3 = obj.y(iin3).*0 - pi;
                        yAdd = [yAdd12;yAdd3];
                        iin = [iin12;iin3];
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd = find(obj.x > (pi/2+eps) | obj.x < (-pi/2-eps));
                        signPos = sign(obj.x(outOfRangeInd));
                        signPos(signPos == 0) = 1;
                        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
                        signPos = sign(obj.y(outOfRangeInd));
                        signPos(signPos == 0) = 1;
                        obj.y(outOfRangeInd) = signPos.*pi - obj.y(outOfRangeInd);
                    end
                elseif type == 180
                    if strcmp(obj.xRangeType,'pos')
                        iout = find(abs(obj.y - 2*pi) < eps | (abs(obj.x - 0) < eps & obj.y > (pi/2+eps) & obj.y <= (3*pi/2+eps)));   % Redundant
                        iin1 = find(abs(obj.y - 3*pi/2) < eps);
                        iin2 = find(abs(obj.y - pi/2) < eps & obj.x > eps);
                        iin = [iin1;iin2];
                        xAdd = [obj.x(iin1);obj.x(iin2) + pi];
                        yAdd = [obj.y(iin1) - 2*pi; pi - obj.y(iin2)];
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd1 = find(obj.y > (3*pi/2 + eps));
                        obj.x(outOfRangeInd1) = obj.x(outOfRangeInd1);
                        obj.y(outOfRangeInd1) = obj.y(outOfRangeInd1) - 2*pi;
                        outOfRangeInd2 = find(obj.y > (pi/2 + eps));
                        obj.x(outOfRangeInd2) = obj.x(outOfRangeInd2) + pi;
                        obj.y(outOfRangeInd2) = pi - obj.y(outOfRangeInd2);
                    elseif strcmp(obj.xRangeType,'sym')
                        %             iout = unique(find(abs(obj.y - pi) < eps | (abs(abs(obj.x) - pi/2) < eps & abs(obj.y) >= (pi/2-eps))));   % Redundant
                        iout = find(abs(obj.y - pi) < eps | (abs(abs(obj.x) - pi/2) < eps & abs(obj.y) >= (pi/2-eps)));   % Redundant
                        iin1 = find(abs(obj.x - 0) < eps & (obj.y >= (pi/2 - eps) | (obj.y <= -(pi/2 - eps) & obj.y > -(pi-eps))));
                        iin2 = find(abs(abs(obj.y) - pi/2) < eps);
                        iin = [iin1;iin2];
                        signPos1 = sign(obj.x(iin1));
                        signPos1(signPos1 == 0) = -1; % Add these to the x = 180 cut...
                        xAdd1 = obj.x(iin1) - signPos1.*pi;
                        signPos2 = sign(obj.x(iin2));
                        signPos2(signPos2 == 0) = 1;
                        xAdd2 = obj.x(iin2) - signPos2.*pi;
                        xAdd = [xAdd1;xAdd2];
                        yAdd = sign(obj.y(iin)).*pi - obj.y(iin);
                        obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        % Apply shift
                        outOfRangeInd = find(obj.y > (pi/2 + eps) | obj.y < (-pi/2 - eps));
                        signPos = sign(obj.x(outOfRangeInd));
                        signPos(signPos == 0) = 1;
                        obj.x(outOfRangeInd) = obj.x(outOfRangeInd) - signPos.*pi;
                        obj.y(outOfRangeInd) = sign(obj.y(outOfRangeInd)).*pi - obj.y(outOfRangeInd); % Never 0 so no need to worry about fixing that
                    end
                end
            else
                warning(['Cant shift a polar grid like ', obj.gridType, ' on a cartesian grid']);
            end
            % Sort
            obj = roundGrid(obj,6);
            obj = obj.sortGrid;
        end
        
        
        %% Coordinate system transformation methods
        function obj = changeCoor(obj,coorTypeString,setStdGrid)
            % CHANGECOOR Change the FarField object coordinate type.
            
            if nargin < 3, setStdGrid = true; end
            mustBeMember(coorTypeString, {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3'});
            handleCoorType = str2func(['coor2',coorTypeString]);
            obj = handleCoorType(obj,setStdGrid); 
        end
        
        function obj = coor2spherical(obj,setStdGrid)
            % COOR2SPHERICAL Change the current FarField object coordinate
            % type to spherical coordinates.
            
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'spherical')
                [obj.E1,obj.E2] = getEspherical(obj);
                obj.coorType = 'spherical';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig1(obj,setStdGrid)
            % COOR2LUDWIG1 Change the current FarField object coordinate
            % type to Ludwig1 coordinates.
            
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig1')
                [obj.E1,obj.E2] = getELudwig1(obj);
                obj.coorType = 'Ludwig1';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig2AE(obj,setStdGrid)
            % COOR2LUDWIG2AE Change the current FarField object coordinate
            % type to Ludwig2AE coordinates.
            
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig2AE')
                [obj.E1,obj.E2] = getELudwig2AE(obj);
                obj.coorType = 'Ludwig2AE';
            end
            if setStdGrid
                obj = obj.grid2AzEl;
            end
        end
        
        function obj = coor2Ludwig2EA(obj,setStdGrid)
            % COOR2LUDWIG2EA Change the current FarField object coordinate
            % type to Ludwig2EA coordinates.
            
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig2EA')
                [obj.E1,obj.E2] = getELudwig2EA(obj);
                obj.coorType = 'Ludwig2EA';
            end
            if setStdGrid
                obj = obj.grid2ElAz;
            end
        end
        
        function obj = coor2Ludwig3(obj,setStdGrid)
            % COOR2LUDWIG2EA Change the current FarField object coordinate
            % type to Ludwig3 coordinates.
            
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig3')
                [obj.E1,obj.E2] = getELudwig3(obj);
                obj.coorType = 'Ludwig3';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        %% polarization type transformation methods
        function obj = changePol(obj,polTypeString)
            % CHANGEPOL Change the FarField object polarization type.
            
            mustBeMember(polTypeString, {'linear','circular','slant'});
            handlePolType = str2func(['pol2',polTypeString]);
            obj = handlePolType(obj);   
        end
        
        function obj = pol2linear(obj)
            % POL2LINEAR Change the FarField object polarization to linear
            % polarization
            
            if ~strcmp(obj.polType,'linear')
                [obj.E1, obj.E2] = getElin(obj);
                obj.polType = 'linear';
            end
        end
        
        function obj = pol2circular(obj)
            % POL2CIRCULAR Change the FarField object polarization to circular
            % polarization
            
            if ~strcmp(obj.polType,'circular')
                [obj.E1,obj.E2] = getEcircular(obj);
                obj.polType = 'circular';
            end
        end
        
        function obj = pol2slant(obj)
            % POL2SLANT Change the FarField object polarization to slant
            % polarization
            
            if ~strcmp(obj.polType,'slant')
                [obj.E1,obj.E2] = getEslant(obj);
                obj.polType = 'slant';
            end
        end
        
        %% Format transformation
        function obj = transformTypes(obj, obj1)
            % TRANSFORMTYPES Transform the properties of obj to that of obj1.
            
            % Function to transform the format of obj to that of obj1 -
            % that is the grid,coor, and pol Types of obj goes to those of
            % obj1.
            objGridType = obj1.gridType;
            objCoorType = obj1.coorType;
            objPolType = obj1.polType;
            handleGridType = str2func(['grid2',objGridType]);
            handleCoorType = str2func(['coor2',objCoorType]);
            handlePolType = str2func(['pol2',objPolType]);
            obj = handleGridType(obj);
            obj = handleCoorType(obj,0);
            obj = handlePolType(obj);
        end
        
        %% Base grid functions
        function obj = reset2Base(obj)
            % RESET2BASE Hard reset the FarField object to its initial format
            % when it was first initiated.
            
            % Hard reset to the base format
            obj.x = obj.xBase;
            obj.y = obj.yBase;
            obj.gridType = obj.gridTypeBase;
            obj.E1 = obj.E1Base;
            obj.E2 = obj.E2Base;
            obj.E3 = obj.E3Base;
            obj.coorType = obj.coorTypeBase;
            obj.polType = obj.polTypeBase;
        end
        
        function obj = grid2Base(obj)
            % GRID2BASE Evaluate the current object (pol and coor) on the
            % base grid.
            
            % Evaluate the current object (pol and coor) on the base grid
            coorTypeIn = obj.coorType;
            polTypeIn = obj.polType;
            coorTypeH = str2func(['coor2',coorTypeIn]);
            polTypeH = str2func(['pol2',polTypeIn]);
            obj = obj.reset2Base;
            % Keep the current coorType and polType
            obj = coorTypeH(obj,false);
            obj = polTypeH(obj);
        end
        
        function obj = currentForm2Base(obj1,stepDeg,xylimsDeg,hemisphere)
            % CURRENTFORM2BASE Sets the FarField object base to the current format.
            
            % Sets the base to the current format. Resamples the field on a
            % regular plaid grid, and makes this the new base grid. This is
            % typically then not where actual samples where, but instead
            % interpolated values.
            % Format of xylimsDeg is [xmin xmax; ymin ymax]
            
            % Due to the interpolation only operating on the y in 180 range
            % from the Direction Cosine transforms...
            assert(~isequal(obj1.yRangeType,'360'),'yRangeType cannot be 360 for the base representation')
            
            % Set defaults
            if strcmp(obj1.gridTypeBase,'DirCos') || strcmp(obj1.gridType,'ArcSin')
                stepX = asin(min(abs(diff(unique(obj1.xBase)))));
                stepY = asin(min(abs(diff(unique(obj1.yBase)))));
                xmin = min(asin(obj1.x));
                xmax = max(asin(obj1.x));
                ymin = min(asin(obj1.y));
                ymax = max(asin(obj1.y));
            else
                % Sort out rounding errors for degrees
                stepX = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.xBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                stepY = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.yBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                xmin = deg2rad(round(rad2deg(min(obj1.x))*10^obj1.nSigDig)./10^obj1.nSigDig);
                xmax = deg2rad(round(rad2deg(max(obj1.x))*10^obj1.nSigDig)./10^obj1.nSigDig);
                ymin = deg2rad(round(rad2deg(min(obj1.y))*10^obj1.nSigDig)./10^obj1.nSigDig);
                ymax = deg2rad(round(rad2deg(max(obj1.y))*10^obj1.nSigDig)./10^obj1.nSigDig);
            end
            
            hem = 'top';
            
            % Overwrite defaults
            if nargin >= 2
                if numel(stepDeg) == 1
                    [stepX,stepY] = deal(deg2rad(stepDeg));
                elseif numel(stepDeg) == 2
                    stepX = deg2rad(stepDeg(1));
                    stepY = deg2rad(stepDeg(2));
                end
                if nargin >= 3
                    xylim = deg2rad(xylimsDeg);
                    xmin = xylim(1,1);
                    xmax = xylim(1,2);
                    ymin = xylim(2,1);
                    ymax = xylim(2,2);
                    if nargin == 4
                        hem = hemisphere;
                    end
                end
            end
            
            if strcmp(obj1.gridType,'DirCos') || strcmp(obj1.gridType,'ArcSin')
                stepX = sin(stepX);
                stepY = sin(stepY);
                xmin = sin(xmin);
                xmax = sin(xmax);
                ymin = sin(ymin);
                ymax = sin(ymax);
            end
            
            % Build the new grid
            Nxi = round((xmax - xmin)/stepX) + 1;
            Nyi = round((ymax - ymin)/stepY) + 1;
            xivect = linspace(xmin,xmax,Nxi);
            yivect = linspace(ymin,ymax,Nyi);
            [Xi,Yi] = meshgrid(xivect,yivect);
            xi = Xi(:);
            yi = Yi(:);
            
            % Interpolate the fields
            [E1grid,E2grid] = deal(zeros(Nxi*Nyi,obj1.Nf));
            for ff = 1:obj1.Nf
                E1grid(:,ff) = interpolateGrid(obj1,'E1',xi,yi,ff,hem);
                E2grid(:,ff) = interpolateGrid(obj1,'E2',xi,yi,ff,hem);
            end
            % Remove the extra phase introduced by the interpolateGrid
            % function - this just keeps the real/imag and phase field
            % consistant with the plotting
            k = 2.*pi.*obj1.freqHz./obj1.c0;
            FFfact = exp(1i.*k.*obj1.r)./obj1.r;
            E1grid = bsxfun(@times,E1grid,FFfact);
            E2grid = bsxfun(@times,E2grid,FFfact);
            
            % Populate the new farField object
            obj = obj1;
            obj.x = xi(:);
            obj.y = yi(:);
            obj.E1 = E1grid;
            obj.E2 = E2grid;
            obj = setBase(obj);
            % Update the current form to the base form
            obj = reset2Base(obj);
        end
        
        
        %% Plotting methods
        function plot(obj,varargin)
            
            %PLOT   Plots a FarField object.
            % plot(obj,varargin) plots a 1-D, 2-D, or 3-D representation
            % of a FarField object.
            %
            % Inputs
            % - obj: FarField object
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   keywords and values are from the sets
            %   -- freqIndex:   The index of the frequency to be plotted (1)
            %   -- plotType:    {('3D') | '2D' | 'polar' | 'cartesian'}
            %   -- output:      {('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio'
            %                   | 'AxialRatioInv', 'CO_XP' | 'XP_CO' | 'W' | 'U'}
            %   -- outputType:  {('mag') | 'phase' | 'real' | 'imag'}
            %
            % Outputs
            % - []
            %
            % Created: 2019-05-06, Dirk de Villiers
            % Updated: 2019-05-06, Dirk de Villiers
            %
            % Tested : Matlab R2018b, Dirk de Villiers
            %  Level : 0
            %   File : \testScripts\testScript_FarField.m
            %
            % Example
            %   FF = FarField;
            %   FF.plot('plotType','2D','step',1)
            
            
            
            % function plot(obj,name,value)
            % Plots a representation of the farfield object in obj
            % name, value are name value pairs, and can be the following:
            %
            % freqIndex is the index of the frequency to be plotted (default 1)
            %
            % plotType can be:
            %   ('3D') | '2D' | 'polar' | 'cartesian'
            %
            % output can be:
            %   ('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio' | 'AxialRatioInv'
            %   'CO_XP' | 'XP_CO' | 'W' | 'U'
            %
            % outputType can be:
            %   ('mag') | 'phase' | 'real' | 'imag' last 3 only used for on E-field plots
            %
            % norm is a boolean (false) to normalize to maximum magnitude
            %
            % dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
            % range (40)
            %
            % scaleMag can be:
            %   ('dB') | 'lin' - only used for magnitude plots
            %
            % scalePhase can be:
            %   ('deg') | 'rad' - only used for phase plots
            %
            % freqUnit can be:
            %   ('GHz') | 'Hz' | 'kHz' | 'MHz' | 'THz'
            %
            % cutConstant can be (used only for polar and cartesian plots):
            %   ('x') | 'y'  (x = [ph|az|ep|u|Xg|asin(u)]; y = [th|el|al|Yg|asin(v)]
            %
            % cutValue can be any value in the available angle range (in rad).
            %
            % step is the plot step size.  Can be empty - then the available data will
            % be used and no surface will be plotted.  If not, a griddata interpolant will be made.
            %
            % plotProperties can be a variety of name, value pairs including:
            %   LineWidth, LineStyle, Color (like '-.')
            %
            % showGrid is a boolean (false) to show the 2D grid where the data is
            % calculated before interpolation.
            %
            % hemisphere is used in gridTypes DirCos and ArcSin and can be:
            %   ('top') | 'bot'
            
            
            
            narginchk(1,40);
            
            %% Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'plot';
            
            typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plot','obj',1);
            addRequired(parseobj,'obj',typeValidationObj);
            
            typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plot','freqIndex');
            addParameter(parseobj,'freqIndex',1,typeValidationFreq);
            
            typeValidationnorm = @(x) validateattributes(x,{'logical','numeric'},{'binary','nonempty','numel',1},'plot','norm');
            addParameter(parseobj,'norm',false,typeValidationnorm );
            
            typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plot','dynamicRange_dB');
            addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );
            
            expectedplotType = {'3D','2D','polar','cartesian'};
            addParameter(parseobj,'plotType','3D', @(x) any(validatestring(x,expectedplotType)));
            
            expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'};
            addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));
            
            expectedoutputType = {'mag','phase','real','imag'};
            addParameter(parseobj,'outputType','mag', @(x) any(validatestring(x,expectedoutputType)));
            
            expectedscaleMag = {'dB','lin'};
            addParameter(parseobj,'scaleMag','dB', @(x) any(validatestring(x,expectedscaleMag)));
            
            expectedscalePhase = {'deg','rad'};
            addParameter(parseobj,'scalePhase','deg', @(x) any(validatestring(x,expectedscalePhase)));
            
            expectedfreqUnit = {'Hz','kHz','MHz','GHz','THz'};
            addParameter(parseobj,'freqUnit','GHz', @(x) any(validatestring(x,expectedfreqUnit)));
            
            expectedcutConstant = {'x','y'};
            addParameter(parseobj,'cutConstant','x', @(x) any(validatestring(x,expectedcutConstant)));
            
            typeValidationcutValue = @(x) validateattributes(x,{'numeric'},{'real'},'plot','cutValue');
            addParameter(parseobj,'cutValue',0,typeValidationcutValue);
            
            typeValidationstep = @(x) validateattributes(x,{'numeric'},{'real'},'plot','step');
            addParameter(parseobj,'step',[],typeValidationstep);     % In degrees
            
            typeValidationLineWidth = @(x) validateattributes(x,{'numeric'},{'real'},'plot','LineWidth');
            addParameter(parseobj,'LineWidth',1,typeValidationLineWidth);
            
            typeValidationLineStyle = @(x) validateattributes(x,{'char'},{'nonempty'},'plot','LineStyle');
            addParameter(parseobj,'LineStyle','-',typeValidationLineStyle);
            
            typeValidationColor = @(x) validateattributes(x,{'numeric','char'},{'nonempty'},'plot','Color');
            addParameter(parseobj,'Color','k',typeValidationColor);
            
            typeValidationshowGrid = @(x) validateattributes(x,{'logical','numeric'},{'binary','nonempty','numel',1},'plot','showGrid');
            addParameter(parseobj,'showGrid',false,typeValidationshowGrid );
            
            expectedhemisphere = {'top','bot'};
            addParameter(parseobj,'hemisphere','top', @(x) any(validatestring(x,expectedhemisphere)));
            
            parse(parseobj, obj, varargin{:});
            
            freqIndex = parseobj.Results.freqIndex;
            output = parseobj.Results.output;
            outputType = parseobj.Results.outputType;
            norm = parseobj.Results.norm;
            dynamicRange_dB = parseobj.Results.dynamicRange_dB;
            plotType = parseobj.Results.plotType;
            scaleMag = parseobj.Results.scaleMag;
            scalePhase = parseobj.Results.scalePhase;
            freqUnitPlot = parseobj.Results.freqUnit;
            cutConstant = parseobj.Results.cutConstant;
            cutValue = parseobj.Results.cutValue;
            step = parseobj.Results.step;
            LineWidth = parseobj.Results.LineWidth;
            LineStyle = parseobj.Results.LineStyle;
            Color = parseobj.Results.Color;
            showGrid = parseobj.Results.showGrid;
            hemisphere = parseobj.Results.hemisphere;
            
            %% Sort out the plot grid and names
            
            % Get valid positions for the plot
            if strcmp(obj.gridType,'DirCos') || strcmp(obj.gridType,'ArcSin')
                % Try to get the direction cosines from the base grid definition - if the
                % base definition is not a direction cosine type it can contain
                % information over the full sphere.
                objBase = obj.grid2Base;
                grid2DirCoshandle = str2func([objBase.gridType,'2DirCos']);
                [~,~,w] = grid2DirCoshandle(objBase.x,objBase.y);
                if strcmp(hemisphere,'top')
                    valAng = w >= 0;
                elseif strcmp(hemisphere,'bot')
                    valAng = w <= 0;
                end
            else
                valAng = ones(obj.Nang,1);
            end
            
            % Get the original grid and output
            X = reshape(obj.x,obj.NyBase,obj.NxBase);
            Y = reshape(obj.y,obj.NyBase,obj.NxBase);
            if strcmp(output,'E1')
                [Z,~,~] = getEfield(obj);
            elseif strcmp(output,'E2')
                [~,Z,~] = getEfield(obj);
            elseif strcmp(output,'E3')
                [~,~,Z] = getEfield(obj);
            else
                outputHandle = str2func(['get',output]);
                Z = outputHandle(obj);
            end
            Z = Z(:,freqIndex);
            
            if isempty(step)
                if strcmp(plotType,'cartesian') || strcmp(plotType,'cartesian')
                    error('Cant do 1D plot with empty step size - this is reserved for making 2D grids etc.')
                else
                    Xi = X;
                    Yi = Y;
                    xi = obj.x;
                    yi = obj.y;
                    NxPlot = obj.Nx;
                    NyPlot = obj.Ny;
                    Zi = Z;
                end
            else
                % Get the interpolated plot points from the step information
                if strcmp(obj.gridType,'DirCos') || strcmp(obj.gridType,'ArcSin')
                    step = sind(step);
                else
                    step = deg2rad(step);
                end
                ximin = min(obj.x);
                ximax = max(obj.x);
                yimin = min(obj.y);
                yimax = max(obj.y);
                NxPlot = round((ximax - ximin)/step) + 1;
                NyPlot = round((yimax - yimin)/step) + 1;
                xivect = linspace(ximin,ximax,NxPlot);
                yivect = linspace(yimin,yimax,NyPlot);
                %     xivect = min(obj.x):step:max(obj.x);
                %     yivect = min(obj.y):step:max(obj.y);
                %     NxPlot = numel(xivect);
                %     NyPlot = numel(yivect);
                [Xi,Yi] = meshgrid(xivect,yivect);
                switch plotType
                    case {'3D','2D'}
                        xi = Xi(:);
                        yi = Yi(:);
                    case {'cartesian','polar'}
                        switch cutConstant
                            case 'x'
                                yi = yivect;
                                xi = ones(size(yi)).*cutValue;
                            case 'y'
                                xi = xivect;
                                yi = ones(size(xi)).*cutValue;
                        end
                end
                [Zi] = interpolateGrid(obj,output,xi,yi,freqIndex,hemisphere);
            end
            
            % Assign axis names
            switch obj.gridType
                case 'DirCos'
                    xiplot = xi;
                    yiplot = yi;
                    xname = [obj.xname, ' = sin(\theta)cos(\phi)'];
                    yname = [obj.yname, ' = sin(\theta)sin(\phi)'];
                    axisUnit = '';
                otherwise
                    X = rad2deg(X);
                    Y = rad2deg(Y);
                    Xi = rad2deg(Xi);
                    Yi = rad2deg(Yi);
                    xiplot = rad2deg(xi);
                    yiplot = rad2deg(yi);
                    axisUnit = '(deg)';
                    xname = [obj.xname, ' ' ,axisUnit];
                    yname = [obj.yname, ' ' ,axisUnit];
            end
            
            %% Condition outputs
            % Phase results
            if (strcmp(output,'E1') || strcmp(output,'E2') || strcmp(output,'E3')) && strcmp(outputType,'phase')
                Zplot = angle(Z);
                Zplot(~valAng) = NaN;
                Ziplot = angle(Zi);
                if norm
                    Zplot = Zplot - max(Zplot);
                    Ziplot = Ziplot - max(Ziplot);
                end
                unit = 'V/m (rad)';
                if strcmp(scalePhase,'deg')
                    Zplot = rad2deg(Zplot);
                    Ziplot = rad2deg(Ziplot);
                    unit = 'V/m (deg)';
                end
                compName = strrep(obj.([output,'name']),'_','\');
            else    % Magnitude/Component results
                switch outputType
                    case 'mag'
                        Zplot = abs(Z);
                        Ziplot = abs(Zi);
                    case 'real'
                        Zplot = real(Z);
                        Ziplot = real(Zi);
                    case 'imag'
                        Zplot = imag(Z);
                        Ziplot = imag(Zi);
                end
                Zplot(~valAng) = NaN;
                switch output
                    case {'Directivity','Gain','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'}
                        dBscale = 10;
                        unit = '';
                        compName = strrep(output,'_','/');
                    case {'E1','E2','E3'}
                        dBscale = 20;
                        unit = 'V/m ';
                        compName = strrep(obj.([output,'name']),'_','\');
                end
                if norm
                    Zplot = Zplot./max(Zplot);
                    Ziplot = Ziplot./max(Ziplot);
                end
                if strcmp(scaleMag,'dB')
                    dBHandle = str2func(['dB',num2str(dBscale)]);
                    Zplot = dBHandle(Zplot);
                    Ziplot = dBHandle(Ziplot);
                    unit = [unit, 'dB'];
                end
            end
            
            %% Make the plots
            switch freqUnitPlot
                case 'Hz'
                    freqMult = 1;
                case 'kHz'
                    freqMult = 1e-3;
                case 'MHz'
                    freqMult = 1e-6;
                case 'GHz'
                    freqMult = 1e-9;
                case 'THz'
                    freqMult = 1e-12;
            end
            freqPlot = obj.freqHz(freqIndex)*freqMult;
            
            switch plotType
                case '3D'
                    % ToDo:
                    % Doesn't work with AzEl/ElAz grids +-90 y-axis breaks the plot
                    
                    % Use the MATLAB antennas toolbox plotting function
                    if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                        % Handle dynamic range here: ToDo
                        if strcmp(outputType,'mag')
                            maxVal = max(Zplot(~isinf(Zplot)));
                            switch scaleMag
                                case 'dB'
                                    if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                        maxVal = (1-norm).*dynamicRange_dB;
                                        minVal = norm.*dynamicRange_dB;
                                    else
                                        minVal = maxVal-dynamicRange_dB;
                                    end
                                    Ziplot(Ziplot<minVal) = minVal;
                                    Ziplot(Ziplot>maxVal) = maxVal;
                                case 'lin'
                                    linHandle = str2func(['lin',num2str(dBscale)]);
                                    if ~norm
                                        dr = linHandle(dynamicRange_dB);
                                        if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                            caxis([0,dr]);
                                        else
                                            caxis([maxVal/dr,maxVal]);
                                        end
                                    end
                            end
                        end
                        iVal = ~isnan(Ziplot);
                        patternCustom(Ziplot(iVal),Yi(iVal),Xi(iVal));
                        title([obj.coorType, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnitPlot])
                    else
                        error(['gridType must be PhTh for 3D plots: found gridType = ', obj.gridType])
                    end
                    
                case '2D'
                    Ziplot = reshape(Ziplot,NyPlot,NxPlot);
                    if ~isempty(step)
                        surf(Xi,Yi,Ziplot,'EdgeColor','Interp','FaceColor','Interp')
                        colorbar
                    end
                    if showGrid || isempty(step)
                        hold on
                        Zplot(valAng == 0) = NaN;
                        plot3(X(:),Y(:),Zplot(:),'k.')
                        hold off
                    end
                    xlabel(xname)
                    ylabel(yname)
                    view([0,90])
                    axis equal
                    xlim([min(xiplot),max(xiplot)])
                    ylim([min(yiplot),max(yiplot)])
                    % Handle dynamic range here
                    if strcmp(outputType,'mag') || strcmp(outputType,'real') || strcmp(outputType,'imag')
                        maxVal = max(Zplot(~isinf(Zplot)));
                        switch scaleMag
                            case 'dB'
                                if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                    rangeZ = [0,dynamicRange_dB] - norm.*dynamicRange_dB;
                                    caxis(rangeZ);
                                    %                         zlim(rangeZ);
                                else
                                    rangeZ = [maxVal-dynamicRange_dB,maxVal];
                                    caxis(rangeZ);
                                    %                         zlim(rangeZ);
                                end
                            case 'lin'
                                linHandle = str2func(['lin',num2str(dBscale)]);
                                if ~norm
                                    dr = linHandle(dynamicRange_dB);
                                    if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                        rangeZ = [0,dr];
                                        caxis(rangeZ);
                                        %                             zlim(rangeZ);
                                    else
                                        rangeZ = [maxVal/dr,maxVal];
                                        if ~(strcmp(outputType,'real') || strcmp(outputType,'imag'))
                                            caxis(rangeZ);
                                        end
                                        %                             zlim([rangeZ]);
                                    end
                                end
                        end
                    end
                    title([obj.coorType, ', ',obj.polType, ' polarisation: ',outputType,'(', compName, ') (',unit,'); Freq = ',num2str(freqPlot),' ', freqUnitPlot])
                case {'cartesian','polar'}
                    % Initial bookkeeping to seperate the two options
                    if strcmp(plotType,'cartesian')
                        plotHandle = str2func(['plot']);
                        limitHandle = str2func(['ylim']);
                        xscale = 1;
                    elseif strcmp(plotType,'polar')
                        assert(~strcmp(obj.gridType,'DirCos'),'Polar plots not supported for DirCos gridType');
                        plotHandle = str2func(['polarplot']);
                        limitHandle = str2func(['rlim']);
                        xscale = pi/180;
                    end
                    switch cutConstant
                        case 'x'
                            plotHandle(yiplot.*xscale,Ziplot,'LineStyle',LineStyle,'LineWidth',LineWidth,'Color',Color), grid on
                            xlab = yname;
                            cutName = obj.xname;
                        case 'y'
                            plotHandle(xiplot.*xscale,Ziplot,'LineStyle',LineStyle,'LineWidth',LineWidth,'Color',Color), grid on
                            xlab = xname;
                            cutName = obj.yname;
                    end
                    % Handle dynamic range here
                    if strcmp(outputType,'mag')
                        maxVal = max(Zplot(:));
                        switch scaleMag
                            case 'dB'
                                if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                    limitHandle([0,dynamicRange_dB] - norm.*dynamicRange_dB);
                                else
                                    limitHandle([maxVal-dynamicRange_dB,maxVal]);
                                end
                            case 'lin'
                                linHandle = str2func(['lin',num2str(dBscale)]);
                                if ~norm
                                    dr = linHandle(dynamicRange_dB);
                                    if strcmp(output,'XP_CO') || strcmp(output,'CO_XP')
                                        limitHandle([0,dr]);
                                    else
                                        limitHandle([maxVal/dr,maxVal]);
                                    end
                                end
                        end
                    end
                    if ~strcmp(obj.gridType,'DirCos')
                        cutValue = rad2deg(cutValue);
                    end
                    titText = [obj.coorType, ', ',obj.polType, ' polarisation; Freq = ',num2str(freqPlot),' ', freqUnitPlot,'; ',cutName, ' = ',num2str(cutValue), ' ',axisUnit];
                    
                    % Final bookkeeping to seperate the two options
                    hold on
                    ax = gca;
                    ylab = [outputType,'(', compName, ') (',unit,')'];
                    if strcmp(plotType,'cartesian')
                        xlabel(xlab)
                        ylabel(ylab)
                        title(titText)
                    elseif strcmp(plotType,'polar')
                        title([ylab,'; ',titText])
                        ax.ThetaZeroLocation = 'top';
                        ax.ThetaDir = 'clockwise';
                    end
            end
            
        end
        
        function plotJones(FF1,FF2,varargin)
            % PLOTJONES Plots the Jones matrix type representation FarFields
            
            % function [] = plotJones(FF1,FF2,varargin)
            % plots the Jones matrix type representation of the farfields specified in
            % FF1 and FF2, where it is assumed they are calculated in the same basis
            % and represent the 2 polarizations
            % name, value are name value pairs, and can be the following:
            %
            % freqIndex is the index of the frequency to be plotted (default 1)
            %
            % dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
            % range (40)
            %
            % step is the plot step size.  Can be empty - then the available data will
            % be used and no surface will be plotted.  If not, a griddata interpolant will be made.
            
            
            %% Parse input
            parseobj = inputParser;
            parseobj.FunctionName = 'plotJones';
            
            typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plotJones','obj',1);
            addRequired(parseobj,'FF1',typeValidationObj);
            addRequired(parseobj,'FF2',typeValidationObj);
            
            typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plotJones','freqIndex');
            addParameter(parseobj,'freqIndex',1,typeValidationFreq);
            
            typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plotJones','dynamicRange_dB');
            addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );
            
            typeValidationstep = @(x) validateattributes(x,{'numeric'},{'real'},'plot','step');
            addParameter(parseobj,'step',[],typeValidationstep);     % In degrees
            
            parse(parseobj, FF1, FF2, varargin{:});
            
            freqIndex = parseobj.Results.freqIndex;
            dynamicRange_dB = parseobj.Results.dynamicRange_dB;
            step = parseobj.Results.step;
            
            
            
            %% Plot the result
            
            if ~isGridEqual(FF1,FF2)
                error('Base grids should be identical for the two input fields');
            else
                subplot(2,2,1)
                plot(FF1,'output','E1','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
                title('J_{11}')
                subplot(2,2,2)
                plot(FF1,'output','E2','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
                title('J_{12}')
                subplot(2,2,3)
                plot(FF2,'output','E1','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
                title('J_{21}')
                subplot(2,2,4)
                plot(FF2,'output','E2','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex)
                title('J_{22}')
            end
        end
        
        function plotPrincipleCuts(FF,varargin)
            % PLOTPRINCIPLECUTS Plots the principle cuts of a FarField in
            % cartesian and polar coordinates.
            
            % [] = plotPrincipleCuts(FF,varargin)
            % plots the principle cuts of the farfield specified in FF
            % name, value are name value pairs, and can be the following:
            %
            % freqIndex is the index of the frequency to be plotted (default 1)
            %
            % dynamicRange_dB is a (positive) dB value for the magnitude plot dynamic
            % range (40)
            %
            % norm is a boolean (false) to normalize to maximum magnitude
            %
            % plotType can be:
            %    ('cartesian') | 'polar'
            %
            % MainPolNr can be
            % (2) | 1 - Indicates which component to plot in solid lines
            %
            % output can be:
            %   ('Directivity') | 'Gain' | 'E1' | 'E2' | 'AxialRatio' | 'AxialRatioInv'
            %   'CO_XP' | 'XP_CO' | 'W' | 'U'
            %   (For Power quantities only one plot is generated, for the rest both will be
            %   plotted, with the first (solid line) determined by the input value)
            
            %% Parse input
            parseobj = inputParser;
            parseobj.FunctionName = 'plotPrincipleCuts';
            
            typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'plotPrincipleCuts','obj',1);
            addRequired(parseobj,'FF',typeValidationObj);
            
            typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'plotPrincipleCuts','freqIndex');
            addParameter(parseobj,'freqIndex',1,typeValidationFreq);
            
            typeValidationDR = @(x) validateattributes(x,{'numeric'},{'real','positive','nonempty','numel',1},'plotPrincipleCuts','dynamicRange_dB');
            addParameter(parseobj,'dynamicRange_dB',40,typeValidationDR );
            
            typeValidationnorm = @(x) validateattributes(x,{'logical','numeric'},{'binary','nonempty','numel',1},'plot','norm');
            addParameter(parseobj,'norm',false,typeValidationnorm );
            
            expectedplotType = {'polar','cartesian'};
            addParameter(parseobj,'plotType','cartesian', @(x) any(validatestring(x,expectedplotType)));
            
            expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','W','U'};
            addParameter(parseobj,'output','Directivity', @(x) any(validatestring(x,expectedoutput)));
            
            parse(parseobj, FF, varargin{:});
            
            freqIndex = parseobj.Results.freqIndex;
            dynamicRange_dB = parseobj.Results.dynamicRange_dB;
            norm = parseobj.Results.norm;
            plotType = parseobj.Results.plotType;
            output = parseobj.Results.output;
            
            
            %% Plot the result
            
            % Estimate a nice step size
            step = median(diff(unique(FF.thBase)));
            if strcmp(FF.gridType,'DirCos') || strcmp(FF.gridType,'ArcSin')
                step = asin(step);
            end
            step = rad2deg(step);
            
            % Output control
            plotSec = true;
            switch output
                case {'Directivity','Gain','W','U'}
                    plotSec = false;
                    Emain = output;
                    ylabText = ['|',output,'| (dB)'];
                case 'E1'
                    Emain = 'E1';
                    Esec = 'E2';
                    ylabText = ['|',FF.E1name,'| (-); |',FF.E2name,'| (--) (dB)' ];
                case 'E2'
                    Emain = 'E2';
                    Esec = 'E1';
                    ylabText = ['|',FF.E2name,'| (-); |',FF.E1name,'| (--) (dB)' ];
                case 'AxialRatio'
                    Emain = 'AxialRatio';
                    plotSec = false;
                    ylabText = ['|AR| (dB)' ];
                case 'AxialRatioInv'
                    Emain = 'AxialRatioInv';
                    plotSec = false;
                    ylabText = ['|AR| (dB)' ];
                case 'CO_XP'
                    Emain = 'CO_XP';
                    plotSec = false;
                    ylabText = ['|CO/XP| (dB)' ];
                case 'XP_CO'
                    Emain = 'XP_CO';
                    plotSec = false;
                    ylabText = ['|XP/CO| (dB)' ];
            end
            
            figure
            switch FF.gridType
                case{'PhTh','AzEl','ElAz'}
                    % Shift the pattern onto a symmetrical grid
                    if ~FF.isGridUniform
                        FF = currentForm2Base(FF,step);
                    else
                        FF = FF.reset2Base;
                    end
                    FF = FF.setXrange('sym');
                    FF = FF.setYrange(360);
                    xVal1 = 0;
                    xVal2 = 90;
                    xVal3 = 45;
                    
                    % Main component
                    plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                        'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal1),...
                        'LineStyle','-','Color','k')
                    plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                        'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal2),...
                        'LineStyle','-','Color','r')
                    plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                        'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal3),...
                        'LineStyle','-','Color','b')
                    
                    if plotSec
                        % 2nd Component
                        plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal1),...
                            'LineStyle','--','Color','k')
                        hold on
                        plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal2),...
                            'LineStyle','--','Color','r')
                        plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal3),...
                            'LineStyle','--','Color','b')
                    end
                    % Add the legend
                    xUnit = '^\circ';
                    legend([FF.xname,'=',num2str(xVal1),xUnit],[FF.xname,'=',num2str(xVal2),xUnit],[FF.xname,'=',num2str(xVal3),xUnit])
                    
                otherwise
                    xVal = 0;
                    yVal = 0;
                    % Main component
                    plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                        'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal),...
                        'cutConstant','x','LineStyle','-','Color','k')
                    plot(FF,'output',Emain,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                        'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(yVal),...
                        'cutConstant','y','LineStyle','-','Color','r')
                    if plotSec
                        % 2nd Component
                        plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(xVal),...
                            'cutConstant','x','LineStyle','--','Color','k')
                        plot(FF,'output',Esec,'outputType','mag','plotType',plotType,'scaleMag','dB','norm',norm,...
                            'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex,'cutValue',deg2rad(yVal),...
                            'cutConstant','y','LineStyle','--','Color','r')
                    end
                    
                    % Add the legend
                    legend([FF.xname,'=',num2str(xVal)],[FF.yname,'=',num2str(yVal)])
            end
            if isequal(plotType,'cartesian')
                ylabel(ylabText)
            end
            
            % Remove the cut value from the title
            h = gca;
            titTextFull = h.Title.String;
            titTextCell = strsplit(titTextFull,'Hz');
            title([titTextCell{1},'Hz'])
            
        end

        function plotGrid(obj,markerStyle)
            % PLOTGRID Plots the FarField object grid.
            
            if nargin < 2
                markerStyle = 'k.';
            end
            switch obj.gridType
                case {'DirCos'}
                    xplot = obj.x;
                    yplot = obj.y;
                    xtext = [obj.xname, ' = sin(\theta)cos(\phi)'];
                    ytext = [obj.yname, ' = sin(\theta)sin(\phi)'];
                    
                otherwise
                    xplot = rad2deg(obj.x);
                    yplot = rad2deg(obj.y);
                    xtext = [obj.xname,' (deg)'];
                    ytext = [obj.yname,' (deg)'];
            end
            plot(xplot,yplot,markerStyle)
            xlabel(xtext)
            ylabel(ytext)
            axis equal
            grid on
            xlim([min(xplot),max(xplot)])
            ylim([min(yplot),max(yplot)])
        end
        
        function plotGridBase(obj)
            % PLOTGRIDBASE Plots the FarField object base grid.
            
            obj = obj.reset2Base;
            obj.plotGrid;
        end
        
        %% Interpolation methods
        function [Zi] = interpolateGrid(obj,output,xi,yi,varargin)
            % INTERPOLATEDGRID Grid interpolation at specified interpolation points.
            
            % Check inputs
            narginchk(4,6);
            
            parseobj = inputParser;
            parseobj.FunctionName = 'interpolateGrid';
            
            typeValidationObj = @(x) validateattributes(x,{'FarField'},{'numel',1},'interpolateGrid','obj',1);
            addRequired(parseobj,'obj',typeValidationObj);
            
            expectedoutput = {'Directivity','Gain','E1','E2','E3','AxialRatio','AxialRatioInv','CO_XP','XP_CO','U','W'};
            addRequired(parseobj,'output', @(x) any(validatestring(x,expectedoutput)));
            
            typeValidationXY = @(x) validateattributes(x,{'numeric'},{'real','nonempty'},'interpolateGrid');
            addRequired(parseobj,'xi',typeValidationXY);
            addRequired(parseobj,'yi',typeValidationXY);
            
            % expectedgridType = {'PhTh','DirCos','AzEl','ElAz','TrueView','ArcSin'};
            % addOptional(parseobj,'gridType',obj.gridType, @(x) any(validatestring(x,expectedgridType)));
            
            typeValidationFreq = @(x) validateattributes(x,{'numeric'},{'real','nonempty','integer'},'interpolateGrid','freqIndex');
            addOptional(parseobj,'freqIndex',1,typeValidationFreq);
            % addParameter(parseobj,'freqIndex',1,typeValidationFreq);
            
            expectedhemisphere = {'top','bot'};
            addOptional(parseobj,'hemisphere','top', @(x) any(validatestring(x,expectedhemisphere)));
            % addParameter(parseobj,'hemisphere','top', @(x) any(validatestring(x,expectedhemisphere)));
            
            parse(parseobj,obj,output,xi,yi,varargin{:});
            
            % gridType = parseobj.Results.gridType;
            freqIndex = parseobj.Results.freqIndex;
            hemisphere = parseobj.Results.hemisphere;
            
            
            % Main code
            
            gridType = obj.gridType;
            
            % Evaluate the field on the base grid - this is where the output function
            % should be best suited for interpolation
            obj = obj.grid2Base;
            % Shift to -180:180 range (if applicable) - this is where the DirCos spits
            % everything out after transforming
            if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                obj = obj.setXrange('sym');
            end
            % Get xi and yi in the base gridType, and on the [-180,180] x-domain for the
            % angular grids
            grid2DirCoshandle = str2func([gridType,'2DirCos']);
            [ui,vi,wi] = grid2DirCoshandle(xi,yi);
            % Check for bottom hemisphere plot - fix the w to be the negative root
            if (strcmp(gridType,'DirCos') || strcmp(gridType,'ArcSin')) && strcmp(hemisphere,'bot')
                wi = -wi;
            end
            DirCos2baseHandle = str2func(['DirCos2',obj.gridType]);
            [xi_bGT,yi_bGT] = DirCos2baseHandle(ui,vi,wi);
            % Find the invalid points included by the external meshgrid
            valAngi = sqrt(ui.^2 + vi.^2) <= 1;
            % Sort out the TrueView special case invalid points
            if strcmp(gridType,'TrueView')
                valAngi = sqrt((xi./pi).^2 + (yi./pi).^2) <= 1;
            end
            
            % Get the valid angle positions - already in baseGrid here, but shifted to
            % +- 180 degrees
            % Get the indexes only, no later reshaping done, different from the grids
            % required for plotting in plot.m
            if strcmp(gridType,'DirCos') || strcmp(gridType,'ArcSin')
                grid2DirCoshandleBase = str2func([obj.gridType,'2DirCos']);
                [~,~,w] = grid2DirCoshandleBase(obj.x,obj.y);
                if strcmp(hemisphere,'top')
                    valAng = find(w >= 0);
                elseif strcmp(hemisphere,'bot')
                    valAng = find(w <= 0);
                end
            else
                valAng = 1:obj.Nang;
            end
            valAng = valAng(:);
            
            % Extract the outputs on the base grid
            if strcmp(output,'E1')
                [Zfreq,~,~] = getEfield(obj);
            elseif strcmp(output,'E2')
                [~,Zfreq,~] = getEfield(obj);
            elseif strcmp(output,'E3')
                [~,~,Zfreq] = getEfield(obj);
            else
                outputHandle = str2func(['get',output]);
                Zfreq = outputHandle(obj);
            end
            
            % Select frequency of interest
            Z = Zfreq(:,freqIndex);
            
            xVal = obj.x(valAng);
            yVal = obj.y(valAng);
            zVal = Z(valAng);
            
            edgeAngExtent_deg = 16;
            % Extend grid past -180 and +180 for interpolation across the axis
            if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                tol = deg2rad(2); % Check for points close to the edges in azimuth
                if abs(min(xVal) + pi) < tol && abs(max(xVal) - pi) < tol
                    edgeAngDeg = 180 - edgeAngExtent_deg;
                    iNeg = find(xVal > deg2rad(edgeAngDeg));
                    iPos = find(xVal < deg2rad(-edgeAngDeg));
                    xVal = [xVal(iNeg)-2*pi;xVal;xVal(iPos)+2*pi];
                    yVal = [yVal(iNeg);yVal;yVal(iPos)];
                    zVal = [zVal(iNeg);zVal;zVal(iPos)];
                end
                %     % Also extend the y-axis
                %     if abs(min(yVal) - 0) < tol
                %         edgeAngDeg = edgeAngExtent_deg;
                %         iNeg = find(xVal < 0 & yVal < deg2rad(edgeAngDeg));
                %         iPos = find(xVal > 0 & yVal < deg2rad(edgeAngDeg));
                %         xVal = [xVal;xVal(iNeg)+pi;xVal(iPos)-pi];
                %         yVal = [yVal;-yVal(iNeg);-yVal(iPos)];
                %         zVal = [zVal;zVal(iNeg);zVal(iPos)];
                %     end
                %     if abs(max(yVal) - pi) < tol
                %         edgeAngDeg = 180 - edgeAngExtent_deg;
                %         iNeg = find(xVal < 0 & yVal > deg2rad(edgeAngDeg));
                %         iPos = find(xVal > 0 & yVal > deg2rad(edgeAngDeg));
                %         xVal = [xVal;xVal(iNeg)+pi;xVal(iPos)-pi];
                %         yVal = [yVal;2*pi-yVal(iNeg);2*pi-yVal(iPos)];
                %         zVal = [zVal;zVal(iNeg);zVal(iPos)];
                %     end
                
            end
            
            % Remove duplicate differing values completely from the set - interpolate
            % over them.  This happens at poles for certain coordinate projections
            % (like at th = 180 in Ludwig3 for instance, of th=0 and 180 for spherical)
            % First find duplicate domain values
            [~,iUnique] = unique([xVal,yVal],'rows');
            removePoints = [];
            if length(iUnique) < length(xVal)
                iRepeated = setdiff((1:length(xVal)).',iUnique);
                repeatedSet = [xVal(iRepeated),yVal(iRepeated),zVal(iRepeated)];
                % Find the repeated values
                repVals = unique(repeatedSet(:,1:2),'rows');
                % Test if different z-values occur
                for ii = 1:length(repVals(:,1))
                    currentValRowIndex = find(ismember(repeatedSet(:,1:2),repVals(ii,:),'rows'));
                    currentZ = repeatedSet(currentValRowIndex,3);
                    if ~all(currentZ == currentZ(1))
                        removePoints = [removePoints;repeatedSet(currentValRowIndex,1:2)];
                    end
                end
            end
            if numel(removePoints) > 0
                iRemove = find(ismember([xVal,yVal],removePoints,'rows'));
                xVal(iRemove) = [];
                yVal(iRemove) = [];
                zVal(iRemove) = [];
            end
            
            % Build the interpolant on the base grid at the valid angles
            if obj.isGridUniform
                try
                    NyVal = length(unique(yVal));
                    NxVal = length(unique(xVal));
                    XVal = reshape(xVal,NyVal,NxVal);
                    YVal = reshape(yVal,NyVal,NxVal);
                    ZVal = reshape(zVal,NyVal,NxVal);
                    Zf = griddedInterpolant(XVal',YVal',ZVal.','linear');
                catch ME
                    % Grid did not work... Go to scatter
                end
            end
            if ~exist('Zf','var')
                warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
                Zf = scatteredInterpolant(xVal,yVal,zVal,'linear');
                warning('on','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
            end
            % Get the values on the scattered set of inputs
            Zi = Zf(xi_bGT,yi_bGT);
            Zi(~valAngi) = NaN;
        end

        %% Phase centre/shifts/rotations of the field
        function [Z, Delta, delta0, eta_pd] = phaseCentreKildal(FF,pol,th_M)
            % PHASECENTREKILDAL Computes the phase centre and approximates
            % the phase efficiency.
            
            % function [PC, eta_phi] = PCoptFunc(pol,freq_vect,th_M,BOR1pattern)
            % computes the phase center and approximate phase efficiency of a given
            % farfield pattern using Kildal 1983 paper
            % Returns:
            % Z - the total phase center (Delta + delta0) position in [m] as a frequency vector
            % Delta - the initial approximation of the PC [m]
            % delta0 - the fine approximation around Delta of the PC [m]
            % eta_pd - the phase efficiency pu as a frequency vector (NaN for large PC displacements)
            % Inputs:
            % FF - FarField class
            % pol - {'x' | 'y' | 'lh' | 'rh'}
            % th_M - subtended angle of the 'reflector' in rad
            
            assert(strcmp(pol,'x')||strcmp(pol,'y')||strcmp(pol,'lh')||strcmp(pol,'rh')||isa('pol','double'),['Error: Unknown parameter for pol: ',pol])
            % Get in BOR1
            if ~FF.symmetryBOR1
                FF = FF.getBOR1pattern;
            end
            freq_vect = FF.freq;
            lambda_vect = FF.c0./freq_vect;
            k_vect = 2*pi./lambda_vect;
            
            th = FF.th(1:FF.Ny);
            
            [A1f,B1f,C1f,D1f] = FF.getBOR1comps;
            
            [phi_0,phi_th_M,k_Delta,Delta,k_delta0,delta0,Z,eta_pd] = deal(zeros(1,FF.Nf));
            for ff = 1:FF.Nf
                
                switch pol
                    case 'x'
                        CO = B1f(:,ff) + D1f(:,ff);
                    case {'y','lh','rh'}
                        CO = A1f(:,ff) + C1f(:,ff);
                    otherwise
                        error(['Unknown pol: ', pol])
                end
                
                % Move the pattern to the approximate PC (Kildal comments 1984)
                phi = unwrap(angle(CO));
                phi_0(ff) = phi(1);
                phi_th_M(ff) = interp1(th,phi,th_M);
                k_Delta(ff) = (phi_0(ff) - phi_th_M(ff))/(1 - cos(th_M));
                phi_Delta = phi - k_Delta(ff).*cos(th);
                Delta(ff) = k_Delta(ff)./k_vect(ff);
                
                % Find the PC from the formulas in Kildal 1983 (maximum eff method)
                % Weighting function
                w = abs(CO).*tan(th./2);
                % Integral constants (change number of points for very sharp patterns...)
                th_int = linspace(0,th_M,501);
                w_int = interp1(th,w,th_int);
                phi_Delta_int = interp1(th,phi_Delta,th_int);
                phi_Delta0 = interp1(th,phi_Delta,0);
                
                Iw = integral1D(th_int,w_int);
                Iwp = integral1D(th_int,w_int.*(phi_Delta_int - phi_Delta0));
                Iwp2 = integral1D(th_int,w_int.*(phi_Delta_int - phi_Delta0).^2);
                Iwc = integral1D(th_int,w_int.*(cos(th_int) - 1));
                Iwc2 = integral1D(th_int,w_int.*(cos(th_int) - 1).^2);
                Iwpc = integral1D(th_int,w_int.*(cos(th_int) - 1).*(phi_Delta_int - phi_Delta0));
                
                % PC pos
                k_delta0(ff) = (Iw.*Iwpc - Iwp.*Iwc)./(Iwc2.*Iw - Iwc.^2);
                delta0(ff) = k_delta0(ff)./k_vect(ff);
                
                Z(ff) = Delta(ff) + delta0(ff);
                
                % Calculate the approximate phase efficiency
                if abs(Z(ff)/lambda_vect(ff)) > pi/4   % Return NaN for large PC errors - approximation not valid
                    eta_pd(ff) = NaN;
                else
                    c = 1 - Iwp2/Iw + (Iwp/Iw)^2;
                    b = Iwpc/Iw - (Iwp*Iwc/Iw^2);
                    a = Iwc2/Iw - (Iwc/Iw)^2;
                    
                    eta_pd(ff) = c + 2*b*k_vect(ff)*Z(ff) - a*(k_vect(ff)*Z(ff)).^2;
                end
            end
        end
        
        function obj = rotate(obj1,rotHandle,rotAng,onlyRotPowerPattern)
            % ROTATE Rotation function for FarField objects.
            
            % General rotation function for FarField objects
            % rotHandle is the function handle for the type of rotation:
            %   rotx3Dsph, roty3Dsph, rotz3Dsph, rotGRASPsph, rotEulersph
            % rotAng is the associated angle in rad. Scalar for rotations
            % around an axis, and [3x1] for GRASP or Euler rotations
            %
            % onlyRotPowerPattern is an optional argument which speeds up
            % the method in the case where only the rotated power pattern
            % is of interest.  The field values will be arbitrary, bu the
            % power pattern (directivity etc.) will be correct.  Used often
            % for noise temeperature calculations.
            
            if nargin < 4
                onlyRotPowerPattern = false;
            end
            
            % This will probably depent on the pattern type which one is
            % best.  Only have spherical and Ludwig 3 implemented for now,
            % so hard-coded.
            baseCoorType = 'Ludwig3';
            coorHandle = str2func(['coor2',baseCoorType]);
            
            % Test if the rotation function handle has the trailing 'sph'
            handleStr = func2str(rotHandle);
            if ~strcmp(handleStr(end-2:end),'sph')
                handleStr = [handleStr,'sph'];  % Add it if not - some user errors fixed at least!
                rotHandle = str2func(handleStr);
            end
            % Transform to sensible grid and coordinate system for rotation
            FFsph = obj1.grid2PhTh;  % Always work in the PhTh coordinate system
            FFsph = FFsph.setXrange('sym'); % Always work in symmetrical xRange
            FFsph = coorHandle(FFsph,false);
            
            %             % Force all the th = 0|180 fields to be identical - fixes pole
            %             % interpolation problems
            %             tol = 10^(-obj1.nSigDig);
            %             i_0_0 = find(abs(FFsph.th) < tol & abs(FFsph.ph) < tol);
            %             i_180_0 = find(abs(FFsph.th - pi) < tol & abs(FFsph.ph) < tol);
            %             if ~isempty(i_0_0)
            %                 E1_0_0 = obj1.E1(i_0_0(1),:);
            %                 E2_0_0 = obj1.E2(i_0_0(1),:);
            %                 i_0 = find(abs(FFsph.th - 0)<tol);
            %                 FFsph.E1(i_0,:) = repmat(E1_0_0,length(i_0),1);
            %                 FFsph.E2(i_0,:) = repmat(E2_0_0,length(i_0),1);
            %             end
            %             if ~isempty(i_180_0)
            %                 E1_180_0 = obj1.E1(i_180_0(1),:);
            %                 E2_180_0 = obj1.E2(i_180_0(1),:);
            %                 i_180 = find(abs(FFsph.th - pi)<tol);
            %                 FFsph.E1(i_180,:) = repmat(E1_180_0,length(i_180),1);
            %                 FFsph.E2(i_180,:) = repmat(E2_180_0,length(i_180),1);
            %             end
            
            % Get the grid step sizes from the original
            %             stepx = (max(FFsph.x) - min(FFsph.x))./(FFsph.Nx-1);
            %             stepy = (max(FFsph.y) - min(FFsph.y))./(FFsph.Ny-1);
            stepx = (max(FFsph.ph) - min(FFsph.ph))./(FFsph.Nx-1);
            stepy = (max(FFsph.th) - min(FFsph.th))./(FFsph.Ny-1);
            stepDeg = rad2deg([stepx,stepy]);
            xmin = min(FFsph.x);
            xmax = max(FFsph.x);
            ymin = min(FFsph.y);
            ymax = max(FFsph.y);
            % Perform the rotation of the grid
            phIn = FFsph.x.';
            thIn = FFsph.y.';
            sphAngIn = [phIn;thIn];
            sphAngRot = rotHandle(sphAngIn,rotAng);
            phOut = sphAngRot(1,:).';
            thOut = sphAngRot(2,:).';
            FFsph.x = phOut;
            FFsph.y = thOut;
            
            % Perform the rotation of the field vectors
            % Coordinate systems
            C0 = CoordinateSystem;
            rotHandleStr = func2str(rotHandle);
            rotHandleCoor = str2func(rotHandleStr(1:end-3));
            Crot = rotHandleCoor(C0,rotAng);
            
            % Vector origin points before rotation
            [OIx,OIy,OIz] = PhTh2DirCos(phIn,thIn);
            origin_In = Pnt3D(OIx(:).',OIy(:).',OIz(:).');
            origin_out = origin_In.changeBase(C0,Crot);
            
            % local unit vector directions before and after rotation
            switch baseCoorType
                case 'spherical'
                    [xHatIn,yHatIn,~] = unitVectorsSpherical(thIn,phIn);
                    [xHatOut,yHatOut,~] = unitVectorsSpherical(thOut,phOut);
                    [E1sign,E2sign] = deal(1,1);
                case 'Ludwig3'
                    [xHatIn,yHatIn] = unitVectorsDirCos(thIn,phIn);
                    [xHatOut,yHatOut] = unitVectorsDirCos(thOut,phOut);
                    [E1sign,E2sign] = deal(-1,-1);
            end
            % Vector tip points before rotation
            xTipIn = origin_In + xHatIn;
            yTipIn = origin_In + yHatIn;
            % Rotate all the points
            xTipOut = xTipIn.changeBase(C0,Crot);
            yTipOut = yTipIn.changeBase(C0,Crot);
            xOut = xTipOut - origin_out;
            yOut = yTipOut - origin_out;
            % Project onto the local unit vectors
            xx = dot(xOut.pointMatrix,xHatOut.pointMatrix);
            xy = dot(xOut.pointMatrix,yHatOut.pointMatrix);
            yx = dot(yOut.pointMatrix,xHatOut.pointMatrix);
            yy = dot(yOut.pointMatrix,yHatOut.pointMatrix);
            E1rot = FFsph.E1.*repmat(xx(:),1,FFsph.Nf) + FFsph.E2.*repmat(yx(:),1,FFsph.Nf);
            E2rot = FFsph.E1.*repmat(xy(:),1,FFsph.Nf) + FFsph.E2.*repmat(yy(:),1,FFsph.Nf);
            FFsph.E1 = E1sign.*E1rot;
            FFsph.E2 = E2sign.*E2rot;
            
            % Set the baseGrid of the rotated object.  This is required
            % since all transformations operate from the base grid
            FFsph = FFsph.sortGrid;
            FFsph = FFsph.setBase;
            if onlyRotPowerPattern
                % Build the new grid
                Nxi = round((xmax - xmin)/stepx) + 1;
                Nyi = round((ymax - ymin)/stepy) + 1;
                xivect = linspace(xmin,xmax,Nxi);
                yivect = linspace(ymin,ymax,Nyi);
                [Xi,Yi] = meshgrid(xivect,yivect);
                xi = Xi(:);
                yi = Yi(:);
                % Interpolate the fields
                [Ugrid] = deal(zeros(Nxi*Nyi,obj1.Nf));
                for ff = 1:obj1.Nf
                    Ugrid(:,ff) = interpolateGrid(FFsph,'U',xi,yi,ff,'top');
                end
                FFsph = FarField.farFieldFromPowerPattern(xi,yi,Ugrid,FFsph.freq,'fieldPol','linearY','freqUnit',FFsph.freqUnit);
            else
                FFsph = FFsph.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
            end
            % Reset the grid and coordinate system, and reset the base back
            % in the original format
            obj = transformTypes(FFsph, obj1);
            % Recall symmetries - only set the symmetries if they are still
            % valid after rotation
            tol = 10^(-obj.nSigDig);
            YZsym = abs((abs(dot(C0.x_axis,Crot.x_axis)) - 1)) < tol;
            XZsym = abs((abs(dot(C0.y_axis,Crot.y_axis)) - 1)) < tol;
            XYsym = abs((abs(dot(C0.z_axis,Crot.z_axis)) - 1)) < tol;
            if YZsym, obj = obj.setSymmetryYZ(obj1.symmetryYZ); end
            if XZsym, obj = obj.setSymmetryXZ(obj1.symmetryXZ); end
            if XYsym, obj = obj.setSymmetryXY(obj1.symmetryXY); end
            obj = obj.setXrange(obj1.xRangeType);
            obj = obj.currentForm2Base();
        end
        
        function obj = shift(obj,shiftVect)
            % SHIFT Shifts the FarField by a specified distance.
            
            % Shifts the FarField by a distance specified in the Pnt3D input
            % shiftVect (only uses the first entry)
            % shiftVect can also be a vector of length 3 with elements
            % [delX,delY,delZ] in m
            if nargin < 2, shiftVect = Pnt3D; end
            if ~isa(shiftVect,'Pnt3D')
                assert(numel(shiftVect)==3,'Expect a vector of length 3 for the shifVect, if not a Pnt3D');
                shiftVect = Pnt3D(shiftVect(1),shiftVect(2),shiftVect(3));

            end
            assert(numel(shiftVect.x)==1,'Only one point allowed when shifting a FarField');
            
            lam = obj.c0./obj.freqHz;
            k = 2.*pi./lam;
            
            r_hat = [sin(obj.th).*cos(obj.ph), sin(obj.th).*sin(obj.ph), cos(obj.th)];
            rmat = repmat(shiftVect.pointMatrix.',obj.Nang,1);
            rdotr = dot(rmat,r_hat,2);
            phase = exp(1i.*bsxfun(@times,k,rdotr));
            obj.E1 = obj.E1.*phase;
            obj.E2 = obj.E2.*phase;
            obj = obj.setBase;
        end
        
        
        %% Maths - all overloaded methods
        function obj = plus(obj1,obj2)
            % PLUS Add two FarFields
            
            obj1base = reset2Base(obj1);
            obj2base = reset2Base(obj2);
            
            if isGridEqual(obj1base,obj2base) && typesAreEqual(obj1base,obj2base)
                obj = obj1base;
                obj.E1 = obj1base.E1 + obj2base.E1;
                obj.E2 = obj1base.E2 + obj2base.E2;
                obj.Prad = obj1base.Prad + obj2base.Prad;
                Pt = obj1base.Prad./obj1base.radEff + obj2base.Prad./obj2base.radEff;
                obj.radEff = obj.Prad./Pt;
                obj = setBase(obj);
            else
                error('Can only add FarFields with equal base grids')
            end
            
            if typesAreEqual(obj1,obj2)
                obj = transformTypes(obj, obj1);
            end
        end
        
        function obj = minus(obj1,obj2)
            % MINUS Subtract two FarFields
            
            obj1base = reset2Base(obj1);
            obj2base = reset2Base(obj2);
            
            if isGridEqual(obj1base,obj2base) && typesAreEqual(obj1base,obj2base)
                obj = obj1base;
                obj.E1 = obj1base.E1 - obj2base.E1;
                obj.E2 = obj1base.E2 - obj2base.E2;
                obj.Prad = obj1base.Prad + obj2base.Prad;
                Pt = obj1base.Prad./obj1base.radEff + obj2base.Prad./obj2base.radEff;
                obj.radEff = obj.Prad./Pt;
                obj = setBase(obj);
            else
                error('Can only subtract FarFields with equal base grids')
            end
            
            if typesAreEqual(obj1,obj2)
                obj = transformTypes(obj, obj1);
            end
        end
        
        function obj = times(obj1,obj2)
            % TIMES Multiply two FarFields
            
            % Don't operate in BaseGrid - straight on the actual values.
            % This is often used with conj, which operates in the current
            % grid.
            if isGridEqual(obj1,obj2) && typesAreEqual(obj1,obj2)
                obj = obj1;
                obj.E1 = obj1.E1.*obj2.E1;
                obj.E2 = obj1.E2.*obj2.E2;
                obj.Prad = obj.pradInt;
                obj.radEff = ones(size(obj.Prad));
                obj = setBase(obj);
            else
                error('Can only multiply FarFields with equal grids')
            end
        end
        
        function obj = conj(obj1)
            % CONJ Get the complex conjugate of E-field values in a FarField.
            
            obj = obj1;
            obj.E1 = conj(obj1.E1);
            obj.E2 = conj(obj1.E2);
        end
        
        function obj = abs(obj1)
            % ABS Get the absolute value of E-field components in a FarField.
            
            obj = obj1;
            obj.E1 = abs(obj1.E1);
            obj.E2 = abs(obj1.E2);
        end
        
        function obj = scale(obj1,scaleFactor)
            % SCALE Scale E-field components by a scaleFactor in a FarField.
            
            % Scale the FarField object E-fields by the scaleFactor
            obj = obj1;
            obj.E1 = obj1.E1.*scaleFactor;
            obj.E2 = obj1.E2.*scaleFactor;
            obj.Prad = obj1.Prad.*(abs(scaleFactor).^2);
        end
        
        function [normE] = norm(obj,Ntype)
            % NORM Calculate the vector norm of the three E-field components.
            
            % Calculate the vector norm of the three E-field components - overloads the MATLAB norm function
            % This is used for error checking mostly - when comparing
            % different fields for instance...
            if nargin == 1
                Ntype = 2;
            end
            nE1 = norm(obj.E1,Ntype);
            nE2 = norm(obj.E2,Ntype);
            nE3 = norm(obj.E3,Ntype);
            normE = [nE1,nE2,nE3];
        end
        
        function [rmsE1,rmsE2,rmsE3] = rms(obj,DIM)
            % RMS Calculate the rms of E-field components in a FarField.
            
            % Calculate the rms of the three E-field components - overloads the MATLAB rms function
            % This is used for error checking mostly - when comparing
            % different fields for instance...
            if nargin < 2
                DIM = 1;
            end
            rmsE1 = rms(obj.E1,DIM);
            rmsE2 = rms(obj.E2,DIM);
            rmsE3 = rms(obj.E3,DIM);
        end
        
        function T = convPower(obj1,obj2)
            % CONVPOWER Convolve the power patterns, over the full sphere,
            % of two FarField objects.
            
            % Convolve the power patterns, over the full sphere, of two
            % FarField objects. Typically used for antenna temperature
            % calculations
            obj1 = reset2Base(obj1);
            obj2 = reset2Base(obj2);
            
            if isGridEqual(obj1,obj2)
                P = obj1.getU.*obj2.getU;
                FF_T = FarField.farFieldFromPowerPattern(obj1.phBase,obj1.thBase,P,obj1.freq);
                T = FF_T.pradInt;
            else
                error('Can only convolve FarFields with equal base grids')
            end
        end
        
        %% Field normalization
        function P = pradInt(obj)
            % PRADINT  Calculates the total power in the field integrated
            % over the full available grid.
            
            % Returns the total power in the field integrated over the
            % full available grid
            obj = reset2Base(obj);
            symFact = 2^(sum(abs([obj.symXY,obj.symXZ,obj.symYZ])));
            assert(obj.isGridUniform,'Must have a plaid, monotonic, uniform grid for power calculation through integration');
            
            
            switch obj.gridType
                case {'PhTh','AzAlt','RAdec','GalLongLat'}
                    if strcmp(obj.gridType,'PhTh')
                        JacFunc = @sin;
                    else
                        JacFunc = @cos;
                    end
                    PH = reshape(obj.x,obj.Ny,obj.Nx);
                    TH = reshape(obj.y,obj.Ny,obj.Nx);
                    U = obj.getU;
                    P = zeros(1,obj.Nf);
                    for ff = 1:obj.Nf
                        if strcmp(obj.symmetryBOR,'none')
                            integrand = reshape(U(:,ff),obj.Ny,obj.Nx).*JacFunc(TH);
                            P(ff) = integral2D(PH,TH,integrand);
                        else
                            Nth = obj.Ny;
                            th_vect = obj.y(1:Nth);
                            if strcmp(obj.symmetryBOR,'BOR0')
                                integrand = 2*U(:,ff).*sin(th_vect);
                            elseif strcmp(obj.symmetryBOR,'BOR1')
                                integrand = (U(1:Nth,ff) + U(Nth+1:end,ff)).*JacFunc(th_vect);
                            end
                            P(ff) = pi*integral1D(th_vect,integrand);
                            symFact = 1;    % Just to be sure...
                        end
                    end
                otherwise
                    error(['pradInt not implemented for gridType = ',obj.gridType,', only for PhTh and astronomical grids'])
            end
            P = P.*symFact;
        end
        
        function obj = setPower(obj1,powerWatt)
            % SETPOWER Normalizes the FarField object to have the total
            % radiated power specified in powerWatt (in Watts)
            
            % Normalizes the FarField object obj1 to have the a total
            % radiated power specified in powerWatt (in Watts, of course)
            % powerWatt can be a vector of length obj.Nf or scalar.
            % The field need not be specified over the full sphere - the
            % total intercepted power in the specified sector will be set
            % to powerWatt.  Default value of 4*pi W will be used for one
            % argument.
            % The grid should be the standard plaid, montonic, uniform grid.
            obj1 = reset2Base(obj1);
            if nargin == 1
                powerWatt = 4*pi;
            end
            if length(powerWatt) == 1
                powerWatt = repmat(powerWatt,1,obj1.Nf);
            elseif length(powerWatt) ~= obj1.Nf
                error('powerWatt should be scalar or of length obj1.Nf');
            end
            P = obj1.Prad;
            Cn = powerWatt./(P);
            obj = obj1;
            obj = scale(obj,sqrt(Cn));
        end
        
        %% Frequency modifications
        function obj = getFi(obj1,freqIndex)
            % GETFI Returns an object only containing the results in
            % freqIndex.
            
            % Returns an object only containing the results in freqIndex
            obj = obj1;
            obj.E1 = obj1.E1(:,freqIndex);
            obj.E2 = obj1.E2(:,freqIndex);
            obj.freq = obj1.freq(freqIndex);
            obj.Prad = obj1.Prad(freqIndex);
            obj.radEff = obj1.radEff(freqIndex);
            obj = setBase(obj);
        end
        
        %% Symmetry handlers
        function obj = setSymmetryXZ(obj,symmetryType)
            % SETSYMMETRYXZ Specify symmetry along the XZ plane.
            
            mustBeMember(symmetryType,{'none','electric','magnetic'})
            if ~strcmp(symmetryType,'none')
                % Test if the input range is valid
                tol = 10^(-obj.nSigDig);
                % Easy to check in TrueView
                obj1 = obj.grid2TrueView;
                assert(all(sign(obj1.y+tol) > 0) || all(sign(obj1.y-tol) < 0),'Invalid range for XZ symmetry')
            end
            obj.symmetryXZ = symmetryType;
        end
        
        function obj = setSymmetryYZ(obj,symmetryType)
            % SETSYMMETRYYZ Specify symmetry along the YZ plane.
            
            mustBeMember(symmetryType,{'none','electric','magnetic'})
            if ~strcmp(symmetryType,'none')
                % Test if the input range is valid
                tol = 10^(-obj.nSigDig);
                % Easy to check in TrueView
                obj1 = obj.grid2TrueView;
                assert(all(sign(obj1.x+tol) > 0) || all(sign(obj1.x-tol) < 0),'Invalid range for YZ symmetry')
            end
            obj.symmetryYZ = symmetryType;
        end
        
        function obj = setSymmetryXY(obj,symmetryType)
            % SETSYMMETRYXY Specify symmetry along the XY plane.
            
            mustBeMember(symmetryType,{'none','electric','magnetic'})
            %             warning('function: setSymmetryXY not implemented yet - unchanged object returned');
        end
        
        function obj = mirrorSymmetricPattern(obj1)
            % MIRRORSYMMETRICPATTERN Returns the full pattern mirrored
            % according to the symmetry definitions
            
            % Returns the full pattern mirrored according to the symmetry
            % definitions
            
            if ~obj1.symXZ && ~obj1.symYZ && ~obj1.symXZ
                obj = obj1;
            else
                gridTypeIn = obj1.gridType;
                coorTypeIn = obj1.coorType;
                if strcmp(obj1.gridTypeBase,'DirCos') || strcmp(obj1.gridType,'ArcSin')
                    stepX = asin(min(abs(diff(unique(obj1.xBase)))));
                    stepY = asin(min(abs(diff(unique(obj1.yBase)))));
                else
                    % Sort out rounding errors for degrees
                    stepX = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.xBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                    stepY = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.yBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                end
                
                gridHandle = str2func(['grid2',gridTypeIn]);
                coorHandle = str2func(['coor2',coorTypeIn]);
                obj1 = obj1.grid2TrueView;
                obj1 = obj1.coor2Ludwig3(false);   % Always work in H/V for symmetry calculations...
                
                % Initialise
                XIn = [obj1.x];
                YIn = [obj1.y];
                E1In = [obj1.E1];
                E2In = [obj1.E2];
                if obj1.symXZ
                    XIn = [XIn;XIn];
                    YIn = [YIn;-YIn];
                    E1In = [E1In;obj1.symXZ.*E1In]; % Mirror according to symmetry
                    E2In = [E2In;-obj1.symXZ.*E2In];  % Mirror according to symmetry
                end
                if obj1.symYZ
                    XIn = [XIn;-XIn];
                    YIn = [YIn;YIn];
                    E1In = [E1In;-obj1.symYZ.*E1In];  % Mirror according to symmetry
                    E2In = [E2In;obj1.symYZ.*E2In]; % Mirror according to symmetry
                end
                % Object for grid/coor transformation
                objD = FarField(XIn,YIn,E1In,E2In,[],...
                    obj1.freq,obj1.Prad.*2,obj1.radEff,obj1.coorType,obj1.polType,obj1.gridType,obj1.freqUnit,obj1.slant);
                obj = coorHandle(objD);
                obj = gridHandle(obj);
                obj = obj.currentForm2Base(rad2deg([stepX,stepY]));
                % Test here for full 4pi grid, and if not, add the missing
                % axis/pole
                if obj1.isGrid4pi && ~obj.isGrid4pi
                    % Initially test for for the azimuth angles in the angle
                    % specified grids
                    if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                        % We are expecting a symmetric x-grid from the TrueView
                        % transformation earlier, so only check for {-180,180}
                        xmin = min(obj.x);
                        xmax = max(obj.x);
                        tol = mean(diff(unique(obj.x)));
                        if abs(xmin + pi) > tol/10
                            % Insert a -pi cut from pi
                            obj = copyAndInsertXcut(obj,pi,-pi,tol);
                        end
                        if abs(xmax - pi) > tol/10
                            % Insert a pi cut from -pi
                            obj = copyAndInsertXcut(obj,-pi,pi,tol);
                        end
                    else
                        warning(['A full 4pi grid is expected, but not achieved.  Fix not yet implemented for gridType: ',obj.gridType]);
                    end
                end
            end
        end
        
        function obj = getBORpattern(obj1,BORcomp)
            % GETBORPATTERN expands the input FarField pattern into its BOR
            % components, BOR0 and BOR1. 
            
            % Function that expands the input FarField pattern into its BOR
            % components, and returns a FarField object which only contains
            % the BOR0 or BOR1 components.  The output field has the same th
            % samples as the input field, but only the principle ph cuts
            
            if nargin == 1
                BORcomp = 1;
            else
                mustBeMember(BORcomp,[0,1]);
            end
            
            tol = 10^(-obj1.nSigDig+1);
            assert(strcmp(obj1.gridType,'PhTh'),'getBORpattern only operates on PhTh grid patterns');
            assert(abs(max(obj1.x) - min(obj1.x)) - 2*pi < tol,'The ph cuts must span 2*pi for BOR expansion');
            assert(obj1.isGridUniform,'A plaid, monotonic, uniform grid is expected for BOR expansion');
            assert(strcmp(obj1.coorType,'spherical'),'getBORpattern only operates on spherical coorType');
            
            Nph = obj1.Nx;
            Nth = obj1.Ny;
            th_vect = unique(obj1.y);
            ph_vect = unique(obj1.x);
            
            % Calculate the DFT in ph to get the BOR components
            % STore th variation in columns and BOR components row-wise
            [An,Bn,Cn,Dn] = deal(zeros(floor(Nph - 1)/2+1,Nth));
            [Ath,Bth,Cth,Dth] = deal(zeros(Nth,obj1.Nf));
            [BORpower] = deal(zeros(1,obj1.Nf));
            for ff = 1:obj1.Nf
                %                 for nn = 0:floor((Nph - 1)/2)
                for nn = 0:BORcomp % Just get what is required for speed - can slot the rest in if more modes are needed later
                    sin_vect = sin(nn*ph_vect);
                    cos_vect = cos(nn*ph_vect);
                    for tt = 1:Nth
                        Gth_vect = obj1.E1((0:(Nph-1))*Nth+tt,ff);
                        Gph_vect = obj1.E2((0:(Nph-1))*Nth+tt,ff);
                        An(nn+1,tt) = 2/Nph.*sum(Gth_vect(:).*sin_vect(:));
                        Bn(nn+1,tt) = 2/Nph.*sum(Gth_vect(:).*cos_vect(:));
                        Cn(nn+1,tt) = 2/Nph.*sum(Gph_vect(:).*cos_vect(:));
                        Dn(nn+1,tt) = -2/Nph.*sum(Gph_vect(:).*sin_vect(:));
                    end
                end
                Ath(:,ff) = An(1+BORcomp,:).';
                Bth(:,ff) = Bn(1+BORcomp,:).';
                Cth(:,ff) = Cn(1+BORcomp,:).';
                Dth(:,ff) = Dn(1+BORcomp,:).';
                
                BORpower_integrand = 1./(2.*obj1.eta0).*(abs(Ath(:,ff)).^2 + abs(Bth(:,ff)).^2 + abs(Cth(:,ff)).^2 + abs(Dth(:,ff)).^2).*sin(th_vect);
                BORpower(ff) = pi.*integral1D(th_vect,BORpower_integrand);
            end
            % Build a suitable FarField class
            if BORcomp == 0
                [PH,TH] = meshgrid(0,th_vect);
                Eth = Bth;
                Eph = Cth;
                symBOR = 'BOR0';
            else
                % For y-pol: A1 -> Gth and D1 -> Gph
                % For x-pol: B1 -> Gth and C1 -> Gph
                [PH,TH] = meshgrid([0,pi/2],th_vect);
                Eth = [Bth;Ath];  % First element corresponds to ph = 0, and second to ph = pi/2
                Eph = [Cth;Dth];
                symBOR = 'BOR1';
            end
            obj = FarField(PH(:),TH(:),Eth,Eph,obj1.freq,BORpower,obj1.radEff,...
                'coorType','spherical','polType',obj1.polType,'gridType','PhTh','freqUnit',obj1.freqUnit,...
                'symmetryBOR',symBOR,'orientation',obj1.orientation,'earthLocation',obj1.earthLocation,'time',obj1.time,'r',obj1.r);
        end
        
        function obj = expandBORpattern(obj1,phStepDeg)
            % EXPANDBORPATTERN Generates a BOR pattern over a 2*pi phi
            % span.
            
            % Expands a BOR pattern, typically generated by
            % FarField.getBORpattern, into a 2*pi ph span
            
            if nargin < 2, phStepDeg = 5; end
            
            assert(~strcmp(obj1.symmetryBOR,'none'),'Input object not BOR symmetric')
            assert(strcmp(obj1.gridType,'PhTh'),'BOR patterns must be specified on a PhTh grid')
            assert(strcmp(obj1.coorType,'spherical'),'BOR patterns must be specified in a spherical coordinate system')
            assert(obj1.isGridUniform,'A plaid, monotonic, uniform grid is expected for BOR field expansion');
            
            assert(strcmp(obj1.polType,'linear'),'BOR patterns must be specified in a linear polarisation: Circular still TODO')
            
            phStep = deg2rad(phStepDeg);
            Nph = 2*pi/phStep + 1;
            Nth = obj1.Ny;
            ph_vect = linspace(0,2*pi,Nph);
            th_vect = obj1.y(1:Nth);
            [PH,TH] = meshgrid(ph_vect,th_vect);
            if strcmp(obj1.symmetryBOR,'BOR0')
                [A,D] = zeros(size(obj1.E1));
                B = obj1.E1;
                C = obj1.E2;
            else
                [A,B,C,D] = obj1.getBOR1comps;
            end
            [Eth,Eph] = deal(zeros(Nth*Nph,obj1.Nf));
            for ff = 1:obj1.Nf
                Gth = bsxfun(@times,sin(PH),A(:,ff)) + bsxfun(@times,cos(PH),B(:,ff));
                Gph = bsxfun(@times,cos(PH),C(:,ff)) - bsxfun(@times,sin(PH),D(:,ff));
                Eth(:,ff) = Gth(:);
                Eph(:,ff) = Gph(:);
            end
            obj = FarField(PH(:),TH(:),Eth,Eph,obj1.freq,obj1.Prad,obj1.radEff,...
                'coorType','spherical','polType',obj1.polType,'gridType','PhTh','freqUnit',obj1.freqUnit,'slant',obj1.slant,...
                'orientation',obj1.orientation,'earthLocation',obj1.earthLocation,'time',obj1.time,'r',obj1.r);
        end
        
        function [A1,B1,C1,D1] = getBOR1comps(obj1)
            % GETBOR1COMPS Returns components of a BOR1 pattern.
            
            assert(strcmp(obj1.symmetryBOR,'BOR1'),'Input object not BOR1 symmetric')
            assert(strcmp(obj1.gridType,'PhTh'),'BOR1 patterns must be specified on a PhTh grid')
            assert(strcmp(obj1.coorType,'spherical'),'BOR1 patterns must be specified in a spherical coordinate system')
            assert(isequal(unique(obj1.x),[0;pi/2]),'Expect ph cuts only at 0 and pi/2')
            assert(obj1.isGridUniform,'A plaid, monotonic, uniform grid is expected for BOR1 field expansion');
            Nth = obj1.Ny;
            A1 = obj1.E1(Nth+1:end,:);
            B1 = obj1.E1(1:Nth,:);
            C1 = obj1.E2(1:Nth,:);
            D1 = obj1.E2(Nth+1:end,:);
        end
        
        %% Format and other testers
        function y = isGridEqual(obj1,obj2)
            % ISGRIDEQUAL Compares the grid between two FarField objects.
            
            % Dont go for formal equality - floating point error just too much...
            tol = 10^(-obj1.nSigDig+1);
            if all(size(obj1.x) == size(obj2.x)) && all(size(obj1.y) == size(obj2.y))
                xEqual = abs(obj1.x - obj2.x) < tol;
                yEqual = abs(obj1.y - obj2.y) < tol;
                gridEqual = strcmp(obj1.gridType,obj2.gridType);
                fEqual = isequal(obj1.freqHz,obj2.freqHz);
                y = all(xEqual) && all(yEqual) && gridEqual && fEqual;
            else
                y = 0;
            end
        end
        
        function y = typesAreEqual(obj1,obj2)
            % TYPESAREEQUAL Compares if the grid, coorTpe and polarization
            % are equal between two FarField objects.
            
            gridEqual = strcmp(obj1.gridType,obj2.gridType);
            coorEqual = strcmp(obj1.coorType,obj2.coorType);
            polEqual = strcmp(obj1.polType,obj2.polType);
            y = gridEqual && coorEqual && polEqual;
        end
        
        function y = isGrid4pi(obj)
            % ISGRID4PI Check if the data is defined over a full sphere.
            
            % Set to the PhTh coordinate system - this is how most data
            % will be generated anyway.
            % Very quick check - necessary but not always sufficient
            phRange = max(obj.phBase) - min(obj.phBase);
            thRange = max(obj.thBase) - min(obj.thBase);
            eps = 1e-4;
            y = ((abs(round(rad2deg(phRange)) - (360/2^(sum(abs([obj.symXZ,obj.symYZ]))))) < eps) & (abs(round(rad2deg(thRange)) - 180/2^abs(obj.symXY)) < eps)) |...
                ((abs(round(rad2deg(phRange)) - 180) < eps) & (abs(round(rad2deg(thRange)) - 360) < eps));
        end
        
        function y = isGridUniform(obj)
            % ISGRIDUNIFORM Checks for a plaid, monotonic, uniform grid format. 
            
            % Test for a plaid, monotonic, uniform grid format
            
            if obj.Nx*obj.Ny == obj.Nang
                X = reshape(obj.x,obj.Ny,obj.Nx);
                Y = reshape(obj.y,obj.Ny,obj.Nx);
                % Test for equal rows in X, and equal columns in Y (plaid)
                rowEq = all(all(bsxfun(@eq,X,X(1,:))));
                colEq = all(all(bsxfun(@eq,Y,Y(:,1))));
                % Test for monotonic
                diffX = diff(X(1,:));
                diffY = diff(Y(:,1));
                monX = all(diffX>0);
                monY = all(diffY>0);
                % Test for uniform
                tol = 10^(-obj.nSigDig+1);
                unX = all(abs(diffX - median(diffX)) < tol);
                unY = all(abs(diffY - median(diffY)) < tol);
                % And all of them
                y = rowEq && colEq && monX && monY && unX && unY;
            else
                y = false;
            end
        end
        
        %% Astronomical methods
        function obj = setOrientation(obj,newOrientation)
            % SETORIENTATION Set the antenna orientation.
            
            currGridType = obj.gridType;
            obj = obj.grid2PhTh;
            obj.orientation = newOrientation;
            obj = obj.grid2AzAlt;
            eval(['obj = grid2',currGridType,'(obj);']);
        end
        
        function obj = setTime(obj,newTime)
            % SETTIME Set time.
            
            currGridType = obj.gridType;
            obj = obj.grid2AzAlt;
            obj.time = newTime;
            eval(['obj = grid2',currGridType,'(obj);']);
        end
        
        function obj = setEarthLocation(obj,newEarthLocation)
            % SETEARTHLOCATION Set the antenna location on Earth.
            
            currGridType = obj.gridType;
            obj = obj.grid2AzAlt;
            obj.earthLocation = newEarthLocation;
            eval(['obj = grid2',currGridType,'(obj);']);
        end
        
        %% File Output methods
        function writeGRASPcut(obj,pathName)
            % WRITEGRASPCUT Write a FarField object to a GRASP cut file
            % format.
            
            % function [] = writeGRASPcut(obj,pathName)
            % Function to write GRASP cut files to the location PathName
            % (.cut appended automatically)
            % th must be sampled uniformly
            
            % Dirk de Villiers
            % Created: 2019-04-22
            
            assert(obj.isGridUniform,'Standard uniform grid required for GRASP cut write')
            assert(strcmp(obj.gridTypeBase,'PhTh'),'Base grid must be PhTh for GRASP cut write')
            
            % Sort out the data
            th_vect = unique(obj.th);
            ph_vect = unique(obj.ph);
            Nth = obj.Ny;
            Nph = obj.Nx;
            
            thD = mean(diff(th_vect));
            
            % Set up file variables
            HEADER = 'Cut file generated from MATLAB FarField object';
            V_INI = rad2deg(th_vect(1));
            V_INC = rad2deg(thD);
            V_NUM = Nth;
            if strcmp(obj.polType,'linear')
                if strcmp(obj.coorType,'spherical')
                    ICOMP = 1;  % linear th ph polarization...
                    E1real = real(obj.E1);
                    E1imag = imag(obj.E1);
                    E2real = real(obj.E2);
                    E2imag = imag(obj.E2);
                elseif strcmp(obj.coorType,'Ludwig3')
                    ICOMP = 3;  % linear Ludwig3 polarization...
                    E1real = real(obj.E2);
                    E1imag = imag(obj.E2);
                    E2real = real(obj.E1);
                    E2imag = imag(obj.E1);
                end
            elseif strcmp(obj.polType,'circular')
                ICOMP = 2;  % circular polarization...
                E1real = real(obj.E2);
                E1imag = imag(obj.E2);
                E2real = real(obj.E1);
                E2imag = imag(obj.E1);
            end
            ICUT = 1; % Standard polar cut where phi is fixed (nr_cuts) and th is varied
            NCOMP = 2;  % for farfields only 2 components needed
            
            % Create and open the file
            PN = [pathName,'.cut'];
            fid = fopen(PN,'wt');
            
            % Check for identical start and end ph values
            if abs((ph_vect(end)-2*pi) - ph_vect(1)) < eps
                phEnd = Nph-1;
            else
                phEnd = Nph;
            end
            
            for ff = 1:obj.Nf
                for pp = 1:phEnd
                    C = rad2deg(ph_vect(pp));
                    
                    % Header lines
                    form = '%s\n';
                    fprintf(fid,form,[HEADER]);
                    
                    form = '%.10E %.10E %i %.10E %i %i %i\n';
                    fprintf(fid,form,[V_INI, V_INC, V_NUM, C, ICOMP, ICUT, NCOMP]);
                    
                    startPos = (pp-1)*Nth + 1;
                    stopPos = pp*Nth;
                    E1r = E1real(startPos:stopPos,ff);
                    E1i = E1imag(startPos:stopPos,ff);
                    E2r = E2real(startPos:stopPos,ff);
                    E2i = E2imag(startPos:stopPos,ff);
                    
                    form = '%.10E\t%.10E\t%.10E\t%.10E\n';
                    for tt = 1:Nth
                        fprintf(fid,form,[E1r(tt), E1i(tt), E2r(tt), E2i(tt)]);
                    end
                end
            end
            
            fclose(fid);
        end
        
    end
    
    methods (Static = true)
        %% Farfield reading methods
        function FF = readGRASPgrd(pathName,varargin)
            % READGRASPGRD Create a FarFiled object from a GRASP .grd file. 
            
            % function FF = readGRASPgrd(filePathName)
            % Reads a GRASP .grd file and returns the FarField object.
            % Not all GRASP functionality supported yet...
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readGRASPgrd';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR0','BOR1'};
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,pi/2],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(pathName,varargin{:})
            
            pathName = parseobj.Results.pathname;
            symmetryXZ = parseobj.Results.symmetryXZ;
            symmetryYZ = parseobj.Results.symmetryYZ;
            symmetryXY = parseobj.Results.symmetryXY;
            symmetryBOR = parseobj.Results.symmetryBOR;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            if ~strcmp(pathName(end-3:end),'.grd')
                pathName = [pathName,'.grd'];
            end
            fid = fopen(pathName);
            
            E1 = [];
            E2 = [];
            E3 = [];
            
            freqMarker = 'FREQUENCIES';
            startMarker = '++++';
            
            % Read the field header
            while 1
                a = fgetl(fid);
                if strncmp(a,freqMarker,11) % Read the frequencies info
                    freqUnit = regexp(a, '(?<=\[)[^)]*(?=\])', 'match', 'once');
                    % Keep reading lines until all frequencies read
                    freq = [];
                    while 1
                        a = fgetl(fid);
                        if strcmp(a,startMarker), break; end
                        freq = [freq,str2num(a)];
                    end
                end
                if strcmp(a,startMarker)
                    a = fgetl(fid);
                    a = fgetl(fid);
                    fieldInfo = str2num(a);
                    NSET = fieldInfo(1);
                    ICOMP = fieldInfo(2);
                    NCOMP = fieldInfo(3);
                    IGRID = fieldInfo(4);
                    [IX,IY] = deal(zeros(1,NSET));
                    for ff = 1:NSET
                        a = fgetl(fid);
                        centreInfo = str2num(a);
                        IX(ff) = centreInfo(1);
                        IY(ff) = centreInfo(2);
                    end
                    break;
                end
            end
            % Check if all beams have the same grid
            cI1 = centreInfo(1,:);
            comp = any(bsxfun(@minus,centreInfo,cI1),2);
            if sum(comp) > 0
                error('All the field sets must have the same grid - here there are different centre points...');
            end
            
            % Front matter done - read the NSET frequency blocks
            for ff = 1:NSET
                a = fgetl(fid);
                gridInfo1 = str2num(a); % Must be str2num
                XS = gridInfo1(1);
                YS = gridInfo1(2);
                XE = gridInfo1(3);
                YE = gridInfo1(4);
                a = fgetl(fid);
                gridInfo2 = str2num(a); % Must be str2num
                NX = gridInfo2(1);
                NY = gridInfo2(2);
                KLIMIT = gridInfo2(3);
                if KLIMIT == 1
                    %ToDo
                    error('KLIMIT = 1 not implemented yet...')
                else
                    IS = 1;
                    JS = 1;
                    IE = NX;
                    JE = NY;
                end
                if ff == 1  % Just build the grid once - assume they are all the same (ToDo: build a check later)
                    DX = (XE-XS)/(NX - 1);
                    DY = (YE-YS)/(NY - 1);
                    XCEN = DX*IX(ff);
                    YCEN = DY*IY(ff);
                    X = XCEN + XS+DX.*((IS:IE) - 1);
                    Y = YCEN + YS+DY.*((JS:JE) - 1);
                end
                
                if NCOMP == 2
                    form = '%f %f %f %f';
                    fieldData = textscan(fid, form, NX*NY);
                    E1 = [E1,(fieldData{1} + 1i.*fieldData{2})];
                    E2 = [E2,(fieldData{3} + 1i.*fieldData{4})];
                elseif NCOMP == 3
                    form = '%f %f %f %f %f %f';
                    fieldData = textscan(fid, form, NX*NY);
                    E1 = [E1,(fieldData{1} + 1i.*fieldData{2})];
                    E2 = [E2,(fieldData{3} + 1i.*fieldData{4})];
                    E3 = [E3,(fieldData{5} + 1i.*fieldData{6})];
                end
                % Dummy read
                a = fgetl(fid);
            end
            fclose(fid);
            
            % Build the object
            [Xmat,Ymat] = ndgrid(X,Y);
            x = Xmat(:);
            y = Ymat(:);
            switch ICOMP
                case 1
                    polType = 'linear';
                    coorType = 'spherical';
                case 2
                    polType = 'circular';
                    coorType = 'spherical';
                case 3
                    polType = 'linear';
                    coorType = 'Ludwig3';
                otherwise
                    error(['ICOMP ',num2str(ICOMP),' case not implemented yet'])
            end
            switch IGRID
                case 1
                    gridType = 'DirCos';
                case 4
                    gridType = 'AzEl';
                case 5
                    gridType = 'TrueView';
                case 6
                    gridType = 'ElAz';
                case 7
                    gridType = 'PhTh';
                otherwise
                    error(['IGRID ',num2str(IGRID),' case not implemented yet'])
            end
           
            % keyboard;
            Prad = ones(size(freq)).*4*pi;
            radEff = ones(size(freq));
            FF = FarField(x,y,E1,E2,freq,Prad,radEff,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'symmetryXZ',symmetryXZ,'symmetryYZ',symmetryYZ,'symmetryXY',symmetryXY,'symmetryBOR',symmetryBOR,'E3',E3,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
        end

        function FF = readFEKOffe(pathName,varargin)
            % READFEKOFFE Create a FarFiled object from a FEKO .ffe file.
            
            %Name: readFEKOffe.m
            %Description:
            %   Function to create a Farfield object from a FEKO .ffe far-field output
            %   file.
            %Inputs:
            % th: column vector [Nang x 1] of th angles in rad
            % ph: column vector [Nang x 1] of ph angles in rad
            %Outputs:
            % --FF: Farfield object containing parameters as derived from target .ffe
            % file.
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readFEKOffe';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR0','BOR1'};
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,pi/2],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(pathName,varargin{:})
            
            pathName = parseobj.Results.pathname;
            symmetryXZ = parseobj.Results.symmetryXZ;
            symmetryYZ = parseobj.Results.symmetryYZ;
            symmetryXY = parseobj.Results.symmetryXY;
            symmetryBOR = parseobj.Results.symmetryBOR;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            eta0 = 3.767303134749689e+02;
            
            %Open the data file
            if ~strcmp(pathName(end-3:end),'.ffe')
                pathName = [pathName,'.ffe'];
            end
            fid = fopen(pathName);
            if (fid==-1)
                error(['Unable to open data file ' fileName '!']);
            end
            
            % Read the main header info
            coorSysMarker = '#Coordinate System:';
            freqMarker = '#Frequency:';
            NthMarker = '#No. of Theta Samples:';
            NphMarker = '#No. of Phi Samples:';
            fieldMarker = '#"Theta""Phi"';
            
            %===================================================================
            % LOAD DATA
            %===================================================================
            
            fCount = 0;
            read = 1;
            while read
                a = fgetl(fid);
                if a == -1
                    read = 0;
                    break;
                end
                if strncmp(a,coorSysMarker,length(coorSysMarker)) % Read the coordinate system type
                    coorSysCell = textscan(a,'%s%s%s');
                    coorType = coorSysCell{3};
                end
                if strncmp(a,freqMarker,length(freqMarker)) % Read the number of frequencies
                    freqCell = textscan(a,'%s%n');
                    freq(fCount+1) = freqCell{2};
                end
                if strncmp(a,NthMarker,length(NthMarker))
                    NthCell = textscan(a,'%s%s%s%s%n');
                    Nth = NthCell{5};
                end
                if strncmp(a,NphMarker,length(NphMarker))
                    NphCell = textscan(a,'%s%s%s%s%n');
                    Nph = NphCell{5};
                end
                
                %     if strncmp(a,fieldMarker,length(fieldMarker))
                aNoSpace = a;
                aNoSpace(ismember(a,' ')) = [];
                if strncmp(aNoSpace,fieldMarker,13)
                    %         keyboard;
                    %         fieldHeader = strsplit(a);
                    %         Ncomp = length(fieldHeader)-1;  % Count how many columns to expect from the header
                    Ncomp = sum(ismember(a,'"'))/2; % Count how many columns to expect from the header
                    getP = false;
                    if Ncomp >=9, getP = true; end  % Assume column 9 is directivity..
                    if fCount == 0
                        % Initialise matrices
                        NfInit = 1; % Use as initial maximum guess...
                        [th,ph,Eth,Eph,Dth,Dph,Dtotclear] = deal(zeros(Nth*Nph,NfInit));
                        fData = zeros(Nth*Nph,Ncomp,NfInit);
                        Prad = zeros(1,NfInit);
                    end
                    fCount = fCount + 1;
                    fDataForm = repmat('%f',1,Ncomp);
                    fData(:,:,fCount) = fscanf(fid,fDataForm,[Ncomp,Nth*Nph])';
                    th(:,fCount) = deg2rad(fData(:,1,fCount));
                    ph(:,fCount) = deg2rad(fData(:,2,fCount));
                    Eth(:,fCount) = fData(:,3,fCount) + 1i.*fData(:,4,fCount);
                    Eph(:,fCount) = fData(:,5,fCount) + 1i.*fData(:,6,fCount);
                    %         Dth(:,fCount) = fData(:,7,fCount);
                    %         Dph(:,fCount) = fData(:,8,fCount);
                    %         Dtot(:,fCount) = fData(:,9,fCount);
                    % Get the power from the directivity/gain
                    if getP
                        U = 1/(2*eta0).*(abs(Eth(:,fCount)).^2 + abs(Eph(:,fCount)).^2);
                        Prad(fCount) = median(4.*pi.*U./10.^(fData(:,9,fCount)./10));
                    else 
                        Prad = [];
                    end
                end
            end
            % Build the object
            x = ph(:,1);
            y = th(:,1);
            
            E1 = Eth;
            E2 = Eph;
           
            coorType = lower(coorType); %coordinate system string fetched from header text of .ffe file
            polType = 'linear'; %It seems that FEKO always outputs linear polarised fields (corresponding to th-ph coordinates)
            gridType = 'PhTh'; %It seems that FEKO always outputs .ffe field values in theta-phi form, so this is hardcoded as such
            %NB: Prad defined earlier
            radEff = ones(size(freq)); %replace with manual radiation efficiency calculation
            freqUnit = 'Hz'; %It seems that FEKO always outputs frequencies in Hz, so this is hardcoded as such
            
            FF = FarField(x,y,E1,E2,freq,Prad,radEff,...
                'coorType',coorType{1},'polType',polType,'gridType',gridType,'freqUnit',freqUnit,'r',1,...
                'symmetryXZ',symmetryXZ,'symmetryYZ',symmetryYZ,'symmetryXY',symmetryXY,'symmetryBOR',symmetryBOR,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
        end

        function FF = readCSTffs(pathName,varargin)
            % READCSTFFS Create a FarFiled object from a CST .ffs file.
            
            % function [FF] = readCSTffs(pathName)
            % Loads a CST generated farfield source file pathName.ffs.
            %
            %
            % Inputs:
            % path_name - Full path and filename (no extension) string
            % Outputs:
            % FF - standard farfield object
            %
            % Dirk de Villiers
            % Created: 2015-02-03
            % Last edit: 2019-18-02
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readCSTffs';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,pi/2],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(pathName,varargin{:})
            
            pathName = parseobj.Results.pathname;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            % Open the data file
            if ~strcmp(pathName(end-3:end),'.ffs')
                pathName = [pathName,'.ffs'];
            end
            fid = fopen(pathName);
            if (fid==-1)
                error(['Unable to open data file ', fileName, '!']);
            end
            
            %===================================================================
            % LOAD DATA
            %===================================================================
            
            % Read the main header info
            freqMarker = '// #Frequencies';
            posMarker = '// Position';
            zAxisMarker = '// zAxis';
            xAxisMarker = '// xAxis';
            powerFreqMarker = '// Radiated/Accepted/Stimulated Power , Frequency';
            NphNthMarker = '// >> Total #phi samples, total #theta samples';
            fieldMarker = '// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi):';
            
            fCount = 0;
            read = 1;
            while read
                a = fgetl(fid);
                if strcmp(a,freqMarker) % Read the number of frequencies
                    Nf = fscanf(fid,'%i',1);
                end
                if strcmp(a,posMarker) % Read the position
                    pos = fscanf(fid,'%f%f%f',[3,1]);
                end
                if strcmp(a,zAxisMarker) % Read the zAxis
                    zAxis = fscanf(fid,'%f%f%f',[3,1]);
                end
                if strcmp(a,xAxisMarker) % Read the xAxis
                    xAxis = fscanf(fid,'%f%f%f',[3,1]);
                end
                if strncmp(a,powerFreqMarker,length(powerFreqMarker))
                    PF = fscanf(fid,'%f',[Nf*4,1]);
                    [Prad,Pacc,Pstim,freq] = deal(zeros(1,Nf));
                    for ii = 1:Nf
                        Prad(ii) = PF((ii-1)*4+1);
                        Pacc(ii) = PF((ii-1)*4+2);
                        Pstim(ii) = PF((ii-1)*4+3);
                        freq(ii) = PF((ii-1)*4+4);
                    end
                end
                if strncmp(a,NphNthMarker,length(NphNthMarker))
                    NphNth = fscanf(fid,'%i %i',[2,1]);
                    Nph = NphNth(1);
                    Nth = NphNth(2);
                    if fCount == 0
                        % Initialise matrices
                        [th,ph,Eth,Eph] = deal(zeros(Nth*Nph,Nf));
                    end
                end
                
                if strncmp(a,fieldMarker,length(fieldMarker))
                    fCount = fCount + 1;
                    fDataForm = '%f %f %f %f %f %f';
                    fData = fscanf(fid,fDataForm,[6,Nth*Nph])';
                    th(:,fCount) = deg2rad(fData(:,2));
                    ph(:,fCount) = deg2rad(fData(:,1));
                    Eth(:,fCount) = fData(:,3) + 1i.*fData(:,4);
                    Eph(:,fCount) = fData(:,5) + 1i.*fData(:,6);
                    if fCount == Nf
                        read = 0;
                    end
                end
            end
            
            fclose(fid);
            
            % Build the FF obj
            x = ph(:,1);
            y = th(:,1);
            E1 = Eth;
            E2 = Eph;
            radEff = Prad./Pacc;
            coorType = 'spherical';
            polType = 'linear';
            gridType = 'PhTh';
            freqUnit = 'Hz';
            
            FF = FarField(x,y,E1,E2,freq,Prad,radEff,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
            
        end

        function FF = readGRASPcut(pathName,nr_freq,nr_cuts,varargin)
            % READGRASPCUT Create a FarFiled object from a GRASP .cut file.
            
            % [FF] = readGRASPcut(pathName)
            % Loads a GRASP generated farfield source file pathName.cut.
            %
            %
            % Inputs:
            % path_name - Full path and filename string
            % nr_freq - number of frequency points
            % nr_cuts - number of cuts taken
            % Outputs:
            % FF - standard farfield object
            %
            % Dirk de Villiers
            % Created: 2019-04-22
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readGRASPcut';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addRequired('nr_freq',typeValidation_scalar);
            parseobj.addRequired('nr_cuts',typeValidation_scalar);
            
            typeValidation_freq = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','increasing','nrows',1},'readGRASPcut');
            parseobj.addOptional('freq',1,typeValidation_freq);
            
            expected_freqUnit = {'Hz','kHz','MHz','GHz','THz'};
            parseobj.addParameter('freqUnit','Hz', @(x) any(validatestring(x,expected_freqUnit)));
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR0','BOR1'};
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,pi/2],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(pathName,nr_freq,nr_cuts,varargin{:})
            
            pathName = parseobj.Results.pathname;
            nr_freq = parseobj.Results.nr_freq;
            nr_cuts = parseobj.Results.nr_cuts;
            freq = parseobj.Results.freq;
            freqUnit = parseobj.Results.freqUnit;
            symmetryXZ = parseobj.Results.symmetryXZ;
            symmetryYZ = parseobj.Results.symmetryYZ;
            symmetryXY = parseobj.Results.symmetryXY;
            symmetryBOR = parseobj.Results.symmetryBOR;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            %Open the data file
            if ~strcmp(pathName(end-3:end),'.cut')
                pathName = [pathName,'.cut'];
            end
%             global fid;
            fid = fopen(pathName,'rt');
%             global fid;
%             [fid, message] = fopen([pathName,'.cut'], 'rt');
            if (fid==-1)
                error(['Unable to open data file ' fileName '!']);
            end
            
            eta0 = 3.767303134749689e+02;
            
            %===================================================================
            % Load data for pre-allocation
            %===================================================================
            % Skip over text line
            dummy = fgetl(fid);
%             form= '%*s %*s %*s %*s';
%             dummy = textscan(fid, form, 1);
            
            % Read info line
            form= '%f %f %f %f %f %f %f';
            cut_info = textscan(fid, form, 1);
            V_INI = cut_info{1};
            V_INC = cut_info{2};
            V_NUM = cut_info{3};
            C = cut_info{4};
            ICOMP = cut_info{5};
            ICUT = cut_info{6};
            % NCOMP = cut_info{7};
            
            % Preallocate
            [th_deg,ph_deg] = deal(zeros(V_NUM*nr_cuts,1));
            [E1,E2] = deal(zeros(V_NUM*nr_cuts,nr_freq));
            
            for ff = 1:nr_freq
                for cc = 1:nr_cuts
                    if ff == 1 % Only do once
                        x_1cut = ones(V_NUM,1).*C;
                        y_1cut = (V_INI:V_INC:(V_INC*(V_NUM - 1) + V_INI)).';
                        
                        if ICUT == 1
                            ph_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = x_1cut;
                            th_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = y_1cut;
                        elseif ICUT == 2
                            th_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = x_1cut;
                            ph_deg(((cc-1)*V_NUM + 1):cc*V_NUM) = y_1cut;
                        end
                    end
                    % Read cut data
                    form= '%f %f %f %f';
                    cut_data = textscan(fid, form, V_NUM);
                    E1(((cc-1)*V_NUM + 1):cc*V_NUM,ff) = cut_data{1} + 1i.*cut_data{2};
                    E2(((cc-1)*V_NUM + 1):cc*V_NUM,ff) = cut_data{3} + 1i.*cut_data{4};
                    
                    dummy = fgetl(fid); % step off previous line
                    dummy = fgetl(fid); % step over header line
                    form= '%f %f %f %f %f %f %f';
                    cut_info = textscan(fid, form, 1);
                    C = cut_info{4};
                end
            end
            fclose(fid);
            
            switch abs(ICOMP)
                case 1
                    polType = 'linear';
                    coorType = 'spherical';
                    E1ff = E1;
                    E2ff = E2;
                case 2
                    polType = 'circular';
                    coorType = 'spherical';
                    E1ff = E2;
                    E2ff = E1;
                case 3
                    polType = 'linear';
                    coorType = 'Ludwig3';
                    E1ff = E2;
                    E2ff = E1;
                otherwise
                    error(['ICOMP ',num2str(ICOMP),' case not implemented yet'])
            end
            gridType = 'PhTh';
            
            % Build the FF obj
            x = deg2rad(ph_deg);
            y = deg2rad(th_deg);
            Prad = ones(1,nr_freq).*4*pi/(2*eta0);
            radEff = ones(1,nr_freq);
            if isscalar(freq)
                freq = repmat(freq,1,nr_freq);
            end
            
            FF = FarField(x,y,E1ff,E2ff,freq,Prad,radEff,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'symmetryXZ',symmetryXZ,'symmetryYZ',symmetryYZ,'symmetryXY',symmetryXY,'symmetryBOR',symmetryBOR,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
        end

        function FF = farFieldFromPowerPattern(x,y,P,freq,varargin)
            % FARFIELDFROMPOWERPATTERN Create a Farfield object from a
            % power pattern.
            
            %Name: farFieldFromPowerPattern.m
            %Description:
            %   Function to create a Farfield object from a power pattern, which must
            %   be provided in the format that is outputted by powerPattern.m.
            %Inputs:
            % --ph: column vector [Nang x 1] of ph angles in rad
            % --th: column vector [Nang x 1] of th angles in rad
            % --P: column vector [Nang x 1] of power pattern values (see Pout output for powerPattern.m for an example of the format)
            % --freq: scalar frequency in Hz
            % --fieldPol: (optional) string denoting how far field should be polarized- options are 'linearX', 'linearY', 'circularLH', 'circularRH'
            % --freqUnit: (optional) string denoting the frequency unit {('Hz'),'kHz','MHz','GHz','THz'}
            %Outputs:
            % --FF: Farfield object containing parameters as determined from inputs
            
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'farFieldFromPowerPattern';
            
            typeValidation_grid = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','ncols',1},'farFieldFromPowerPattern');
            parseobj.addRequired('x',typeValidation_grid);
            parseobj.addRequired('y',typeValidation_grid);
            
            typeValidation_P = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan'},'farFieldFromPowerPattern');
            parseobj.addRequired('P',typeValidation_P);
            
            typeValidation_freq = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','increasing','nrows',1},'farFieldFromPowerPattern');
            parseobj.addRequired('freq',typeValidation_freq);
            
            expected_freqUnit = {'Hz','kHz','MHz','GHz','THz'};
            parseobj.addParameter('freqUnit','Hz', @(x) any(validatestring(x,expected_freqUnit)));
            
            expected_fieldPol = {'linearX','linearY','circularLH','circularRH'};
            parseobj.addParameter('fieldPol','linearY', @(x) any(validatestring(x,expected_fieldPol)));
            
            expected_gridType = {'PhTh','DirCos','AzEl','ElAz','AzAlt','TrueView','ArcSin','Mollweide','RAdec','GalLongLat'};
            parseobj.addParameter('gridType','PhTh', @(x) any(validatestring(x,expected_gridType)));
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR1'};  % Can't do BOR0 symmetry
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,pi/2],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(x,y,P,freq,varargin{:})
            
            x = parseobj.Results.x;
            y = parseobj.Results.y;
            P = parseobj.Results.P;
            freq = parseobj.Results.freq;
            fieldPol = parseobj.Results.fieldPol;
            gridType = parseobj.Results.gridType;
            freqUnit = parseobj.Results.freqUnit;
            symmetryXZ = parseobj.Results.symmetryXZ;
            symmetryYZ = parseobj.Results.symmetryYZ;
            symmetryXY = parseobj.Results.symmetryXY;
            symmetryBOR = parseobj.Results.symmetryBOR;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            %constants
            eta0 = 3.767303134749689e+02;
            tol = 1e-15 ;
            
            if isscalar(P)
                P = P./(4.*pi).*ones(size(ph));
            end
            %From power pattern and polarization parameters, generate E1 and E2 accordingly
            if strcmp(symmetryBOR,'none')
                coorType = 'Ludwig3';
                switch fieldPol
                    case 'linearX' % linearly polarised along X-axis
                        polType = 'linear';
                        E1  = sqrt(P.*2*eta0);
                        E2  = zeros(size(P));
                    case 'linearY' % linearly polarised along Y-axis
                        polType = 'linear';
                        E1  = zeros(size(P));
                        E2  = sqrt(P.*2*eta0);
                    case 'circularLH'  % Lefthand Circular polarization
                        polType = 'circular';
                        E1  = sqrt(P.*2*eta0);
                        E2  = zeros(size(P));
                    case 'circularRH'  % Righthand Circular polarization
                        polType = 'circular';
                        E1  = zeros(size(P));
                        E2  = sqrt(P.*2*eta0);
                    otherwise
                        error('fieldPol input string unrecognised')
                end
            else
                assert(strcmp(gridType,'PhTh'),'PhTh grid required for BOR1 definition')
                assert(all(abs(unique(x) - [0;pi/2]) < tol),'Invalid range for BOR1 symmetry (E-plane and H-plane required)')
                coorType = 'spherical';
                iph0 = find(abs(x - 0) < tol);
                iph90 = find(abs(x - pi/2) < tol);
                [E1,E2] = deal(zeros(size(P)));
                switch fieldPol
                    case 'linearX' 
                        polType = 'linear';
                        E1(iph0,:)  = sqrt(P(iph0,:).*2*eta0);
                        E2(iph90,:)  = -sqrt(P(iph90,:).*2*eta0);
                    case 'linearY'
                        polType = 'linear';
                        E1(iph90,:)  = sqrt(P(iph90,:).*2*eta0);
                        E2(iph0,:)  = sqrt(P(iph0,:).*2*eta0);
                    otherwise
                        error('BOR1 functionality for CP not yet implemented')
                end
            end
            
            % Build the object - use default Prad and radEff
            FF = FarField(x,y,E1,E2,freq,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'symmetryXZ',symmetryXZ,'symmetryYZ',symmetryYZ,'symmetryXY',symmetryXY,'symmetryBOR',symmetryBOR,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
        end
        
    end
    
    %% Internal helper functions
    methods (Access = private)
        %% Grid getters
        function [u, v, w] = getDirCos(obj)
            % GETDIRCOS Get DirCos grid.
            
            handle2DirCos = str2func([obj.gridTypeBase,'2DirCos']);
            [u,v,w] = handle2DirCos(obj.xBase,obj.yBase);
        end
        
        function [ph, th] = getPhTh(obj)
            % GETPHTH Get PhTh grid.
            
            switch obj.gridTypeBase
                case 'PhTh'
                    ph = obj.xBase;
                    th = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [ph,th] = DirCos2PhTh(u,v,w);
            end
        end
        
        function [az, el] = getAzEl(obj)
            % GETAZEL Get AzEl grid.
            
            switch obj.gridTypeBase
                case 'AzEl'
                    el = obj.yBase;
                    az = obj.xBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [az,el] = DirCos2AzEl(u,v,w);
            end
        end
        
        function [ep, al] = getElAz(obj)
            % GETELAZ Get ElAz grid.
            
            switch obj.gridTypeBase
                case 'ElAz'
                    ep = obj.xBase;
                    al = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [ep,al] = DirCos2ElAz(u,v,w);
            end
        end
        
        function [Xg, Yg] = getTrueView(obj)
            % GETTRUEVIEW Get TrueView grid.
            
            switch obj.gridTypeBase
                case 'TrueView'
                    Xg = obj.xBase;
                    Yg = obj.xBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [Xg,Yg] = DirCos2TrueView(u,v,w);
            end
        end
        
        function [asinu, asinv] = getArcSin(obj)
            % GETARCSIN Get ArcSin grid.
            
            switch obj.gridTypeBase
                case 'ArcSin'
                    asinu = obj.xBase;
                    asinv = obj.yBase;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [asinu,asinv] = DirCos2ArcSin(u,v,w);
            end
        end
        
        function [az,alt] = getAzAlt(obj)
            % GETAZALT Get AzAlt grid.
            
            switch obj.gridTypeBase
                case 'AzAlt'
                    az = obj.xBase;
                    alt = obj.yBase;
                case 'RAdec'
                    horzCoords= wrap2pi(celestial.coo.horiz_coo([obj.x obj.y],juliandate(obj.time),obj.earthLocation(1:2),'h'));
                    az = horzCoords(:,1); %right ascension
                    alt = horzCoords(:,2); %declination
                case 'GalLongLat'
                    obj1 = obj.grid2RAdec;
                    horzCoords= wrap2pi(celestial.coo.horiz_coo([obj1.x obj1.y],juliandate(obj1.time),obj1.earthLocation(1:2),'h'));
                    az = horzCoords(:,1); %right ascension
                    alt = horzCoords(:,2); %declination
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [az,alt] = DirCos2AzAlt(u,v,w);
            end
        end
        
        function [RA, dec] = getRAdec(obj)
            % GETRADEC Get RAdec grid.
            
            switch obj.gridTypeBase
                case 'AzAlt'
                    equCoords= wrap2pi(celestial.coo.horiz_coo([obj.x obj.y],juliandate(obj.time),obj.earthLocation(1:2),'e'));
                    RA = equCoords(:,1); %right ascension
                    dec = equCoords(:,2); %declination
                case 'RAdec'
                    RA = obj.xBase;
                    dec = obj.yBase;
                case 'GalLongLat'
                    [equCoords,~] = celestial.coo.coco([obj.x obj.y],'g','j2000.0','r','r');
                    RA = wrap2pi(equCoords(:,1)); %right ascension
                    dec = wrap2pi(equCoords(:,2)); %declination
                case obj.projectionGrids
                    [u,v,w] = getDirCos(obj);
                    [RA,dec] = DirCos2AzAlt(u,v,w);
                otherwise
                    error('Grid type must either be AzAlt, GalLongLat or a projection grid')
            end
        end
        
        function [long, lat] = getGalLongLat(obj)
            % GETGALLINGLAT Get GalLongLat grid.
            
            switch obj.gridTypeBase
                case 'PhTh'
                    obj1 = obj.grid2AzAlt;
                    obj1 = obj1.grid2RAdec;
                    [galCoords,~] = celestial.coo.coco([obj1.x obj1.y],'j2000.0','g','r','r');
                    long = wrap2pi(galCoords(:,1)); %galactic longitude
                    lat = galCoords(:,2); %galactic latitude
                case 'AzAlt'
                    obj1 = obj.grid2RAdec;
                    [galCoords,~] = celestial.coo.coco([obj1.x obj1.y],'j2000.0','g','r','r');
                    long = wrap2pi(galCoords(:,1)); %galactic longitude
                    lat = galCoords(:,2); %galactic latitude
                case 'RAdec'
                    [galCoords,~] = celestial.coo.coco([obj.x obj.y],'j2000.0','g','r','r');
                    long = wrap2pi(galCoords(:,1)); %galactic longitude
                    lat = galCoords(:,2); %galactic latitude
                case 'GalLongLat'
                    long = obj.xBase;
                    lat = obj.yBase;
                case obj.projectionGrids
                    [u,v,w] = getDirCos(obj);
                    [long,lat] = DirCos2AzAlt(u,v,w);
                otherwise
                    error('Grid type must either be AzAlt, GalLongLat or a projection grid')
                    
            end
        end
        
        %% Coordinate system getters
        function [Eth, Eph, Er] = getEspherical(obj)
            % GETESPHERICAL Get Espherical coordinates.
            
            [Ph,Th] = getPhTh(obj);
            TH = repmat(Th(:,1),1,obj.Nf);
            PH = repmat(Ph(:,1),1,obj.Nf);
            % Change to the Base values here....
            switch obj.coorTypeBase
                case 'spherical'
                    Eth = obj.E1Base;
                    Eph = obj.E2Base;
                case 'Ludwig1'
                    Eth = cos(TH).*cos(PH).*obj.E1Base + cos(TH).*sin(PH).*obj.E2Base - sin(TH).*obj.E3Base;
                    Eph = -sin(PH).*obj.E1Base + cos(PH).*obj.E2Base;
                case 'Ludwig2AE'
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Del = cos(PH).^2 + cos(TH).^2.*sin(PH).^2;
                    Eth = (cosEl./Del).*(cos(PH).*obj.E1Base + cos(TH).*sin(PH).*obj.E2Base);
                    Eph = (cosEl./Del).*(-cos(TH).*sin(PH).*obj.E1Base + cos(PH).*obj.E2Base);
                case 'Ludwig2EA'
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Del = cos(TH).^2.*cos(PH).^2 + sin(PH).^2;
                    Eth = (cosAl./Del).*(cos(TH).*cos(PH).*obj.E1Base + sin(PH).*obj.E2Base);
                    Eph = (cosAl./Del).*(-sin(PH).*obj.E1Base + cos(TH).*cos(PH).*obj.E2Base);
                case 'Ludwig3'
                    Del = 1;
                    Eth = (1./Del).*(cos(PH).*obj.E1Base + sin(PH).*obj.E2Base);
                    Eph = (1./Del).*(-sin(PH).*obj.E1Base + cos(PH).*obj.E2Base);
            end
            Er = zeros(size(Eth));
        end
        
        function [Ex, Ey, Ez] = getELudwig1(obj)
            % GETELUDWIG1 Get ELudwig1 coordinates.
            
            switch obj.coorTypeBase
                case 'Ludwig1'
                    Ex = obj.E1Base;
                    Ey = obj.E2Base;
                    Ez = obj.E3Base;
                otherwise
                    [Eth, Eph, ~] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    % Assume farfield so no radial E-field
                    Ex = cos(TH).*cos(PH).*Eth - sin(PH).*Eph;
                    Ey = cos(TH).*sin(PH).*Eth + cos(PH).*Eph;
                    Ez = zeros(size(Ex));   % Strict definition in the Ludwig paper
            end
        end
        
        function [Eaz, Eel, E3] = getELudwig2AE(obj)
            % GETELUDWIG2AE Get ELudwig2AE coordinates.
            
            switch obj.coorTypeBase
                case 'Ludwig2AE'
                    Eaz = obj.E1Base;
                    Eel = obj.E2Base;
                    E3 = obj.E3Base;
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Eaz = (1./cosEl).*(cos(PH).*Eth - cos(TH).*sin(PH).*Eph);
                    Eel = (1./cosEl).*(cos(TH).*sin(PH).*Eth + cos(PH).*Eph);
                    E3 = zeros(size(Eaz));
                    % Sort out singularities poles - ToDo
                    phPoles = deg2rad([-270,-90,90,270].');
                    poleMat = [ones(4,1).*deg2rad(90),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[obj.th,obj.ph],'rows');
                    iPole = iPole(iPole>0);
                    [Eaz(iPole,:),Eel(iPole,:),E3(iPole,:)] = deal(0);
                    %                     Eaz(iPole,:) = Eth(iPole,:) + Eph(iPole,:);
            end
        end
        
        function [Eal, Eep, E3] = getELudwig2EA(obj)
            % GETELUDWIG2EA Get ELudwig2EA coordinates.
            
            switch obj.coorTypeBase
                case 'Ludwig2EA'
                    Eal = obj.E1Base;
                    Eep = obj.E2Base;
                    E3 = obj.E3Base;
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Eal = (1./cosAl).*(cos(TH).*cos(PH).*Eth - sin(PH).*Eph);
                    Eep = (1./cosAl).*(sin(PH).*Eth + cos(TH).*cos(PH).*Eph);
                    E3 = zeros(size(Eal));
                    % Sort out singularities poles - ToDo
                    phPoles = deg2rad([-360,-180,0,180,360].');
                    poleMat = [ones(5,1).*deg2rad(90),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[obj.th,obj.ph],'rows');
                    iPole = iPole(iPole>0);
                    [Eal(iPole,:),Eep(iPole,:),E3(iPole,:)] = deal(0);
            end
        end
        
        function [Eh, Ev, E3] = getELudwig3(obj)
            % GETELUDWIG3 Get ELudwig3 coordinates.
            
            switch obj.coorTypeBase
                case 'Ludwig3'
                    Eh = obj.E1Base;
                    Ev = obj.E2Base;
                    E3 = obj.E3Base;
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    Eh = cos(PH).*Eth - sin(PH).*Eph;
                    Ev = sin(PH).*Eth + cos(PH).*Eph;
                    E3 = zeros(size(Eh));
            end
        end
        
        function obj = setBase(obj)
            % SETBASE Set base grid.
            
            obj.xBase = obj.x;
            obj.yBase = obj.y;
            obj.NxBase = obj.Nx;
            obj.NyBase = obj.Ny;
            obj.phBase = obj.ph;
            obj.thBase = obj.th;
            obj.E1Base = obj.E1;
            obj.E2Base = obj.E2;
            obj.E3Base = obj.E3;
            obj.gridTypeBase = obj.gridType;
            obj.coorTypeBase = obj.coorType;
            obj.polTypeBase = obj.polType;
        end
        
        function [ph,th] = getPhThCurrent(obj)
            % GETPHTHCURRENT Get current PhTh grid.
            
            % This does not operate on the base grid, but instead the
            % current grid
            if strcmp(obj.gridType,'PhTh')
                ph = obj.x;
                th = obj.y;
            else
                handle2DirCos = str2func([obj.gridType,'2DirCos']);
                [u,v,w] = handle2DirCos(obj.x,obj.y);
                [ph,th] = DirCos2PhTh(u,v,w);
            end
        end
        
        %% Polarization type getters
        function [E1lin, E2lin, E3lin] = getElin(obj)
            % GETELIN Get linear polarization.
            
            % Start at the base, transform to correct coordinate system
            coorTypeIn = obj.coorType;
            coorTypeH = str2func(['coor2',coorTypeIn]);
            obj1 = obj.reset2Base;
            obj1 = coorTypeH(obj1,false);
            switch obj.polTypeBase % Should be the same as the transformed object - can use obj or obj1
                case 'linear'
                    E1lin = obj1.E1;
                    E2lin = obj1.E2;
                case 'circular'
                    Del = 2*1i;
                    E1lin = sqrt(2)./Del.*(1i.*obj1.E1 + 1i.*obj1.E2);
                    E2lin = sqrt(2)./Del.*(-obj1.E1 + obj1.E2);
                case 'slant'
                    PSI = ones(size(obj1.E1)).*obj.slant; % Slant of input object
                    Del = 1;
                    E1lin = 1./Del.*(cos(PSI).*obj1.E1 + sin(PSI).*obj1.E2);
                    E2lin = 1./Del.*(-sin(PSI).*obj1.E1 + cos(PSI).*obj1.E2);
            end
            E3lin = zeros(size(E1lin));
        end
        
        function [Elh,Erh,E3circ] = getEcircular(obj)
            % GETECIRCULAR Get circular polarization.
            
            switch obj.polType
                case 'circular'
                    Elh = obj.E1;
                    Erh = obj.E2;
                    E3circ = obj.E3;
                otherwise
                    [E1lin, E2lin] = getElin(obj);
                    Elh = 1/sqrt(2).*(E1lin - 1i.*E2lin);
                    Erh = 1/sqrt(2).*(E1lin + 1i.*E2lin);
                    E3circ = zeros(size(Elh));
            end
        end
        
        function [Exp,Eco,E3slant] = getEslant(obj)
            % GETESLANT Get slant polarization.
            
            switch obj.polType
                case 'slant'
                    Exp = obj.E1;
                    Eco = obj.E2;
                    E3slant = obj.E3;
                otherwise
                    [E1lin, E2lin] = getElin(obj);
                    PSI = ones(size(obj.E1)).*obj.slant;
                    Exp = cos(PSI).*E1lin - sin(PSI).*E2lin;
                    Eco = sin(PSI).*E1lin + cos(PSI).*E2lin;
                    E3slant = zeros(size(Exp));
            end
        end
        
        % Set the names of the 2 grid components
        function [xname,yname] = setXYnames(obj)
            % SETXYNAMES Set x and y names
            
            switch obj.gridType
                case 'PhTh'
                    xname = '\phi';
                    yname = '\theta';
                case 'DirCos'
                    xname = 'u';
                    yname = 'v';
                case 'AzEl'
                    xname = 'az';
                    yname = 'el';
                case 'ElAz'
                    xname = '\epsilon';
                    yname = '\alpha';
                case 'TrueView'
                    xname = 'Xg=\theta cos(\phi)';
                    yname = 'Yg=\theta sin(\phi)';
                case 'ArcSin'
                    xname = 'Xg=asin(u)';
                    yname = 'Yg=asin(v)';
                case 'Mollweide'
                    xname = 'Xg';
                    yname = 'Yg';
                case 'AzAlt'
                    xname = 'North-az';
                    yname = 'alt';
                case 'RAdec'
                    xname = 'RA';
                    yname = 'dec';
                case 'GalLongLat'
                    xname = 'long';
                    yname = 'lat';
            end
        end
        
        % Set the names of the 2 farfield components based on the
        % polarization type.  Names used for info and plotting.
        function [E1name,E2name] = setEnames(obj)
            % SETENAMES Set E-field component names
            
            switch obj.polType
                case 'circular'
                    E1name = 'Elh';
                    E2name = 'Erh';
                case 'slant'
                    E1name = 'Exp';
                    E2name = 'Eco';
                case 'linear'
                    switch obj.coorType
                        case 'spherical'
                            E1name = 'Eth';
                            E2name = 'Eph';
                        case 'Ludwig1'
                            E1name = 'Ex';
                            E2name = 'Ey';
                        case 'Ludwig2AE'
                            E1name = 'Eaz';
                            E2name = 'Eel';
                        case 'Ludwig2EA'
                            E1name = 'Eal';
                            E2name = 'Eep';
                        case 'Ludwig3'
                            E1name = 'Eh';
                            E2name = 'Ev';
                    end
                    
                otherwise
                    error(['Unknown polType property: ', obj.polType]);
            end
        end
        
        function [xRangeType,yRangeType] = setRangeTypes(obj)
            % SETRANGETYPES Returns current rangeType
            
            % Try to figure out what the current rangeType is.
            % Not much error checking is done - assume somewhat
            % sensible inputs are provided most of the time.
            xRangeType = 'sym';
            if (strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz') || strcmp(obj.gridType,'AzAlt') || strcmp(obj.gridType,'RAdec') || strcmp(obj.gridType,'GalLongLat') )
                if min(obj.x) >= 0 && obj.symXZ == 0
                    xRangeType = 'pos';
                end
                if max(obj.y) - min(obj.y) <= pi+median(diff(unique(obj.y)))/2
                    yRangeType = '180';
                else
                    yRangeType = '360';
                end
            else
                yRangeType = [];
            end
        end
        
        function objNew = shiftRedun(obj,iout,iin,xAdd,yAdd)
            % SHIFTREDUN Shift redundant points.
            
            objNew = obj;
            Nredun = min(numel(iout),numel(iin));    % How many we have to remove/add
            
            % Remove old redundant points
            if Nredun > 0
                objNew.x(iout(1:Nredun)) = [];
                objNew.y(iout(1:Nredun)) = [];
                objNew.E1(iout(1:Nredun),:) = [];
                objNew.E2(iout(1:Nredun),:) = [];
            end
            
            % Add in new redundant points
            E1Add = obj.E1(iin,:);
            E2Add = obj.E2(iin,:);
            objNew.x = [objNew.x;xAdd(1:Nredun)];
            objNew.y = [objNew.y;yAdd(1:Nredun)];
            objNew.E1 = [objNew.E1;E1Add(1:Nredun,:)];
            objNew.E2 = [objNew.E2;E2Add(1:Nredun,:)];
        end
        
        
    end
    
end