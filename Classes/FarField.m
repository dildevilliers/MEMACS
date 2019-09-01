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
        coorType    % Coordinate system type {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3','power'}
        polType     % polarization type {'linear','circular','slant'}
        gridType    % Grid type {'PhTh','DirCos','AzEl','ElAz','Horiz','TrueView','ArcSin','Mollweide','RAdec','GalLongLat'}
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
        NxBase      % Number of unique x points in the base grid
        NyBase      % Number of unique y points in the base grid
        xRange  % [1x2] vector of the xRange
        yRange  % [1x2] vector of the yRange
        symXZ   % XZ plane symmetry
        symYZ   % YZ plane symmetry
        symXY   % XY plane symmetry
        julDate % Julian Date
        sphereGrids % Set of all spherical grids
        baseTypeDifferent % Flag that indicates base grid and current grid are local and astro (or changed astro) (true) or both local/both astro (false)
        coorOrientation % Coordinate system indicating the orientation
        angBackRotate   % Required Euler rotation angles to align the antenna to to North,Zenith system
    end
    
    properties (SetAccess = private, Hidden = true)
        E3 = []         % Radial E-field component
        time(1,1) datetime                              % Time as a datetime object (used, for instance, in astronomical observation)
        orientation(1,3) double {mustBeReal, mustBeFinite}      % Antenna orientation - in Euler angles (alpha,beta,gamma) in radians relative to an x->North, y->West, z->Zenith system  ([0,0,0])
        earthLocation(1,3)  double {mustBeReal, mustBeFinite}   % Antenna location on the Earth longitude and latitude in radians, and height above sea level in meters - defaults to the roof of the Stellenbosch University E&E Engineering Dept. :)

        % Keep the input data here to not lose some info when going through
        % transformations and back...
        xBase
        yBase
        phBase
        thBase
        gridTypeBase
        E1Base
        E2Base
        E3Base
        coorTypeBase
        polTypeBase
    end
    
    properties (Constant = true, Hidden = true)
        c0 = 299792458;
        eps0 = 8.854187817000001e-12;
        mu0 = 1.256637061435917e-06;
        eta0 = 3.767303134749689e+02;
        nSigDig = 3;
        projectionGrids = {'DirCos','TrueView','Arcsin','Mollweide'};
        astroGrids = {'Horiz','RAdec','GalLongLat'};
        localGrids = {'PhTh','AzEl','ElAz'};
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
            % - freq: Frequencies where the fields are defined in Hz (1e9), [1 x Nf]
            % - Prad: Radiated power at each frequency in W (4*pi/(2.*obj.eta0)), [1 x Nf]
            % - radEff: Radiation efficiency at each frequency (1), [1 x Nf]
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   keywords and values are from the sets 
            %   -- E3:          Radial component of E-field ([]), [Nang x 1]
            %   -- coorType:    {('spherical')|'Ludwig1'|'Ludwig2AE'|'Ludwig2EA'|'Ludwig3'}
            %   -- polType:     {('linear')|'circular'|'slant'}
            %   -- gridType:    {('PhTh')|'DirCos'|'AzEl'|'ElAz'|'Horiz'|'TrueView'|'ArcSin'|'Mollweide'|'RAdec'|'GalLongLat'}
            %   -- freqUnit:    {('Hz'),'kHz','MHz','GHz','THz'}
            %   -- symmetryXZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryYZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryXY:  {('none')|'electric'|'magnetic'}
            %   -- symBOR:      {('none')|'BOR0'|'BOR1'}
            %   -- r:           Radius where the E-Field is evaluated in m, (1)
            %   -- slant:       slant angle in rad for polType=slant, (pi/4)
            %   -- orientation: Antenna orientation - in Euler angles (alpha,beta,gamma) in radians relative to an x->North, y->West, z->Zenith system  ([0,0,0])
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
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   F.plot
            
            % Set up defaults: z-directed incremental dipole
            [ph0,th0] = PhThGrid;
            f0 = 1e9;
            Eth0 = sqrt(3*obj.eta0).*sin(th0);
            Eph0 = zeros(size(Eth0));
            % Limit the zeros to a finite number
            Eth0(Eth0 == 0) = lin20(-200);
            Eph0(Eph0 == 0) = lin20(-200);
            Prad0 = [];
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'FarField';
            
            % Optional inputs
            typeValidation_grid = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','ncols',1},'FarField');
            parseobj.addOptional('x',ph0,typeValidation_grid);
            parseobj.addOptional('y',th0,typeValidation_grid);
            
%             typeValidation_fields = @(x) validateattributes(x,{'numeric'},{'finite','nonnan'},'FarField');
            typeValidation_fields = @(x) validateattributes(x,{'numeric'},{},'FarField');
            parseobj.addOptional('E1',Eth0,typeValidation_fields);
            parseobj.addOptional('E2',Eph0,typeValidation_fields);
            
            typeValidation_freq = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nondecreasing','nrows',1},'FarField');
            parseobj.addOptional('freq',f0,typeValidation_freq);
            
%             typeValidation_power = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nrows',1},'FarField');
            typeValidation_power = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan'},'FarField');
            parseobj.addOptional('Prad',Prad0,typeValidation_power);
            parseobj.addOptional('radEff',1,typeValidation_power);
            
            % Name-value pairs
            parseobj.addParameter('E3',[],typeValidation_fields);
            
            expected_coorType = {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3','power'};
            parseobj.addParameter('coorType','spherical', @(x) any(validatestring(x,expected_coorType)));
            
            expected_polType = {'linear','circular','slant','none'};
            parseobj.addParameter('polType','linear', @(x) any(validatestring(x,expected_polType)));
            
            expected_gridType = {'PhTh','DirCos','AzEl','ElAz','Horiz','TrueView','ArcSin','Mollweide','RAdec','GalLongLat'};
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
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'FarField');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
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
            
            % Check input sizes
            Nang = size(obj.x,1);
            Nf = numel(obj.freq);
            assert(size(obj.y,1)==Nang && size(obj.E1,1)==Nang,'x,y and E1 must have the same number of rows')
            if ~strcmp(obj.coorType,'power')
                assert(size(obj.y,1)==Nang && size(obj.E1,1)==Nang && size(obj.E2,1)==Nang,'x,y,E1 and E2 must have the same number of rows')
            end
            % Integrate power if none provided, but a full sphere of fields is 
            if isempty(obj.Prad)
                if obj.isGrid4pi
                    obj.Prad = obj.pradInt;
                else
                    obj.Prad = 4*pi/(2.*obj.eta0);
                end
            end   
            if isscalar(obj.Prad) == 1, obj.Prad = repmat(obj.Prad,1,Nf); end
            if isscalar(obj.radEff) == 1, obj.radEff = repmat(obj.radEff,1,Nf); end
            assert(size(obj.E1,2)==Nf && size(obj.Prad,2)==Nf && size(obj.radEff,2)==Nf,'E1,freq, Prad and radEff must have the same number of columns')
            if ~strcmp(obj.coorType,'power')
                assert(size(obj.E1,2)==Nf && size(obj.E2,2)==Nf && size(obj.Prad,2)==Nf && size(obj.radEff,2)==Nf,'E1, E2, freq, Prad and radEff must have the same number of columns')
            end
            % Pol type meaningless for real, power only patterns
            if strcmp(obj.coorType,'power')
%                 assert(isreal(obj.E1),'Real values expected for E1 for power only pattern');
                obj.E2 = [];
                obj.polType = 'none'; 
            end
            
            tol = 10^(-obj.nSigDig);
            % Check input symmetry validity
            if sum(abs([obj.symXZ,obj.symYZ,obj.symXY])) > 0 && ~strcmp(obj.symmetryBOR,'none')
                error('Cannot specify a plane- and BOR-symmetry for the same field')
            end
            if ~strcmp(obj.symmetryXZ,'none')
                [~,v] = getDirCos(obj);
                assert(all(sign(v+tol) > 0) || all(sign(v-tol) < 0),'Invalid range for XZ symmetry')
            end
            if ~strcmp(obj.symmetryYZ,'none')
                [u] = getDirCos(obj);
                assert(all(sign(u+tol) > 0) || all(sign(u-tol) < 0),'Invalid range for YZ symmetry')
            end
            if ~strcmp(obj.symmetryXY,'none')
                error('function: setSymmetryXY not implemented yet - please redefine with the full field and symmetryXY = none');
            end
            if strcmp(obj.symmetryBOR,'BOR0')
                assert(any(obj.gridType,'PhTh'),'BOR0 symmetry only defined for PhTh grids here')
                assert(obj.Nx == 1,'Invalid range for BOR0 symmetry (one x-cut maximum)')
            elseif strcmp(obj.symmetryBOR,'BOR1')
                assert(strcmp(obj.gridType,'PhTh'),'BOR1 symmetry only defined for PhTh grids at this stage')
                assert(all(abs(unique(obj.ph) - [0;pi/2]) < tol),'Invalid range for BOR1 symmetry (E-plane and H-plane required)')
            end
            
        end
        
        %% Dependency-based Setters
        function ph = get.ph(obj)
            if ~any(strcmp(obj.gridType,obj.astroGrids))
                [ph,~] = getPhTh(obj);
            else
                ph = [];
            end
        end
        
        function th = get.th(obj)
            if ~any(strcmp(obj.gridType,obj.astroGrids))
                [~,th] = getPhTh(obj);
            else
                th = [];
            end
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
        
        function NxBase = get.NxBase(obj)
            NxBase = length(unique(obj.xBase));
        end
        
        function NyBase = get.NyBase(obj)
            NyBase = length(unique(obj.yBase));
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
        
        function xRange = get.xRange(obj)
            xVect = unique(obj.x);
            xRange = [min(xVect),max(xVect)];
        end
        
        function yRange = get.yRange(obj)
            yVect = unique(obj.y);
            yRange = [min(yVect),max(yVect)];
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
        
        function julDate = get.julDate(obj)
            julDate = convert.date2jd([obj.time.Day,obj.time.Month,obj.time.Year,obj.time.Hour,obj.time.Minute,obj.time.Second]);
        end
        
        function sphereGrids = get.sphereGrids(obj)
            sphereGrids = [obj.astroGrids,obj.localGrids];
        end
        
        function baseTypeDifferent = get.baseTypeDifferent(obj)
            baseTypeDifferent = (any(strcmp(obj.gridType,obj.astroGrids)) && any(strcmp(obj.gridTypeBase,obj.localGrids))) || ...
                    (any(strcmp(obj.gridTypeBase,obj.astroGrids)) && any(strcmp(obj.gridType,obj.localGrids))) || ...
                    (any(strcmp(obj.gridType,obj.astroGrids)) && any(strcmp(obj.gridTypeBase,obj.astroGrids)) && ~strcmp(obj.gridType,obj.gridTypeBase));
        end
        
        function coorOrientation = get.coorOrientation(obj)
            C0 = CoordinateSystem;
            coorOrientation = C0.rotEuler(obj.orientation);
        end
        
        function angBackRotate = get.angBackRotate(obj)
            angBackRotate = getEulerangBetweenCoors(CoordinateSystem,obj.coorOrientation);
        end
        
        %% Field and frequency setters
        function obj = setFreq(obj,freq,freqUnit)
            % SETFREQ sets the frequency of the object
            
            if nargin > 1
                assert(numel(freq) == size(obj.E1,2),'Error, freq must be the same length as the number of columns in E1')
                obj.freq = freq;
            end
            if nargin > 2
                obj.freqUnit = freqUnit;
            end
        end
        
        function obj = setEfield(obj,iSet,E1in,E2in,E3in,calcP)
            % SETEFIELD sets the E-fields of the object
            
            % Sets the E-fields after creation of the object at the
            % positions iSet.  The size of iSet must correspond the the
            % size of obj.E1 = [Nang x Nf], or be a single column of length
            % Nang (changes will be applied to all frequencies then).  Use
            % with caution as it can change the power in the pattern and 
            % the associated directivities.
            % It changes the base and tries to calculate a new power if the
            % pattern exists over a full sphere
            
            assert(size(iSet,1) == obj.Nang,'iSet must have Nang rows')
            assert(size(iSet,2) == obj.Nf || size(iSet,2) == 1,'iSet must have Nf or 1 columns')
            
            if nargin < 6, calcP = true; end
            if size(iSet,2) == 1, iSet = repmat(iSet,1,obj.Nf); end
            if nargin > 2, obj.E1(iSet) = E1in; end
            if nargin > 3 && ~isempty(E2in)
                obj.E2(iSet) = E2in; 
            end
            if nargin > 4 && ~isempty(E3in)
                obj.E3(iSet) = E3in;
            end
            if ~isempty(obj.xBase) || ~isempty(obj.E1Base)
                warning('FarField:baseRemoveWarning','setEfield will remove the base grid from the object, since the operation is performed on the current grid, and the base typically becomes invalid')
            end
            obj = obj.clearBase;
            
            if obj.isGrid4pi && calcP
                obj.Prad = obj.pradInt; 
            end
        end
        
        %% Pattern getters
        function FFpattern = getFarFieldStruct(obj)
            % GETFARFIELDSTRUCT Returns the legacy FarField struct data format.
            % FFpattern = getFarFieldStruct(obj) Converts a FarField
            % object into a struct data type. The struct format is compatible
            % with old script-based code.  In some cases the th and ph
            % fields must be expanded to Nf columns, but this depends on
            % the MATLAB version and the script implementation.
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
            % Updated: 2019-08-19, Dirk de Villiers
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
            FFpattern.th = obj.y;
            FFpattern.ph = obj.x;
            FFpattern.Eth = obj.E1;
            FFpattern.Eph = obj.E2;
            FFpattern.freq = obj.freqHz;
            FFpattern.Nth = obj.Ny;
            FFpattern.Nph = obj.Nx;
            FFpattern.Nf = obj.Nf;
            FFpattern.Prad = obj.Prad;
            FFpattern.radEff = obj.radEff;
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
            
            if strcmp(obj.coorType,'power')
                E1field = obj.E1;
                [E2field,E3field] = deal([]);
            else
                k = 2.*pi.*obj.freqHz./obj.c0;
                FFfact = exp(-1i.*k.*obj.r)./obj.r;
                E1field = bsxfun(@times,obj.E1,FFfact);
                if ~isempty(obj.E2)
                    E2field = bsxfun(@times,obj.E2,FFfact);
                else
                    E2field = [];
                end
                if ~isempty(obj.E3)
                    E3field = bsxfun(@times,obj.E3,FFfact);
                else
                    E3field = [];
                end
            end
        end
        
        function [W] = getW(obj)
            % GETW Returns the radiation density in W/m2.
            
            % function [W] = getW(obj)
            % returns the radiation density in W/m2 [Nang x Nf]
            [E1f, E2f] = getEfield(obj);    % Can use any orthogonal pair
            if ~isempty(E2f)
                W = 1./(2.*obj.eta0).*(abs(E1f).^2 + abs(E2f).^2);
            else
                W = 1./(2.*obj.eta0).*(abs(E1f).^2);
            end
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
            assert(~strcmp(obj.coorType,'power'),'AR undefined for power only patterns')
            AR = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        function [ARinv] = getAxialRatioInv(obj)
            % GETAXIALRATIOINV Returns the inverted Axial Ration (linear).
            
            % function [ARinv] = getAxialRatioInv(obj)
            % returns the inverted Axial Ratio in ARinv [Nang x Nf]
            assert(~strcmp(obj.coorType,'power'),'AR undefined for power only patterns')
            ARinv = sqrt((abs(obj.E1).^2 + abs(obj.E2).^2 - abs(obj.E1.^2 + obj.E2.^2))./(abs(obj.E1).^2 + abs(obj.E2).^2 + abs(obj.E1.^2 + obj.E2.^2)));
        end
        
        function [Xpol] = getCO_XP(obj)
            % GETCO_XP Returns the CO/XP ratio (linear).
            
            % function [Xpol] = getCO_XP(obj)
            % returns the CO/XP ratio (linear) [Nang x Nf]
            assert(~strcmp(obj.coorType,'power'),'Polarisation undefined for power only patterns')
            Xpol = (abs(obj.E2)./abs(obj.E1)).^2;
        end
        
        function [Xpol] = getXP_CO(obj)
            % GETXP_CO Returns the XP/CO ratio (linear)
            
            % function [Xpol] = getXP_CO(obj)
            % returns the XP/CO ratio (linear) [Nang x Nf]
            assert(~strcmp(obj.coorType,'power'),'Polarisation undefined for power only patterns')
            Xpol = (abs(obj.E1)./abs(obj.E2)).^2;
        end
        
        %% Performance metrics
        function [SLL1,SLL2,SLLstruct] = getSLL(obj)
            % GETSLL Get the sidelobe level of the beam
            %
            % [SLL1,SLL2,SLLstruct] = getSLL(obj) calculates the the
            % sidelone level of a well-defined beam pattern pointing in the
            % z-direction.  
            % 
            % Inputs
            % - obj: FarField object (defined on regular grid)
            %
            % Outputs
            % - SLL1:  Maximum first SLL (over x-grid) in dB ([1 x obj.Nf])
            % - SLL2:  Maximum second SLL (over x-grid) in dB ([1 x obj.Nf])
            % - SLLstruct: Full information over cut angles and frequency ([obj.Nx x obj.Nf])
            %   -- SLL1: first SLL (in dB)
            %   -- SLL2: second SLL (in dB)
            %   -- ang1: angle (y-value) of first SLL (in unit of obj.y)
            %   -- ang2: angle (y-value) of second SLL (in unit of obj.y)
            %
            % Dependencies
            % -
            %
            % Created: 2019-08-17, Dirk de Villiers
            % Updated: 2019-08-18, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField.readGRASPgrd;
            %   [SLL1,SLL2] = F.getSLL;
            %   iF = 1;
            %   F.plot('plotType','2D','showGrid',1,'norm',true,'freqIndex',iF), axis normal, hold on, view([140,45])
            %   plot3(unique(rad2deg(F.x)),rad2deg(sllStruct.ang1(:,iF)),sllStruct.SLL1(:,iF),'or')
            %   plot3(unique(rad2deg(F.x)),rad2deg(sllStruct.ang2(:,iF)),sllStruct.SLL2(:,iF),'or')

            
            assert(obj.isGridUniform,'Uniform grid expected for SLL calculation')
            
            [SLL1,SLL2] = deal(NaN(1,obj.Nf));
            [SLL1full,SLL2full,ang1full,ang2full] = deal(NaN(obj.Nx,obj.Nf)); 
            % Get valid th angles from the -3dB beamwidth
            [~,BW] = obj.getBeamwidth(-3);
            BWmin = min(BW);     % Minimum -3dB beamwidth
            for ff = 1:obj.Nf
                % Get normalized dB directivity
                D = dB10(obj.getFi(ff).getDirectivity) - obj.Directivity_dBi(ff);
                
                yVect = unique(obj.y);
                AF = [obj.y(1:obj.Ny),reshape(D,obj.Ny,obj.Nx)];
                AF = AF(yVect > BWmin(ff),:);
                
                try
                    % Not the most stable function - so catch errors and
                    % return nothing
                    [pk,ps] = SLL(AF);
                catch 
                    % Just return the NaNs
                    return;
                end
                SLL1(ff) = max(pk(:,1));
                if nargout > 1
                    SLL2(ff) = max(pk(:,2));
                end
                if nargout > 2
                   SLL1full(:,ff) = pk(:,1);
                   SLL2full(:,ff) = pk(:,2);
                   ang1full(:,ff) = ps(:,1);
                   ang2full(:,ff) = ps(:,2);
                end
            end
            if nargout > 2
                SLLstruct = struct('SLL1',SLL1full,'SLL2',SLL2full,'ang1',ang1full,'ang2',ang2full);
            end
        end
        
        function [BW,BWfull] = getBeamwidth(obj,dBlevel)
            % GETBEAMWIDTH Get the beamwidth of the beam
            %
            % [BW,BWfull] = getBeamwidth(obj,dBlevel) calculates the the
            % beamwidth of a well-defined beam pattern pointing in the
            % z-direction.  A vector of requested beamwidth levels can be provided
            % 
            % Inputs
            % - obj: FarField object (defined on regular grid)
            % - dBlevel: The (vector) value of the beam taper where the
            %            beamwidth is required (negative dB). A value of
            %            -inf is interpreted as requesting the first null.
            %            Default -3 dB
            %
            % Outputs
            % - BW:  Mean beamwidth (over x-grid) in the unit of the y-grid ([length(dBlevel) x obj.Nf])
            % - BWfull: Full beamwidth over cut angles and frequency ([obj.Nx x obj.Nf x length(dBlevel)])
            %
            % Dependencies
            % -
            %
            % Created: 2019-08-18, Dirk de Villiers
            % Updated: 2019-08-18, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField.readGRASPgrd;
            %   [BW,BWfull] = F.getBeamwidth([-3,-inf]);
            %   iF = 3;
            %   F.plot('plotType','2D','showGrid',1,'norm',true,'freqIndex',iF), axis normal, hold on
            %   plot(rad2deg(unique(F.x)),rad2deg(BWfull(:,iF,1)),'k','linewidth',2)
            %   plot(rad2deg(unique(F.x)),rad2deg(BWfull(:,iF,2)),'k','linewidth',2)
            
            % TODO: get more accurate first null position
            
            assert(obj.isGridUniform,'Uniform grid expected for SLL calculation')
            
            if nargin < 2
                dBlevel = -3;
            else
                if any(dBlevel >= 0)
                    warning('dBlevel should typically be a negative number')
                end
            end
            
            BWfull = NaN(obj.Nx,obj.Nf,length(dBlevel));
            BW = NaN(length(dBlevel),obj.Nf);
            nullPos = obj.yRange(2); % Deafault the null position to the edge of the y-range
            try  % Can sometimes crash at a variety of places for strange beams - just return NaNs
                for ff = 1:obj.Nf
                    % Get normalized dB directivity
                    D = reshape(dB10(obj.getFi(ff).getDirectivity),obj.Ny,obj.Nx) - obj.Directivity_dBi(ff);
                    for xx = 1:obj.Nx
                        % Get a single cut
                        dvect = D(:,xx);
                        % Find the first null - the difference in pattern must be positive
                        % after a negative difference...
                        ddif = diff(sign(diff(dvect)));
                        nP = find(ddif == 2,1);
                        if ~isempty(nP), nullPos = nP; end
                        % Interpolate up to the null (ang as function of level)
                        angVect = obj.y(1:nullPos);
                        for dd = 1:length(dBlevel)
                            if isinf(dBlevel(dd))
                                BWfull(xx,ff,dd) = angVect(nullPos);
                            else
                                
                                BWfull(xx,ff,dd) = interp1(dvect(1:nullPos),angVect,dBlevel(dd),'spline');
                                
                            end
                        end
                    end
                end
                for dd = 1:length(dBlevel)
                    BW(dd,:) = mean(BWfull(:,:,dd));
                end
            catch
                return
            end
        end
        
        %% Field normalization
        function P = pradInt(obj)
            % PRADINT  Calculates the total power in the field
            %
            % P = pradInt(obj) calculates the total power in the field over
            % the full available grid.  Will not return the actual radiated
            % power if the grid does not cover the full sphere.
            % 
            % Inputs
            % - obj: FarField object
            %
            % Outputs
            % - P:  Radiated power in Watt
            %
            % Dependencies
            % -
            %
            % Created: 2019, Dirk de Villiers
            % Updated: 2019-08-09, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   P = F.pradInt
            
            symFact = 2^(sum(abs([obj.symXY,obj.symXZ,obj.symYZ])));
            
            % Can be slow for huge objects
            ny = obj.Ny;
            nx = obj.Nx;
            assert(ny*nx == obj.Nang,'Can only integrate power on regular, plaid, monotonic grids. Please provide the power with the constructor when making this object.')
            
            switch obj.gridType
                case 'PhTh'
                        JacFuncTh = @sin;
                case {'AzEl','ElAz'} 
                        JacFuncTh = @cos;
                otherwise
                    error(['pradInt not implemented for gridType = ',obj.gridType])
            end
            PH = reshape(obj.x,ny,nx);
            TH = reshape(obj.y,ny,nx);
            U = obj.getU;
            P = zeros(1,obj.Nf);
            for ff = 1:obj.Nf
                if strcmp(obj.symmetryBOR,'none')
                    integrand = reshape(U(:,ff),ny,nx).*abs(JacFuncTh(TH));
                    P(ff) = integral2D(PH,TH,integrand);
                else
                    Nth = ny;
                    th_vect = obj.y(1:Nth);
                    if strcmp(obj.symmetryBOR,'BOR0')
                        integrand = 2*U(:,ff).*sin(th_vect);
                    elseif strcmp(obj.symmetryBOR,'BOR1')
                        integrand = (U(1:Nth,ff) + U(Nth+1:end,ff)).*abs(JacFuncTh(th_vect));
                    end
                    P(ff) = pi*integral1D(th_vect,integrand);
                    symFact = 1;    % Just to be sure...
                end
            end
            P = P.*symFact;
        end
        
        function obj = setPower(obj,powerWatt)
            % SETPOWER Normalizes the FarField object to power level
            %
            % obj = setPower(obj,powerWatt) Normalizes the FarField object 
            % to have the a total radiated power specified in powerWatt
            % The field need not be specified over the full sphere - the
            % total intercepted power in the specified sector will be set
            % to powerWatt.  
            % 
            % Inputs
            % - obj: FarField object
            % - powerWatt: desired radiated power in Watt scalar or 
            %              [1 x obj.Nf] (default = 4*pi/(2*377) W)
            %
            % Outputs
            % - obj: FarField object
            %
            % Dependencies
            % -
            %
            % Created: 2019, Dirk de Villiers
            % Updated: 2019-08-09, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   F = F.setPower(1);
            %   P = F.pradInt
            
            if nargin == 1
                powerWatt = 4*pi/(2.*obj.eta0);
            end
            if length(powerWatt) == 1
                powerWatt = repmat(powerWatt,1,obj.Nf);
            elseif length(powerWatt) ~= obj.Nf
                error('powerWatt should be scalar or of length obj1.Nf');
            end
            P = obj.Prad;
            Cn = powerWatt./(P);
            obj.Prad = P;
            obj = scale(obj,sqrt(Cn));
        end
        
        %% Grid transformation setters
        function obj = changeGrid(obj,gridTypeString)
            % CHANGEGRID Change the current FarField object grid.  
            %
            % obj = changeGrid(obj,gridTypeString) transforms the current
            % grid to that specified in gridTypeString.
            % This function calls the appropriate grid2* function to
            % excecute. Note that all the grid2* functions can be directly
            % called with the object as argument to do the same job.
            % 
            % In all cases, if no base grid exists in the object
            % the base will be set before the transformation to maintain
            % the original provided grid values for later use. If a base
            % grid is present, the object will first be set into the base
            % grid form before transformation - transformations thus always
            % happens using the original provided data in the base.
            %
            % If the object is already in the specified grid, nothing is
            % done.
            %
            % When a projection type grid is provided, and transformation
            % is done to one of the astronomical grids, the transformation
            % is done locally.  That means no orientation, position or time
            % information is used - it is assumed that the provided
            % projection is of the astronomical grid.  The opposite is also
            % true - astronomical grids will be projected locally and not
            % shifted and rotated first.
            % 
            % Inputs
            % - obj: FarField object
            % - gridTypeString: 'PhTh'|'DirCos'|'AzEl'|'ElAz'|'TrueView'|'ArcSin'|'Horiz'|'RAdec'|'GalLongLat'
            %
            % Outputs
            % - obj: FarField object - possibly with added base grid
            %
            % Dependencies
            % -
            %
            % Created: 2019, Ridalise Louw
            % Updated: 2019-08-13, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   F = F.changeGrid('AzEl');
            %   F.plot('plotType','2D','showGrid',1)
            
            mustBeMember(gridTypeString, {'PhTh','DirCos','AzEl','ElAz','TrueView','ArcSin','Horiz','RAdec','GalLongLat'});
            handleGridType = str2func(['grid2',gridTypeString]);
            obj = handleGridType(obj);  
        end
        
        % Local grids
        function obj = grid2PhTh(obj)
            % GRID2PHTH Change the current grid to a PhTh grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'PhTh')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                mustRotate = ~all(obj.orientation == [0,0,0]) && any(strcmp(obj.gridType,obj.astroGrids));
                [obj.x,obj.y] = getPhTh(obj);
                obj.gridType = 'PhTh';
                % Rotate to sort out orientation
                if mustRotate, obj = obj.rotate(@rotEulersph,obj.angBackRotate); end
            end
        end
        
        function obj = grid2AzEl(obj)
            % GRID2AZEL Change the current grid to a AzEl grid.
            %
            % See help changeGrid for details
            
            % First get in PhTh from the base, then transform the current
            % system to AzEl.  This way the rotation is done on the PhTh
            % grid only...
            
            if ~strcmp(obj.gridType,'AzEl')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                mustRotate = ~all(obj.orientation == [0,0,0]) && any(strcmp(obj.gridType,obj.astroGrids));
                if mustRotate, obj = obj.grid2PhTh; end
                [obj.x,obj.y] = getAzEl(obj);
                obj.gridType = 'AzEl';
            end
        end
        
        function obj = grid2ElAz(obj)
            % GRID2ELAZ Change the current grid to a ElAz grid.
            %
            % See help changeGrid for details
                        
            % First get in PhTh from the base, then transform the current
            % system to AzEl.  This way the rotation is done on the PhTh
            % grid only...
            
            if ~strcmp(obj.gridType,'ElAz')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                mustRotate = ~all(obj.orientation == [0,0,0]) && any(strcmp(obj.gridType,obj.astroGrids));
                if mustRotate, obj = obj.grid2PhTh; end
                [obj.x,obj.y] = getElAz(obj);
                obj.gridType = 'ElAz';
            end
        end
        
        % Projections
        function obj = grid2DirCos(obj)
            % GRID2DIRCOS Change the current grid to a DirCos grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'DirCos')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                [obj.x,obj.y] = getDirCos(obj);
                obj.gridType = 'DirCos';
            end
        end
        
        function obj = grid2TrueView(obj)
            % GRID2TRUEVIEW Change the current grid to a TrueView grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'TrueView')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                [obj.x,obj.y] = getTrueView(obj);
                obj.gridType = 'TrueView';
            end
        end
        
        function obj = grid2ArcSin(obj)
            % GRID2ARCSIN Change the current grid to an ArcSin grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'ArcSin')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                [obj.x,obj.y] = getArcSin(obj);
                obj.gridType = 'ArcSin';
            end
        end
        
        % Astro grids
        function obj = grid2Horiz(obj)
            % GRID2Horiz Change the current FarField object grid to an Horiz grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'Horiz')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                obj = obj.coor2power(false);
                mustRotate = ~all(obj.orientation == [0,0,0]) && any(strcmp(obj.gridType,obj.localGrids));
                if mustRotate, obj = obj.rotate(@rotEulersph,obj.orientation); end
                [obj.x,obj.y] = getHoriz(obj);
                obj.gridType = 'Horiz';
            end
        end
        
        function obj = grid2RAdec(obj)
            % GRID2RAdec Change the current FarField object grid to a RAdec grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'RAdec')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                obj = obj.coor2power(false);
                mustRotate = ~all(obj.orientation == [0,0,0]) && any(strcmp(obj.gridType,obj.localGrids));
                if mustRotate, obj = obj.rotate(@rotEulersph,obj.orientation); end
                [obj.x,obj.y] = getRAdec(obj);
                obj.gridType = 'RAdec';
            end
        end
        
        function obj = grid2GalLongLat(obj)
            % GRID2GALLONGLAT Change the current FarField object grid to a GalLongLat grid.
            %
            % See help changeGrid for details
            
            if ~strcmp(obj.gridType,'GalLongLat')
                if isempty(obj.xBase)
                    obj = obj.setBaseGrid;
                else
                    obj = obj.grid2Base;
                end
                obj = obj.coor2power(false);
                mustRotate = ~all(obj.orientation == [0,0,0]) && any(strcmp(obj.gridType,obj.localGrids));
                if mustRotate, obj = obj.rotate(@rotEulersph,obj.orientation); end
                [obj.x,obj.y] = getGalLongLat(obj);
                obj.gridType = 'GalLongLat';
            end
        end
        
        %% Grid range shifters
        function obj = sortGrid(obj,nSigDig)
            % SORTGRID Sort grid with corresponding E-field values 
            % in ascending order, according to x then y.
            
            % Round to some significant digits for sorting (some issues can
            % arise in deg2rad and rad2deg
            if nargin < 2
                nSigDig = obj.nSigDig;
            end
            obj = roundGrid(obj,nSigDig);
%             [~,iSort] = sortrows([obj.x,obj.y],[1 2]);
            [~,iSort] = unique([obj.x,obj.y],'rows');
            obj.x = obj.x(iSort);
            obj.y = obj.y(iSort);
            obj.E1 = obj.E1(iSort,:);
            if ~isempty(obj.E2), obj.E2 = obj.E2(iSort,:); end
            if ~isempty(obj.E3), obj.E3 = obj.E3(iSort,:); end
            
            % make sure the base grid keeps the same order in order to associate fields and base values properly
            if ~isempty(obj.xBase)
                obj.xBase = obj.xBase(iSort);
                obj.yBase = obj.yBase(iSort);
                [obj.phBase,obj.thBase] = obj.getPhTh; % Dependent, so no need to sort
            end
            if ~isempty(obj.E1Base)
                obj.E1Base = obj.E1Base(iSort,:);
                if ~isempty(obj.E2Base), obj.E2Base = obj.E2Base(iSort,:); end
                if ~isempty(obj.E3Base), obj.E3Base = obj.E3Base(iSort,:); end
            end
        end
        
        function obj = roundGrid(obj,nSigDig)
            % ROUNDGRID Round grid entries to some significant number of digits.
            
            % Round to some significant digits for sorting (some issues can
            % arise in deg2rad and rad2deg
            if nargin < 2
                nSigDig = obj.nSigDig;
            end
%             xRound = round(obj.x*10^nSigDig)/10^nSigDig;
%             yRound = round(obj.y*10^nSigDig)/10^nSigDig;
            xRound = round(rad2deg(obj.x)*10^nSigDig)/10^nSigDig;
            yRound = round(rad2deg(obj.y)*10^nSigDig)/10^nSigDig;
            obj.x = deg2rad(xRound);
            obj.y = deg2rad(yRound);
        end
        
        function obj = copyAndInsertXcut(obj1,xvalCopy,xvalPaste,tol)
            % COPYANDINSERTXCUT Copy a FarField x-cut into another position.
            
            % Use this to copy an X cut into another position.  Typically
            % handy when some transformation does not include the closing
            % cut - that is the 0 and 360 or -180 and 180 cuts.  Can in
            % principle be used to do random stuff - so careful.
            
            if nargin < 4
                tol = mean(diff(unique(obj1.x)));
            end
            % Make a whole new object to initialise the base
            % correctly - that is, no base after this change since it
            % changed the grid size
            inInd = find(abs(obj1.x - xvalCopy) < tol);
            xNew = [obj1.x;xvalPaste.*ones(size(inInd))];
            yNew = [obj1.y;obj1.y(inInd)];
            E1New = [obj1.E1;obj1.E1(inInd,:)];
            if ~isempty(obj1.E2)
                E2New = [obj1.E2;obj1.E2(inInd,:)];
            else
                E2New = [];
            end
            if ~isempty(obj1.E3)
                E3New = [obj1.E3;obj1.E3(inInd,:)];
            else
                E3New = [];
            end
            obj = FarField(xNew,yNew,E1New,E2New,obj1.freq,obj1.Prad,obj1.radEff,...
                'coorType',obj1.coorType,'polType',obj1.polType,'gridType',obj1.gridType,'freqUnit',obj1.freqUnit,'r',1,...
                'slant',obj1.slant,'orientation',obj1.orientation,'earthLocation',obj1.earthLocation,'time',obj1.time);
            obj = obj.sortGrid;
        end
        
        function obj = copyAndInsertYcut(obj1,yvalCopy,yvalPaste,tol)
            % COPYANDINSERTYCUT Copy a FarField y-cut into another position.
            
            % Use this to copy an X cut into another position.  Typically
            % handy when some transformation does not include the closing
            % cut - that is the 0 and 360 or -180 and 180 cuts.  Can in
            % principle be used to do random stuff - so careful.
            
            if nargin < 4
                tol = mean(diff(unique(obj1.y)));
            end
            % Make a whole new object to initialise the base
            % correctly - that is, no base after this change since it
            % changed the grid size
            inInd = find(abs(obj1.y - yvalCopy) < tol);
            xNew = [obj1.x;obj1.x(inInd)];
            yNew = [obj1.y;yvalPaste.*ones(size(inInd))];
            E1New = [obj1.E1;obj1.E1(inInd,:)];
            if ~isempty(obj1.E2)
                E2New = [obj1.E2;obj1.E2(inInd,:)];
            else
                E2New = [];
            end
            if ~isempty(obj1.E3)
                E3New = [obj1.E3;obj1.E3(inInd,:)];
            else
                E3New = [];
            end
            obj = FarField(xNew,yNew,E1New,E2New,obj1.freq,obj1.Prad,obj1.radEff,...
                'coorType',obj1.coorType,'polType',obj1.polType,'gridType',obj1.gridType,'freqUnit',obj1.freqUnit,'r',1,...
                'slant',obj1.slant,'orientation',obj1.orientation,'earthLocation',obj1.earthLocation,'time',obj1.time);
            obj = obj.sortGrid;
        end
        
        function obj = setYrange(obj,type)
            % SETYRANGE Set the y-range (th, el, or al) for the angular gridTypes
            % in the FarField object.
            % Function deprecated by setRangeSph. This version is unstable in some
            % cases and not being maintained - will be removed in future
            % release.
            
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
            
            warning('Method FarField.setYrange is deprecated by FarField.setRangeSph.  Please change code accordingly.  This function will be removed in a future release.')
            
            assert(type == 180 || type == 360,'Unknown type: Should be 180 or 360');
            % Do nothing if the range is already what is requested
            if strcmp(obj.yRangeType,num2str(type)), return; end
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
                        % Fix pole fields sign change
                        % TODO
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
                        % Fix pole fields sign change
                        % TODO
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
                        % Fix pole fields sign change
                        % TODO
                    elseif strcmp(obj.xRangeType,'sym')
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
                        % Fix pole fields sign change
                        iPole = intersect(find(abs(obj.x) >= pi/2-eps),find(abs(obj.y) < eps));
                        obj.E1(iPole,:) = -obj.E1(iPole,:);
                        obj.E2(iPole,:) = -obj.E2(iPole,:);
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
%                         obj = shiftRedun(obj,iout,iin,xAdd,yAdd);
                        obj = insertMissingCuts(obj,iin,xAdd,yAdd);
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
%             obj = roundGrid(obj);
            obj = obj.sortGrid;
        end
        
        function obj = setXrange(obj, type)
            % SETXRANGESPH Set the x-range (ph, az or ep) for the spherical gridTypes
            % Function deprecated by setRangeSph. This version is unstable in some
            % cases and not being maintained - will be removed in future
            % release.
            
            warning('Method FarField.setXrange is deprecated by FarField.setRangeSph.  Please change code accordingly.  This function will be removed in a future release.')
            
            mustBeMember(type,{'pos','sym'})
            tol = 10^(-obj.nSigDig+2);
            if any(strcmp(obj.gridType,obj.sphereGrids))
                if ~strcmp(type,obj.xRangeType) % Do nothing if not required
                    xSpanOrig = max(obj.x) - min(obj.x);
                    switch type
                        case 'sym'
                            switch obj.yRangeType
                                case '180'
                                    obj.x = wrap2pi(obj.x);
                                    [val1,val2] = deal(-pi,pi);
                                case '360'
                                    obj.y = wrap2pi(obj.y);
                                    iShift = find(obj.x > pi/2);
                                    xshift = @(x) x - pi;
                                    
                                    obj.x = wrap2pi(obj.x);
                                    [val1,val2] = deal(-pi,pi);
                            end
                        case 'pos'
                            obj.x = wrap22pi(obj.x);
                            [val1,val2] = deal(0,2*pi);
                    end
                    
                    % Insert missing cuts
                    if abs(xSpanOrig - 2*pi) < tol
                        i1 = find(abs(obj.x - val1) < tol);
                        i2 = find(abs(obj.x - val2) < tol);
                        iin = [i1;i2];
                        if numel(iin) > 0
                            xAdd = [ones(size(i1)).*val2;ones(size(i2)).*val1];
                            yAdd = obj.y(iin);
                            obj = insertMissingCuts(obj,iin,xAdd,yAdd);
                        end
                    end
                    % Sort
                    obj = obj.sortGrid(obj.nSigDig-1);
                end
            else
                warning('setXrange only operates on the spherical grid types')
            end
        end
        
        function obj = setRangeSph(obj,xType,yType)
            % SETRANGESPH sets the x and y ranges of spherical grids
            %
            % obj = setRangeSph(obj,xType,yType) sets the x- and y-ranges
            % of the object to xType = {('sym') | 'pos'}; yType = {('180') | '360'}
            % Depending on the gridType, different representations of the y
            % range result.  In some cases the base grid is changed -
            % specifically when the required transformation results in a
            % different number of grid points to accommodate the new poles.
            % This is avoided as far as possible, and only happens when the
            % full sphere is represented and when changing between yRange
            % 180 and 360 (and back);
            % 
            % Inputs
            % - obj: FarField object
            % - xType: char array {('sym') | 'pos'}
            % - yType: char array {('180') | '360'}
            %
            % Outputs
            % - obj:    Farfield object with new range and possible new base
            %
            % Dependencies
            % -
            %
            % Created: 2019-06-10, Dirk de Villiers
            % Updated: 2019-08-17, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   F.plot('plotType','2D','showGrid',true)
            %   Fs180 = F.setRangeSph;
            %   figure, Fs180.plot('plotType','2D','showGrid',true)
            %   Fp360 = F.setRangeSph('pos','360');
            %   figure, Fp360.plot('plotType','2D','showGrid',true)
            %   Fs360 = F.setRangeSph('sym','360');
            %   figure, Fs360.plot('plotType','2D','showGrid',true)
            
            
            % Strategy is to have functions to go to and from the standard
            % sym|180 case
            if nargin < 3 
                yType = '180';
            end
            if nargin < 2
                xType = 'sym';
            end
            if isnumeric(yType), yType = num2str(yType); end
            mustBeMember(xType,['sym','pos'])
            mustBeMember(yType,['180','360'])
            if any(strcmp(obj.gridType,obj.sphereGrids))
                rangeHandle = str2func(['range',yType,xType]);
                if ~isempty(obj.xBase) || ~isempty(obj.E1Base)
                    warning('FarField:baseRemoveWarning','setRangeSph will remove the base grid from the object, since the operation is performed on the current grid, and the grid size might change')
                end
                obj = obj.clearBase;
                obj = rangeHandle(obj);
            else
                warning(['Cant set the range of a non-spherical grid like: ',obj.gridType,'; setRangeSph did nothing'])
            end
        end
        
        function obj = getRange(obj,xRangeNew,yRangeNew)
            % GETRANGE returns a FarField over a reduced angular range
            %
            % obj = getRange(obj,xRangeNew,yRangeNew) reduces the range of
            % the object to that specified in xRangeNew and yRangeNew.
            % Operates on the current grid, and will return strange results
            % if ranges are specified as projections or angles if they
            % should not be.  The output object base is changed to the new
            % range specified.
            % 
            % Inputs
            % - obj: FarField object
            % - xRangeNew: 2 element vector with [xmin,xmax] in rad
            % - yRangeNew: 2 element vector with [xmin,xmax] in rad (if empty: obj.yRange)
            %
            % Outputs
            % - obj:    Farfield object with new base
            %
            % Dependencies
            % -
            %
            % Created: 2019-06-10, Dirk de Villiers
            % Updated: 2019-06-10, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   [xR,yR] = deal([pi/2,3*pi/2],[0,pi/2])
            %   Fn = F.getRange(xR,yR);
            %   Fn.plot('plotType','2D','step',5,'showGrid',true)
            
            if nargin < 3
                yRangeNew = obj.yRange;
            end
            assert(numel(xRangeNew)==2,'Expecting a 2 element vector for xRangeNew')
            assert(numel(yRangeNew)==2,'Expecting a 2 element vector for yRangeNew')
            
            tolX = diff(obj.xRange)./obj.Nx/10;
            tolY = diff(obj.yRange)./obj.Ny/10;
            iX = find(obj.x >= xRangeNew(1)-tolX & obj.x <= xRangeNew(2)+tolX);
            iY = find(obj.y >= yRangeNew(1)-tolY & obj.y <= yRangeNew(2)+tolY);
            iKeep = intersect(iX,iY);
            obj.x = obj.x(iKeep);
            obj.y = obj.y(iKeep);
            obj.E1 = obj.E1(iKeep,:);
            if numel(obj.E2) > 0
                obj.E2 = obj.E2(iKeep,:);
            end
            if numel(obj.E3) > 0
                obj.E3 = obj.E3(iKeep,:);
            end
            if ~isempty(obj.xBase)
                obj.xBase = obj.xBase(iKeep);
                obj.yBase = obj.yBase(iKeep);
            end
            if ~isempty(obj.E1Base)
                obj.E1Base = obj.E1Base(iKeep);
            end
            if ~isempty(obj.E2Base)
                obj.E2Base = obj.E2Base(iKeep);
            end
            if ~isempty(obj.E3Base)
                obj.E3Base = obj.E3Base(iKeep);
            end
        end
        
        %% Coordinate system transformation methods
        function obj = changeCoor(obj,coorTypeString,setStdGrid)
            % CHANGECOOR Change the FarField object coordinate type.
            %
            % obj = changeCoor(obj,coorTypeString,setStdGrid) transforms 
            % the current coordinate system to that specified in 
            % coorTypeString. The grid is also changed to the standard grid
            % associated with each coordinate system type if requested by
            % the setStdGrid flag.
            % This function calls the appropriate coor2* function to
            % excecute. Note that all the coor2* functions can be directly
            % called with the object as argument to do the same job.
            % 
            % In all cases (except power), if no base field exists in the object
            % the base will be set before the transformation to maintain
            % the original provided field values for later use. If a base
            % field is present, the object will first be set into the base
            % field form before transformation - transformations thus always
            % happens using the original provided data in the base.
            %
            % If the object is already in the specified coordinate type, nothing is
            % done.
            %
            % If changing to the power type, the base field will be
            % cleared. The idea is to store a field with smaller memory
            % requirements, so only one component (magnitude) is kept, and
            % therefore transforming back to another coorType from power is
            % not possible.
            %
            % Inputs
            % - obj: FarField object
            % - coorTypeString: 'spherical'|'Ludwig1'|'Ludwig2AE'|'Ludwig2EA'|'Ludwig3'|'power'
            % - setStdGrid:     Logical that determines if the grid should
            %                   be kept (false) or changed to the stadard
            %                   assiciated with the coorType (true
            %                   - default)
            %
            % Outputs
            % - obj: FarField object - possibly with added base grid
            %
            % Dependencies
            % -
            %
            % Created: 2019, Ridalise Louw
            % Updated: 2019-08-13, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   F = F.changeCoor('Ludwig3');
            %   F.plot('plotType','2D','showGrid',1,'output','E1')
            
            if nargin < 3, setStdGrid = true; end
            mustBeMember(coorTypeString, {'spherical','Ludwig1','Ludwig2AE','Ludwig2EA','Ludwig3','power'});
            handleCoorType = str2func(['coor2',coorTypeString]);
            obj = handleCoorType(obj,setStdGrid); 
        end
        
        function obj = coor2spherical(obj,setStdGrid)
            % COOR2SPHERICAL Change the current coorType to spherical
            %
            % See help changeCoor for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the coordinate system type of a power only pattern')
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'spherical')
                obj = obj.reset2Base;   % Make sure the grid is in the base format, and the E-fields if they can
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1,obj.E2] = getEspherical(obj);
                obj.coorType = 'spherical';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig1(obj,setStdGrid)
            % COOR2LUDWIG1 Change the current coorType to Ludwig1
            %
            % See help changeCoor for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the coordinate system type of a power only pattern')
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig1')
                obj = obj.reset2Base;   % Make sure the grid is in the base format, and the E-fields if they can
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1,obj.E2] = getELudwig1(obj);
                obj.coorType = 'Ludwig1';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2Ludwig2AE(obj,setStdGrid)
            % COOR2LUDWIG2AE Change the current coorType to Ludwig2 Az/El
            %
            % See help changeCoor for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the coordinate system type of a power only pattern')
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig2AE')
                obj = obj.reset2Base;   % Make sure the grid is in the base format, and the E-fields if they can
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1,obj.E2] = getELudwig2AE(obj);
                obj.coorType = 'Ludwig2AE';
            end
            if setStdGrid
                obj = obj.grid2AzEl;
            end
        end
        
        function obj = coor2Ludwig2EA(obj,setStdGrid)
            % COOR2LUDWIG2EA Change the current coorType to Ludwig2 El/Az
            %
            % See help changeCoor for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the coordinate system type of a power only pattern')
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig2EA')
                obj = obj.reset2Base;   % Make sure the grid is in the base format, and the E-fields if they can
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1,obj.E2] = getELudwig2EA(obj);
                obj.coorType = 'Ludwig2EA';
            end
            if setStdGrid
                obj = obj.grid2ElAz;
            end
        end
        
        function obj = coor2Ludwig3(obj,setStdGrid)
            % COOR2LUDWIG3 Change the current coorType to Ludwig3
            %
            % See help changeCoor for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the coordinate system type of a power only pattern')
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'Ludwig3')
                obj = obj.reset2Base;   % Make sure the grid is in the base format, and the E-fields if they can
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1,obj.E2] = getELudwig3(obj);
                obj.coorType = 'Ludwig3';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        function obj = coor2power(obj,setStdGrid)
            % COOR2POWER Change the current coorType to power
            %
            % See help changeCoor for details
            
            if nargin < 2, setStdGrid = true; end
            if ~strcmp(obj.coorType,'power')
                obj = obj.reset2Base;   % Make sure the grid is in the base format, and the E-fields if they can
                if ~isempty(obj.E1Base)
                    obj = obj.clearBaseFields;    % Clean the base - no going back from change to power
                end
                obj.E1 = sqrt(2*obj.eta0.*obj.getW);
                obj.E2 = [];
                obj.coorType = 'power';
                obj.polType = 'none';
            end
            if setStdGrid
                obj = obj.grid2PhTh;
            end
        end
        
        %% polarization type transformation methods
        function obj = changePol(obj,polTypeString)
            % CHANGEPOL Change the FarField object polarization type.
            %
            % obj = changePol(obj,polTypeString) transforms 
            % the current polarization type to that specified in 
            % polTypeString. 
            % This function calls the appropriate pol2* function to
            % excecute. Note that all the pol2* functions can be directly
            % called with the object as argument to do the same job.
            % 
            % In all cases, if no base field exists in the object
            % the base will be set before the transformation to maintain
            % the original provided field values for later use. If a base
            % field is present, the object will first be set into the base
            % field form before transformation - transformations thus always
            % happens using the original provided data in the base.
            %
            % If the object is already in the specified polarization type, nothing is
            % done.
            %
            % When a slant type is requested, the public property 'slant' is
            % used to determine the polarization angle.  A default of pi/4
            % radians is set. The parameter should be updated before
            % calling the pol2slant method.
            %
            % Inputs
            % - obj: FarField object
            % - polTypeString: 'linear'|'circular'|'slant'
            %
            % Outputs
            % - obj: FarField object - possibly with added base grid
            %
            % Dependencies
            % -
            %
            % Created: 2019, Ridalise Louw
            % Updated: 2019-08-16, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField;
            %   F = F.changePol('circular');
            %   F.plot('plotType','2D','showGrid',1,'output','E1')
            
            mustBeMember(polTypeString, {'linear','circular','slant'});
            handlePolType = str2func(['pol2',polTypeString]);
            obj = handlePolType(obj);   
        end
        
        function obj = pol2linear(obj)
            % POL2LINEAR Change the polarization to linear
            %
            % See help changePol for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the polarisation system type of a power only pattern')
            if ~strcmp(obj.polType,'linear')
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1, obj.E2] = getElin(obj);
                obj.polType = 'linear';
            end
        end
        
        function obj = pol2circular(obj)
            % POL2CIRCULAR Change the polarization to circular
            %
            % See help changePol for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the polarisation system type of a power only pattern')
            if ~strcmp(obj.polType,'circular')
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
                [obj.E1,obj.E2] = getEcircular(obj);
                obj.polType = 'circular';
            end
        end
        
        function obj = pol2slant(obj)
            % POL2SLANT Change the polarization to slant
            %
            % See help changePol for details
            
            assert(~strcmp(obj.coorType,'power'),'Cannot change the polarisation system type of a power only pattern')
            if ~strcmp(obj.polType,'slant')
                if isempty(obj.E1Base)
                    obj = obj.setBaseFields;    % Set E-fields base if none is present 
                end
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
            obj = handleCoorType(obj,false);
            if ~strcmp(objPolType,'none'), obj = handlePolType(obj); end
        end
        
        %% Base grid functions
        function obj = reset2Base(obj)
            % RESET2BASE Hard reset the FarField object to the base format
            
            if ~isempty(obj.xBase)
                obj.x = obj.xBase;
                obj.y = obj.yBase;
                obj.gridType = obj.gridTypeBase;
            end
            if ~isempty(obj.E1Base)
                obj.E1 = obj.E1Base;
                obj.E2 = obj.E2Base;
                obj.E3 = obj.E3Base;
                obj.coorType = obj.coorTypeBase;
                obj.polType = obj.polTypeBase;
            end
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
            if strcmp(coorTypeIn,'power')
                obj.polType = 'none';
            else
                obj = polTypeH(obj);
            end
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
            xR = obj1.xRange;
            yR = obj1.yRange;
            [stepX,stepY] = obj1.gridStep;
            if any(strcmp(obj1.gridType,{'DirCos','ArcSin'}))
                stepX = asin(stepX);
                stepY = asin(stepY);
                xR = asin(xR);
                yR = asin(yR);
            end
            [xmin,xmax] = deal(xR(1),xR(2));
            [ymin,ymax] = deal(yR(1),yR(2));

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
            
            if any(strcmp(obj1.gridType,{'DirCos','ArcSin'}))
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
            E1grid = zeros(Nxi*Nyi,obj1.Nf);
            if ~isempty(obj1.E2), E2grid = E1grid; end
            for ff = 1:obj1.Nf
                FFtemp = obj1.getFi(ff); % Much faster than letting interpolateGrid handle the frequency selection
                E1grid(:,ff) = interpolateGrid(FFtemp,'E1',xi,yi,1,hem);
                if ~isempty(obj1.E2), E2grid(:,ff) = interpolateGrid(FFtemp,'E2',xi,yi,1,hem); end
%                 E1grid(:,ff) = interpolateGrid(obj1,'E1',xi,yi,ff,hem);
%                 if ~isempty(obj1.E2), E2grid(:,ff) = interpolateGrid(obj1,'E2',xi,yi,ff,hem); end
            end
            % Remove the extra phase introduced by the interpolateGrid
            % function - this just keeps the real/imag and phase field
            % consistant with the plotting
            k = 2.*pi.*obj1.freqHz./obj1.c0;
            FFfact = exp(1i.*k.*obj1.r)./obj1.r;
            E1grid = bsxfun(@times,E1grid,FFfact);
            if ~isempty(obj1.E2)
                E2grid = bsxfun(@times,E2grid,FFfact);
            else
                E2grid = [];
            end
            
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
        function [] = plot(obj,varargin)
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
            %   -- norm:        Boolean (false) to normalise output
            %   -- dynamicRange_dB: A (positive) dB value for the magnitude plot dynamic range (40)
            %   -- scaleMag:    {('dB') | 'lin'}, only used for magnitude plots
            %   -- scalePhase:  {('deg') | 'rad'}, only used for phase plots
            %   -- freqUnit:    {('GHz') | 'Hz' | 'kHz' | 'MHz' | 'THz'}
            %   -- cutConstant: determines along which axis to cut in 1D plots {('x') | 'y'} 
            %   -- cutValue:    Any value in the available angle range for 1D cuts range in rad (0)
            %   -- step:        Plot step size in deg. Can be empty - then the available data will
            %                   be used and no surface will be plotted.  If
            %                   not, a griddata interpolant will be made.
            %                   Can be a 2 element vector, [xStep,yStep],
            %                   and if step == 0 an automatic step size
            %                   will be selected (0).
            %   -- plotProperties: can be a variety of name, value pairs
            %                      including LineWidth, LineStyle, Color (like '-.')
            %   -- showGrid:   Boolean (false) to show the 2D grid where the data is
            %                  calculated before interpolation.
            %   -- hemisphere: Used in gridTypes DirCos and ArcSin {('top')
            %                  | 'bot'}
            %
            % Outputs
            % - 
            %
            % Dependencies
            % - MATLAB Antennas Toolbox for 3D plot
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
            %   FF.plot('plotType','2D')
            
            
            narginchk(1,40);
            
            % Parsing through the inputs
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
            addParameter(parseobj,'step',0,typeValidationstep);     % In degrees
            
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
            
            % Sort out the plot grid and names
            
            % Get valid positions for the plot
            if any(strcmp(obj.gridType,{'DirCos','ArcSin'}))
                % Try to get the direction cosines from the base grid definition - if the
                % base definition is not a direction cosine type it can contain
                % information over the full sphere.
                if ~isempty(obj.xBase)
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
%             elseif strcmp(obj.gridType,'TrueView')
%                 valAng = sqrt((obj.x./pi).^2 + (obj.y./pi).^2) <= max(obj.th)./pi;
            else
                valAng = ones(obj.Nang,1);
            end
            
            % Get the original grid and output
            X = obj.x;
            Y = obj.y;
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
            
            % Interpolate or not? Use the supplied step to decide. Empty
            % will use what we have on the grid - may crash.  Zero will see
            % if the grid is uniform and then not interpolate, if not, it
            % will.
            interpOut = true;
            if isempty(step)
                if any(strcmp(plotType,{'cartesian','polar'}))
                    error('Cant do 1D plot with empty step size - this is reserved for making 2D grids etc.')
                else
                    NxPlot = obj.Nx;
                    NyPlot = obj.Ny;
                    assert(NxPlot*NyPlot == obj.Nang,'Cannot pass empty step when the grid is not uniform')
                    interpOut = false;
                end
            else
                if any(strcmp(obj.gridType,{'DirCos','ArcSin'}))
                    step = sind(step);
                else
                    step = deg2rad(step);
                end
                if isscalar(step) == 1
                    if step ~= 0
                        step = ones(1,2).*step;
                    else
                        NxPlot = obj.Nx;
                        NyPlot = obj.Ny;
                        if NxPlot*NyPlot == obj.Nang
                            interpOut = false;
                        else
                            [step(1),step(2)] = obj.gridStep;
                        end
                    end
                end
            end
                
            if interpOut
                [ximin,ximax] = deal(obj.xRange(1),obj.xRange(2));
                [yimin,yimax] = deal(obj.yRange(1),obj.yRange(2));
                NxPlot = round((ximax - ximin)/step(1)) + 1;
                NyPlot = round((yimax - yimin)/step(2)) + 1;
                xivect = linspace(ximin,ximax,NxPlot);
                yivect = linspace(yimin,yimax,NyPlot);
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
            else
                Xi = reshape(obj.x,NyPlot,NxPlot);
                Yi = reshape(obj.y,NyPlot,NxPlot);
                xi = obj.x;
                yi = obj.y;
                Zi = Z;
                Zi(~valAng) = NaN;
            end
            
            % Assign axis names
            switch obj.gridType
                case 'DirCos'
                    xiplot = xi;
                    yiplot = yi;
                    xnamePlot = [obj.xname, ' = sin(\theta)cos(\phi)'];
                    ynamePlot = [obj.yname, ' = sin(\theta)sin(\phi)'];
                    axisUnit = '';
                otherwise
                    X = rad2deg(X);
                    Y = rad2deg(Y);
                    Xi = rad2deg(Xi);
                    Yi = rad2deg(Yi);
                    xiplot = rad2deg(xi);
                    yiplot = rad2deg(yi);
                    axisUnit = '(deg)';
                    xnamePlot = [obj.xname, ' ' ,axisUnit];
                    ynamePlot = [obj.yname, ' ' ,axisUnit];
            end
            
            % Condition outputs
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
                        unit = 'V/m';
                        compName = strrep(obj.([output,'name']),'_','\');
                    otherwise
                        error(['output: ' output,' not implemented in plot function'])
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
            
            % Make the plots
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
                    xlabel(xnamePlot)
                    ylabel(ynamePlot)
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
                            xlab = ynamePlot;
                            cutName = obj.xname;
                        case 'y'
                            plotHandle(xiplot.*xscale,Ziplot,'LineStyle',LineStyle,'LineWidth',LineWidth,'Color',Color), grid on
                            xlab = xnamePlot;
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
        
        function [ax] = plotJones(FF1,FF2,varargin)
            % PLOTJONES Plots the Jones matrix type representation FarFields
            
            % function [ax] = plotJones(FF1,FF2,varargin)
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
            %
            % Returns a vector of axis handles
            
            
            % Parse input
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
            addParameter(parseobj,'step',0,typeValidationstep);     % In degrees
            
            parse(parseobj, FF1, FF2, varargin{:});
            
            freqIndex = parseobj.Results.freqIndex;
            dynamicRange_dB = parseobj.Results.dynamicRange_dB;
            step = parseobj.Results.step;
            
            % Plot the result
            [w,h] = deal(0.4);
            botTop = 0.56;
            botBot = 0.08;
            leftLeft = 0.1;
            leftRight = 0.52;
            figPos = [300,110,[4,3].*220];
%             ax(1:4) = axis;
            if ~isGridEqual(FF1,FF2)
                error('Base grids should be identical for the two input fields');
            else
%                 figure
                ax(1) = subplot('position',[leftLeft botTop w h]);
%                 ax(1) = subplot(2,2,1,'align');
                plot(FF1,'output','E1','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex);
                title('J_{11}')
                xlabel('')
                ax(1).XLabel.Visible = 'off';
%                 ax(1).XTickLabel = [];

                ax(2) = subplot('position',[leftRight botTop w h]);
%                 ax(2) = subplot(2,2,2,'align');
                plot(FF1,'output','E2','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex);
                title('J_{12}')
                ax(2).XLabel.Visible = 'off';
%                 ax(2).XTickLabel = [];
                ax(2).YLabel.Visible = 'off';
%                 ax(2).YTickLabel = [];

                ax(3) = subplot('position',[leftLeft botBot w h]);
%                 ax(3) = subplot(2,2,3,'align');
                plot(FF2,'output','E1','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex);
                title('J_{21}')
                
                ax(4) = subplot('position',[leftRight botBot w h]);
%                 ax(4) = subplot(2,2,4,'align');
                plot(FF2,'output','E2','outputType','mag','plotType','2D','scaleMag','dB','norm',0,'step',step,'dynamicRange_dB',dynamicRange_dB,'freqIndex',freqIndex);
                title('J_{22}')
                ax(4).YLabel.Visible = 'off';
%                 ax(4).YTickLabel = [];
                
                for aa = 1:4
                    ax(aa).ActivePositionProperty = 'position';
                    ax(aa).XMinorTick = 'on';
                    ax(aa).YMinorTick = 'on';
                    ax(aa).Box = 'on';
                    set(ax(aa),'fontsize',11);
                end
                set(gcf,'Position',figPos)
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
            
            % Parse input
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
            
            
            % Plot the result
            
            % Estimate a nice step size
            thRange = (max(FF.y) - min(FF.y));
            Nstep = round(sqrt(FF.Nang));
            Nstep = Nstep + 1-mod(Nstep,2);
            step = thRange/(Nstep-1);
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
                    FF = FF.setRangeSph('sym','360');
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
            
            gridTypeIn = obj.gridType;
            
            % Evaluate the field on the base grid - this is where the output function
            % should be best suited for interpolation.  Can't do this for
            % astroGrids, or cases where we have transformed from local to
            % astro or vice versa, or when the base is not defined yet
            if ~obj.baseTypeDifferent && ~isempty(obj.xBase) 
                obj = obj.grid2Base; 
            end
            
            % Shift to x = sym and y = 180 range (if applicable) - this is where the DirCos spits
            % everything out after transforming, and it is standard for
            % extending the range for smooth interpolants
            if any(strcmp(obj.gridType,obj.localGrids))
                obj = obj.setRangeSph('sym','180');
            end
            % Get xi and yi in the base gridType, and on the [-180,180] x-domain for the
            % angular grids
            grid2DirCoshandle = str2func([gridTypeIn,'2DirCos']);
            [ui,vi,wi] = grid2DirCoshandle(xi,yi);
            valAngi = true(size(ui));
            if any(strcmp(gridTypeIn,{'DirCos','ArcSin'}))
                % Find the invalid points included by the external meshgrid
                valAngi = sqrt(ui.^2 + vi.^2) <= max(sin(obj.th));
                % Check for bottom hemisphere plot - fix the w to be the negative root
                if strcmp(hemisphere,'bot')
                    wi = -wi;
                end
            elseif strcmp(gridTypeIn,'TrueView')
                % Sort out the TrueView special case invalid points
                valAngi = sqrt((xi./pi).^2 + (yi./pi).^2) <= max(obj.th)./pi;
            end
            DirCos2baseHandle = str2func(['DirCos2',obj.gridType]);
            [xi_bGT,yi_bGT] = DirCos2baseHandle(ui,vi,wi);
            
            % Get the valid angle positions - already in baseGrid here, but shifted to
            % +- 180 degrees
            % Get the indexes only, no later reshaping done, different from the grids
            % required for plotting in FarField.plot
            if strcmp(gridTypeIn,'DirCos') || strcmp(gridTypeIn,'ArcSin')
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
            if any(strcmp(obj.gridType,obj.sphereGrids))
%             if strcmp(obj.gridType,'PhTh') || strcmp(obj.gridType,'AzEl') || strcmp(obj.gridType,'ElAz')
                tol = deg2rad(2); % Check for points close to the edges in azimuth
                if abs(min(xVal) + pi) < tol && abs(max(xVal) - pi) < tol
                    edgeAngDeg = 180 - edgeAngExtent_deg;
                    iNeg = find(xVal > deg2rad(edgeAngDeg));
                    iPos = find(xVal < deg2rad(-edgeAngDeg));
                    xVal = [xVal(iNeg)-2*pi;xVal;xVal(iPos)+2*pi];
                    yVal = [yVal(iNeg);yVal;yVal(iPos)];
                    zVal = [zVal(iNeg);zVal;zVal(iPos)];
                end
%                     % Also extend the y-axis
%                     if abs(min(yVal) - 0) < tol
%                         edgeAngDeg = edgeAngExtent_deg;
%                         iNeg = find(xVal < 0 & yVal < deg2rad(edgeAngDeg));
%                         iPos = find(xVal > 0 & yVal < deg2rad(edgeAngDeg));
%                         xVal = [xVal;xVal(iNeg)+pi;xVal(iPos)-pi];
%                         yVal = [yVal;-yVal(iNeg);-yVal(iPos)];
%                         zVal = [zVal;zVal(iNeg);zVal(iPos)];
%                     end
%                     if abs(max(yVal) - pi) < tol
%                         edgeAngDeg = 180 - edgeAngExtent_deg;
%                         iNeg = find(xVal < 0 & yVal > deg2rad(edgeAngDeg));
%                         iPos = find(xVal > 0 & yVal > deg2rad(edgeAngDeg));
%                         xVal = [xVal;xVal(iNeg)+pi;xVal(iPos)-pi];
%                         yVal = [yVal;2*pi-yVal(iNeg);2*pi-yVal(iPos)];
%                         zVal = [zVal;zVal(iNeg);zVal(iPos)];
%                     end
                
            end
            
            % Remove duplicate differing values completely from the set - interpolate
            % over them.  This happens at poles for certain coordinate projections
            % (like at th = 180 in Ludwig3 for instance, of th=0 and 180 for spherical)
            % First find duplicate domain values
            [~,iUnique] = unique([xVal,yVal],'rows');
            removePoints = [];
            if length(iUnique) < length(xVal) && 0
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
            
            % Fix the poles on a case-by-case basis - mostly a futile
            % excercise so far...
            % Reinterpolate the points in the pole on the new grid with the
            % interpolated function which is already spread out around the
            % pole
            tol = deg2rad(2);
            if strcmp(obj.coorTypeBase,'Ludwig2AE') && strcmp(obj.coorType,'spherical')
                % Pole at th = 0:
                % Check if there are multiple points there
                if length(find(yi == 0)) > 1
                    iPole = find(abs(yi - 0) < tol);
                    iClose = find(abs(yi) < deg2rad(edgeAngExtent_deg)); 
                    iBase = setxor(iPole,iClose); % Remove the pole region from close in points for interpolation
                    iPole0 = find(yi == 0 & xi == 0);
%                     % Mirror across pole, and include the actual pole
%                     xiBase = [xi(iBase);wrap2pi(xi(iBase) + pi);xi(iPole0)]; 
%                     yiBase = [yi(iBase);-yi(iBase);yi(iPole0)];
%                     ZiBase = [Zi(iBase);Zi(iBase);Zi(iPole0)];
                    xiBase = [xi(iBase);xi(iPole0)];
                    yiBase = [yi(iBase);yi(iPole0)];
                    ZiBase = [Zi(iBase);Zi(iPole0)];
                    Zf2 = scatteredInterpolant(xiBase,yiBase,ZiBase,'linear');
                    Zi(iPole) = Zf2(xi(iPole),yi(iPole));
                end
                
            end
            
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
            if strcmp(FF.symmetryBOR,'none')
                FF = FF.getBORpattern;
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
            %   rotxsph, rotysph, rotzsph, rotGRASPsph, rotEulersph
            % rotAng is the associated angle in rad. Scalar for rotations
            % around an axis, and [3x1] for GRASP or Euler rotations
            %
            % onlyRotPowerPattern is an optional argument which speeds up
            % the method in the case where only the rotated power pattern
            % is of interest.  The field values will be arbitrary, but the
            % power pattern (directivity etc.) will be correct.  Used often
            % for noise temeperature calculations.
            
            if nargin == 5
                warning('onlyRotPowerPattern deprecated - use a coorType = power pattern instead, and pass 4 arguments only')
            end
            if nargin < 4
                onlyRotPowerPattern = false;
            end
            
            % This will probably depend on the pattern type which one is
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
%             FFsph = FFsph.setXrange('sym'); % Always work in symmetrical xRange
            FFsph = FFsph.setRangeSph('sym');
            if ~strcmp(obj1.coorType,'power'), FFsph = coorHandle(FFsph,false); end % No coordinate transform required for power only
            
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
            % Get a sensible grid size
            Nxin = FFsph.Nx;
            Nyin = FFsph.Ny;
            if Nxin*Nyin > FFsph.Nang
                Nxin = sqrt(FFsph.Nang/2);
                Nyin = Nxn/2;
                Nxin = floor(Nxin/2)*2 + 1;
                Nyin = floor(Nyin/2)*2 + 1;
            end
            stepx = (max(FFsph.ph) - min(FFsph.ph))./(Nxin-1);
            stepy = (max(FFsph.th) - min(FFsph.th))./(Nyin-1);
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
            FFsph = FFsph.roundGrid;  % Found some rare, strange behaviour in the interpolant if this is not done
            
            % Coordinate systems
            C0 = CoordinateSystem;
            rotHandleStr = func2str(rotHandle);
            rotHandleCoor = str2func(rotHandleStr(1:end-3));
            Crot = rotHandleCoor(C0,rotAng);
            
            if strcmp(obj1.coorType,'power')
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
                FFsph = FFsph.setBase;
                for ff = 1:obj1.Nf
                    Ugrid(:,ff) = interpolateGrid(FFsph,'U',xi,yi,ff,'top');
                end
                FFsph = FarField.farFieldFromPowerPattern(xi,yi,Ugrid,FFsph.freq,'fieldPol','power','freqUnit',FFsph.freqUnit,'symmetryXZ',obj1.symmetryXZ,'symmetryYZ',obj1.symmetryYZ,'symmetryXY',obj1.symmetryXY,'symmetryBOR',obj1.symmetryBOR);
            else
                % Perform the rotation of the field vectors
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
                    FFsph = FarField.farFieldFromPowerPattern(xi,yi,Ugrid,FFsph.freq,'fieldPol','linearY','freqUnit',FFsph.freqUnit,'symmetryXZ',obj1.symmetryXZ,'symmetryYZ',obj1.symmetryYZ,'symmetryXY',obj1.symmetryXY,'symmetryBOR',obj1.symmetryBOR);
                else
                    FFsph = FFsph.currentForm2Base(stepDeg,rad2deg([xmin,xmax;ymin,ymax]));
                end
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
%             obj = obj.setXrange(obj1.xRangeType);
            obj = obj.setRangeSph(obj1.xRangeType);
            obj = obj.currentForm2Base();
        end
        
        function obj = shift(obj,shiftVect)
            % SHIFT Shifts the FarField by a specified distance.
            
            % Shifts the FarField by a distance specified in the Pnt3D input
            % shiftVect (only uses the first entry)
            % shiftVect can also be a vector of length 3 with elements
            % [delX,delY,delZ] in m
            
            assert(~strcmp(obj.coorType,'power'),'Shift meaningless on power only patterns')
            assert(~any(strcmp(obj.gridType,obj.astroGrids)),'Shift only operates on local grids/projections.')
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
            
            [obj1base,obj2base] = mathSetup(obj1,obj2);
            
            obj = obj1base;
            obj.E1 = obj1base.E1 + obj2base.E1;
            if ~isempty(obj1.E2) && ~isempty(obj1.E2)
                obj.E2 = obj1base.E2 + obj2base.E2;
            else
                obj.E2 = [];
            end
            if ~isempty(obj1.E3) && ~isempty(obj1.E3)
                obj.E3 = obj1base.E3 + obj2base.E3;
            else
                obj.E3 = [];
            end
            try
                obj.Prad = obj.pradInt;
            catch ME
                obj.Prad = ones(size(obj.freqHz)).*4*pi/(2.*obj.eta0);
            end
            obj.radEff = ones(size(obj.freqHz));
            obj = setBase(obj);
        end
        
        function obj = minus(obj1,obj2)
            % MINUS Subtract two FarFields
            
            [obj1base,obj2base] = mathSetup(obj1,obj2);
            
            obj = obj1base;
            obj.E1 = obj1base.E1 - obj2base.E1;
            if ~isempty(obj1.E2) && ~isempty(obj1.E2)
                obj.E2 = obj1base.E2 - obj2base.E2;
            else
                obj.E2 = [];
            end
            if ~isempty(obj1.E3) && ~isempty(obj1.E3)
                obj.E3 = obj1base.E3 - obj2base.E3;
            else
                obj.E3 = [];
            end
            try
                obj.Prad = obj.pradInt;
            catch ME
                obj.Prad = ones(size(obj.freqHz)).*4*pi/(2.*obj.eta0);
            end
            obj.radEff = ones(size(obj.freqHz));
            obj = setBase(obj);
        end
        
        function obj = times(obj1,obj2)
            % TIMES Multiply two FarFields
            
            [obj1base,obj2base] = mathSetup(obj1,obj2);
            
            obj = obj1base;
            obj.E1 = obj1base.E1 .* obj2base.E1;
            if ~isempty(obj1.E2) && ~isempty(obj1.E2)
                obj.E2 = obj1base.E2 .* obj2base.E2;
            else
                obj.E2 = [];
            end
            if ~isempty(obj1.E3) && ~isempty(obj1.E3)
                obj.E3 = obj1base.E3 .* obj2base.E3;
            else
                obj.E3 = [];
            end
            obj.Prad = obj.pradInt;
            obj.radEff = ones(size(obj.freqHz));
            obj = setBase(obj);
        end
        
        function obj = conj(obj1)
            % CONJ Get the complex conjugate of E-field values in a FarField.
            
            obj = obj1;
            obj.E1 = conj(obj1.E1);
            obj.E2 = conj(obj1.E2);
            obj = setBase(obj);
        end
        
        function obj = abs(obj1)
            % ABS Get the absolute value of E-field components in a FarField.
            
            obj = obj1;
            obj.E1 = abs(obj1.E1);
            obj.E2 = abs(obj1.E2);
            obj = setBase(obj);
        end
        
        function obj = scale(obj1,scaleFactor)
            % SCALE Scale E-field components by a scaleFactor in a FarField.
            
            % Scale the FarField object E-fields by the scaleFactor
            obj = obj1;
            obj.E1 = obj1.E1.*scaleFactor;
            if ~isempty(obj1.E2)
                obj.E2 = obj1.E2.*scaleFactor;
            else
                obj.E2 = obj1.E2;
            end
            obj.Prad = obj1.Prad.*(abs(scaleFactor).^2);
            obj = setBase(obj);
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
            % Set possible NaN's to 0
            iNaN1 = isnan(obj.E1);
            obj.E1(iNaN1) = 0;
            rmsE1 = rms(obj.E1,DIM);
            
            iNaN2 = isnan(obj.E2);
            obj.E2(iNaN2) = 0;
            rmsE2 = rms(obj.E2,DIM);
            
            iNaN3 = isnan(obj.E3);
            obj.E3(iNaN3) = 0;
            rmsE3 = rms(obj.E3,DIM);
        end
        
        function T = convPower(obj1,obj2)
            % CONVPOWER Convolve the power patterns, over the full sphere,
            % of two FarField objects.
            
            % Convolve the power patterns, over the full sphere, of two
            % FarField objects. Typically used for antenna temperature
            % calculations
            [obj1,obj2] = mathSetup(obj1,obj2);
            
            P = obj1.getU.*obj2.getU;
            E1U = sqrt(P./obj1.r^2.*2.*obj1.eta0);
%             PradDummy = 1; % Dummy so the constructor does not integrate
%             FF_T = FarField.farFieldFromPowerPattern(obj1.x,obj1.y,P,obj1.freqHz,PradDummy,'gridType',obj1.gridType,'fieldPol','power',...
%                 'symmetryXZ',obj1.symmetryXZ,'symmetryYZ',obj1.symmetryYZ,'symmetryXY',obj1.symmetryXY,'symmetryBOR',obj1.symmetryBOR);
%             FF_T = FarField.farFieldFromPowerPattern(obj1.x,obj1.y,P,obj1.freqHz,'gridType',obj1.gridType,'fieldPol','power');
            % Make dummy Prad and radEff to formally calculate with
            % intergal - this way the constructor does not integrate
            FF_T = FarField(obj1.x,obj1.y,E1U,[],obj1.freq,1,1,...
                'coorType','power','polType',obj1.polType,'gridType',obj1.gridType,'freqUnit',obj1.freqUnit,'r',obj1.r,...
                'symmetryXZ',obj1.symmetryXZ,'symmetryYZ',obj1.symmetryYZ,'symmetryXY',obj1.symmetryXY,'symmetryBOR',obj1.symmetryBOR);
            
            T = FF_T.pradInt; 
%             T = FF_T.Prad;
        end
        
        %% Frequency and field modifications
        function obj = getFi(obj1,freqIndex)
            % GETFI Returns an object only containing the results in
            % freqIndex.
            
            % Returns an object only containing the results in freqIndex
            obj = obj1;
            obj.E1 = obj1.E1(:,freqIndex);
            if ~isempty(obj1.E2), obj.E2 = obj1.E2(:,freqIndex); end
            if ~isempty(obj1.E3), obj.E3 = obj1.E3(:,freqIndex); end
            obj.freq = obj1.freq(freqIndex);
            obj.Prad = obj1.Prad(freqIndex);
            obj.radEff = obj1.radEff(freqIndex);
            if ~isempty(obj1.xBase)
                obj = obj.setBaseGrid;
            end
            if ~isempty(obj1.E1Base)
                obj = obj.setBaseFields;
            end
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
                [stepX,stepY] = obj1.gridStep;
                if any(strcmp(obj1.gridTypeBase,{'DirCos','ArcSin'}))
%                     stepX = asin(min(abs(diff(unique(obj1.xBase)))));
%                     stepY = asin(min(abs(diff(unique(obj1.yBase)))));
                    stepX = asin(stepX);
                    stepY = asin(stepY);
%                 else
%                     % Sort out rounding errors for degrees
%                     stepX = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.xBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
%                     stepY = deg2rad(round(rad2deg(min(abs(diff(unique(obj1.yBase)))))*10^obj1.nSigDig)/10^obj1.nSigDig);
                end
                
                gridHandle = str2func(['grid2',gridTypeIn]);
                coorHandle = str2func(['coor2',coorTypeIn]);
                obj1 = obj1.grid2TrueView;
                if ~strcmp(obj1.coorType,'power')
                    obj1 = obj1.coor2Ludwig3(false);   % Always work in H/V for symmetry calculations...
                end
                
                % Initialise
                XIn = [obj1.x];
                YIn = [obj1.y];
                E1In = [obj1.E1];
                E2In = [obj1.E2];
                if obj1.symXZ
                    XIn = [XIn;XIn];
                    YIn = [YIn;-YIn];
                    E1In = [E1In;obj1.symXZ.*E1In]; % Mirror according to symmetry
                    if ~isempty(E2In), E2In = [E2In;-obj1.symXZ.*E2In]; end  % Mirror according to symmetry
                end
                if obj1.symYZ
                    XIn = [XIn;-XIn];
                    YIn = [YIn;YIn];
                    E1In = [E1In;-obj1.symYZ.*E1In];  % Mirror according to symmetry
                    if ~isempty(E2In), E2In = [E2In;obj1.symYZ.*E2In]; end % Mirror according to symmetry
                end
                % Object for grid/coor transformation
                objD = FarField(XIn,YIn,E1In,E2In,obj1.freq,obj1.Prad.*2,obj1.radEff,...
                'coorType',obj1.coorType,'polType',obj1.polType,'gridType',obj1.gridType,'freqUnit',obj1.freqUnit,...
                'E3',obj1.E3,'r',obj1.r,'orientation',obj1.orientation,'earthLocation',obj1.earthLocation,'time',obj1.time);
            
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
            phRange = max(obj.ph) - min(obj.ph);
            thRange = max(obj.th) - min(obj.th);
            tol = 10^(-obj.nSigDig);
            y = ((abs(round(rad2deg(phRange)) - (360/2^(sum(abs([obj.symXZ,obj.symYZ]))))) < tol) & (abs(round(rad2deg(thRange)) - 180/2^abs(obj.symXY)) < tol)) |...
                ((abs(round(rad2deg(phRange)) - 180) < tol) & (abs(round(rad2deg(thRange)) - 360) < tol));
        end
        
        function y = isGridUniform(obj)
            % ISGRIDUNIFORM Checks for a plaid, monotonic, uniform grid format. 
            
            % Test for a plaid, monotonic, uniform grid format
            
            % Can be slow for huge objects
            ny = obj.Ny;
            nx = obj.Nx;
            
            if nx*ny == obj.Nang
                X = reshape(obj.x,ny,nx);
                Y = reshape(obj.y,ny,nx);
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
            if any(strcmp(currGridType,obj.astroGrids))
                obj = obj.grid2PhTh; % Do in the old orientation
                obj.orientation = newOrientation;
                eval(['obj = grid2',currGridType,'(obj);']);
            else
                obj.orientation = newOrientation; % Just change the parameter
            end
        end
        
        function obj = setTime(obj,newTime)
            % SETTIME Set time.
            
            currGridType = obj.gridType;
            if any(strcmp(currGridType,obj.astroGrids))
                obj = obj.grid2Horiz; % Do in the old time
                obj.time = newTime;
                eval(['obj = grid2',currGridType,'(obj);']);
            else
                obj.time = newTime; % Just change the parameter
            end
        end
        
        function obj = setEarthLocation(obj,newEarthLocation)
            % SETEARTHLOCATION Set the antenna location on Earth.
            
            currGridType = obj.gridType;
            if any(strcmp(currGridType,obj.astroGrids))
                obj = obj.grid2Horiz; % Do at the old location
                obj.earthLocation = newEarthLocation;
                eval(['obj = grid2',currGridType,'(obj);']);
            else
                obj.earthLocation = newEarthLocation; % Just change the parameter
            end
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
            assert(strcmp(obj.gridType,'PhTh'),'Grid type must be PhTh for GRASP cut write')
            
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
            % READGRASPGRD Create a FarField object from a GRASP .grd file. 
            %
            % FF = readGRASPgrd(pathName,varargin) loads a FarField object
            % from the GRASP .grd file at pathName. Can have several optional
            % arguments describing the local field as name value pairs.
            % Not all GRASP functionality supported yet...
            % 
            % Inputs
            % - pathName: Full path and filename string. Can be empty -
            %               then gui will request an grd file
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   acceptable keywords are  
            %   -- symmetryXZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryYZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryXY:  {('none')|'electric'|'magnetic'}
            %   -- symBOR:      {('none')|'BOR0'|'BOR1'}
            %   -- r:           See FarField constructor help for details
            %   -- orientation: See FarField constructor help for details
            %   -- earthLocation: See FarField constructor help for details
            %   -- time:        See FarField constructor help for details
            %
            % Outputs
            % - FF:    Farfield object
            %
            % Dependencies
            % -
            %
            % Created: 2019-06-10, Dirk de Villiers
            % Updated: 2019-08-09, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField.readGRASPgrd;
            %   F.plot
            
            
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
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            if nargin == 0
                [name,path] = uigetfile('*.grd');
                pathName = [path,name];
            end
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
                    [DX,DY] = deal(0);
                    if NX > 1
                        DX = (XE-XS)/(NX - 1);
                    end
                    if NY > 1
                        DY = (YE-YS)/(NY - 1);
                    end
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
%             [Xmat,Ymat] = ndgrid(X,Y);
            [Xmat,Ymat] = meshgrid(X,Y);
            x = Xmat(:);
            y = Ymat(:);
            swopE1E2 = false;
            switch ICOMP
                case 1
                    polType = 'linear';
                    coorType = 'spherical';
                case 2
                    swopE1E2 = true;
                    polType = 'circular';
                    coorType = 'spherical';
                case 3
                    polType = 'linear';
                    coorType = 'Ludwig3';
                otherwise
                    error(['ICOMP ',num2str(ICOMP),' case not implemented yet'])
            end
            swopXY = false;
            switch IGRID
                case 1
                    gridType = 'DirCos';
                    [xScale,yScale] = deal(1);
                case 9
                    gridType = 'AzEl';
                    [xScale,yScale] = deal(pi/180);
                case 5
                    gridType = 'TrueView';
                    xScale = -pi/180;
                    yScale = -xScale;
                case 10
                    swopXY = true;
                    gridType = 'ElAz';
                    [xScale,yScale] = deal(pi/180);
                case 7
                    gridType = 'PhTh';
                    [xScale,yScale] = deal(pi/180);
                otherwise
                    error(['IGRID ',num2str(IGRID),' case not implemented yet'])
            end
            
            if swopXY
                % Swop definitions to be consistent with GRASP
                xt = x;
                x = y;
                y = xt;
            end
            
            % Fix the shapes of the field data
            for ff = 1:NSET
                e1 = reshape(E1(:,ff),NX,NY).';
                E1(:,ff) = e1(:);
                e2 = reshape(E2(:,ff),NX,NY).';
                E2(:,ff) = e2(:);
                if NCOMP == 3
                    e3 = reshape(E3(:,ff),NX,NY).';
                    E3(:,ff) = e3(:);
                end
            end
            if swopE1E2
                % Swop definitions to be consistent with GRASP
                Et = E1;
                E1 = E2;
                E2 = Et;
            end
                
            % keyboard;
            eta0 = 3.767303134749689e+02;
            Prad = ones(size(freq)).*4*pi./(2*eta0);
            radEff = ones(size(freq));
            FF = FarField(x.*xScale,y.*yScale,E1,E2,freq,Prad,radEff,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'symmetryXZ',symmetryXZ,'symmetryYZ',symmetryYZ,'symmetryXY',symmetryXY,'symmetryBOR',symmetryBOR,'E3',E3,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
            if swopXY
                FF = FF.sortGrid;
            end
        end

        function FF = readGRASPcut(pathName,nr_freq,nr_cuts,varargin)
            % READGRASPCUT Create a FarField object from a GRASP .cut file.
            %
            % FF = readGRASPcut(pathName,nr_freq,nr_cuts,varargin) loads a FarField object
            % from the GRASP .cut file at pathName. Can have several optional
            % arguments describing the local field as name value pairs.
            % Not all GRASP functionality supported yet - can only handle
            % polar ph-th grids (most often used anyway)
            % 
            % Inputs
            % - pathName: Full path and filename string. Can be empty -
            %               then gui will request a cut file
            % - nr_freq:  Number of frequencies in the file
            % - nr_cuts:  Number of ph cuts in the file
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   acceptable keywords are  
            %   -- freq:        See FarField constructor help for details
            %   -- freqUnit:    {('Hz')|'kHz'|'MHz'|'GHz'|'THz'}
            %   -- symmetryXZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryYZ:  {('none')|'electric'|'magnetic'}
            %   -- symmetryXY:  {('none')|'electric'|'magnetic'}
            %   -- symBOR:      {('none')|'BOR0'|'BOR1'}
            %   -- r:           See FarField constructor help for details
            %   -- orientation: See FarField constructor help for details
            %   -- earthLocation: See FarField constructor help for details
            %   -- time:        See FarField constructor help for details
            %
            % Outputs
            % - FF:    Farfield object
            %
            % Dependencies
            % -
            %
            % Created: 2019-04-22, Dirk de Villiers
            % Updated: 2019-08-09, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField.readGRASPcut;
            %   F.plot
            
            
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
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            if nargin == 0
                [name,path] = uigetfile('*.cut');
                pathName = [path,name];
                nr_freq = input('Enter number of frequencies:');
                nr_cuts = input('Enter number of cuts:');
            end
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
            fid = fopen(pathName,'rt');
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
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            if nargin == 0
                [name,path] = uigetfile('*.ffe');
                pathName = [path,name];
            end
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
            % READCSTFFS Create a FarField object from a CST .ffs file.
            %
            % FF = readCSTffs(pathName,varargin) loads a FarField object
            % from the CST .ffs file at pathName. Can have several optional
            % arguments describing the local field as name value pairs
            % 
            % Inputs
            % - pathName: Full path and filename string. Can be empty -
            %               then gui will request an ffs file
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   acceptable keywords are  
            %   -- r:           See FarField constructor help for details
            %   -- orientation: See FarField constructor help for details
            %   -- earthLocation: See FarField constructor help for details
            %   -- time:        See FarField constructor help for details
            %
            % Outputs
            % - FF:    Farfield object
            %
            % Dependencies
            % -
            %
            % Created: 2019-06-10, Dirk de Villiers
            % Updated: 2019-06-10, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField.readCSTffs;
            %   F.plot
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readCSTffs';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            if nargin == 0
                [name,path] = uigetfile('*.ffs');
                pathName = [path,name];
            end
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

        function FF = readCSTtxt(pathName,varargin)
            % READCSTTXT Create a FarField object from a CST .txt file
            %
            % FF = readCSTtxt(pathName,freq,Prad,radEff,varargin) loads a FarField object
            % from the CST .txt file at pathName. Can have several optional
            % arguments describing the local field as name value pairs. The
            % file should be generated by exporting the ASCII plot data
            % from a 2D FarField plot in CST. This format only allows a
            % single frequency at a time, and has no information about the
            % power, efficiency or frequency, so these can be optionally
            % provided
            % 
            % Inputs
            % - pathName: Full path and filename string. Can be empty -
            %               then gui will request an ffs file
            % - freq: Frequencies where the fields are defined in Hz, [1 x Nf]
            % - Prad: Radiated power at each frequency in W, [1 x Nf]
            % - radEff: Radiation efficiency at each frequency, [1 x Nf]
            % * Arbitrary number of pairs of arguments: ...,keyword,value,... where
            %   acceptable keywords are
            %   -- freq:        See FarField constructor help for details
            %   -- freqUnit:    {('Hz')|'kHz'|'MHz'|'GHz'|'THz'}
            %   -- r:           See FarField constructor help for details
            %   -- orientation: See FarField constructor help for details
            %   -- earthLocation: See FarField constructor help for details
            %   -- time:        See FarField constructor help for details
            %
            % Outputs
            % - FF:    Farfield object
            %
            % Dependencies
            % -
            %
            % Created: 2019-08-13, Dirk de Villiers
            % Updated: 2019-08-13, Dirk de Villiers
            %
            % Tested : Matlab R2018b
            %  Level : 2
            %   File : testScript_FarField.m
            %
            % Example
            %   F = FarField.readCSTtxt;
            %   F.plot('plotType','2D','showGrid',1)
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readCSTtxt';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            typeValidation_freq = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','nondecreasing','nrows',1},'FarField');
            parseobj.addOptional('freq',1e9,typeValidation_freq);
            
            typeValidation_power = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan'},'FarField');
            parseobj.addOptional('Prad',[],typeValidation_power);
            parseobj.addOptional('radEff',1,typeValidation_power);

            expected_freqUnit = {'Hz','kHz','MHz','GHz','THz'};
            parseobj.addParameter('freqUnit','Hz', @(x) any(validatestring(x,expected_freqUnit)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            if nargin == 0
                [name,path] = uigetfile('*.txt');
                pathName = [path,name];
            end
            parseobj.parse(pathName,varargin{:})
            
            pathName = parseobj.Results.pathname;
            freq = parseobj.Results.freq;
            Prad = parseobj.Results.Prad;
            radEff = parseobj.Results.radEff;
            freqUnit = parseobj.Results.freqUnit;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            % Open the data file
            if ~strcmp(pathName(end-3:end),'.txt')
                pathName = [pathName,'.txt'];
            end
            fid = fopen(pathName);
            if (fid==-1)
                error(['Unable to open data file ', fileName, '!']);
            end
            
            header = fgetl(fid);
            fgetl(fid);
            Ncol = 8;
            data = fscanf(fid,'%f%f%f%f%f%f%f%f');
            fclose(fid);
            
            % Reshape the data into a matrix
            Nrow = length(data)/Ncol;
            data = reshape(data,Ncol,Nrow).';
            xIn = deg2rad(data(:,2));
            yIn = deg2rad(data(:,1));
            
            headerCell = split(header,']');
            
            % figure out the grid type
            gridMark = strtok(headerCell{1},' [');
            switch gridMark
                case 'Theta'
                    gridTypeIn = 'PhTh';
                    coorTypeIn = [];
                case 'Elev.'
                    gridTypeIn = 'AzEl';
                    coorTypeIn = 'Ludwig2AE';
                case 'Alpha'
                    gridTypeIn = 'ElAz';
                    coorTypeIn = 'Ludwig2EA';
            end
            
            % Figure out the coor and pol Type
            coorMark = strtok(headerCell{4},' [');
            if contains(coorMark,'Theta')
                coorTypeIn = 'spherical';
                polTypeIn = 'linear';
            elseif contains(coorMark,'Horiz')
                coorTypeIn = 'Ludwig3';
                polTypeIn = 'linear';
            elseif contains(coorMark,'Left')
                polTypeIn = 'circular';
                % In this case we are actually not sure if it is due to
                % Ludwig3 or spherical - so just select spherical as
                % default
                if isempty(coorTypeIn)
                    coorTypeIn = 'spherical';
                end
            elseif contains(coorMark,'Azimu')
                coorTypeIn = 'Ludwig2AE';
                polTypeIn = 'linear';
            elseif contains(coorMark,'Alpha')
                coorTypeIn = 'Ludwig2EA';
                polTypeIn = 'linear';
            else
                error('Only linear and circular polarizations implemented at this stage')
            end
            
            % Figure out the output type
            [typeMark,scaleMark] = strtok(headerCell{3},' [');
            if contains(typeMark,'Abs(V')
                outTypeIn = 'E-pattern';
            else
                error('Only the E-pattern output type supported at this stage')
            end
            
            % Figure out the scale
            if contains(scaleMark,'dB')
                E1In = lin20(data(:,4));
                E2In = lin20(data(:,6));
            else
                E1In = data(:,4);
                E2In = data(:,6);
            end
            E1In = E1In.*exp(1i.*deg2rad(data(:,5)));
            E2In = E2In.*exp(1i.*deg2rad(data(:,7)));
            
            FF = FarField(xIn,yIn,E1In,E2In,freq,Prad,radEff,...
                'coorType',coorTypeIn,'polType',polTypeIn,'gridType',gridTypeIn,'freqUnit',freqUnit,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
            
            % Copy in the missing redundant edge cut
            if strcmp(FF.yRangeType,'180')
                xvalCopy = min(FF.xRange);
                xvalPaste = xvalCopy + 2*pi;
                FF = FF.copyAndInsertXcut(xvalCopy,xvalPaste);
            elseif strcmp(FF.yRangeType,'360')
                yvalCopy = min(FF.yRange);
                yvalPaste = yvalCopy + 2*pi;
                FF = FF.copyAndInsertYcut(yvalCopy,yvalPaste);
            end
            
            FF.Prad = FF.pradInt;
            FF.radEff = FF.Prad./0.5;   % Default CST power is 0.5 W - this could be wrong in some cases if the CST power is changed in the simulation
            
        end
        
        function FF = readNFSscan(pathName,varargin)
            % READNFSSCAN Create a FarField object from a NFS spherical range scan .txt file.
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readNFSscan';
            
            typeValidator_pathName = @(x) isa(x,'char');
            parseobj.addRequired('pathname',typeValidator_pathName);
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            if nargin == 0
                [name,path] = uigetfile('*.txt');
                pathName = [path,name];
            end
            parseobj.parse(pathName,varargin{:})
            
            pathName = parseobj.Results.pathname;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            
            % Open the data file
            if ~strcmp(pathName(end-3:end),'.txt')
                pathName = [pathName,'.txt'];
            end
            fid = fopen(pathName);
            if (fid==-1)
                error(['Unable to open data file ', fileName, '!']);
            end
            
            % Set up markers in the text
            directivityMarker = 'Directivity';
            farFieldInfoMarker = 'Far-field display setup';
            freqInfoMarker = 'BeamFrequency';
            
            % Read the field header
            while 1
                a = fgetl(fid);
                if strncmp(a,directivityMarker,11), break; end
            end
            D = textscan(a,'%s%s%f%s');
            directivity_dB = D{3}; 
            % Skip over a few to the farfield data
            while 1
                a = fgetl(fid);
                if strncmp(a,farFieldInfoMarker,24), break; end
            end
            % Get the grid information
            fgetl(fid); % Skip th/az heading line
            a = fgetl(fid); % Read the first dimension data line1
            D1 = textscan(a,'%s%s%f%s%s%s%f%s%s%s%f');
            a = fgetl(fid); % Read the first dimension data line2
            D2 = textscan(a,'%s%f%s%s%s%f%s%s%s%f%s');
            [x1span,x1center,x1N] = deal(D1{3},D1{7},D1{11});
            [x1start,x1stop,x1delta] = deal(D2{2},D2{6},D2{10});
            
            fgetl(fid); % Skip ph/el heading line
            a = fgetl(fid); % Read the second dimension data line1
            D1 = textscan(a,'%s%s%f%s%s%s%f%s%s%s%f');
            a = fgetl(fid); % Read the second dimension data line2
            D2 = textscan(a,'%s%f%s%s%s%f%s%s%s%f%s');
            [x2span,x2center,x2N] = deal(D1{3},D1{7},D1{11});
            [x2start,x2stop,x2delta] = deal(D2{2},D2{6},D2{10});
            
            fgetl(fid); % Skip rotation line
            fgetl(fid); % Skip interpolation line
            a = fgetl(fid); % Get coordinate system line
            D = textscan(a,'%s%s%s%s%s%s%s%s%s');
            switch D{3}{1}(1:end-1)
                case 'Th-Phi'
                    gridType = 'PhTh';
                    dataMarker = 'Theta(deg)Phi(deg)';
                case 'Az/El'
                    gridType = 'AzEl';
                    dataMarker = 'Azimuth(deg)Elevation(deg)';
                otherwise
                    error(['Unknown grid type found: ',D{3}{1}])
            end
            switch D{5}{1}
                case 'L2'
                    switch D{6}{1}(1:end-1)
                        case 'Az/El'
                            coorType = 'Ludwig2AE';
                        case 'Eth-Eph'
                            coorType = 'spherical';
                    end
                otherwise
                    error(['Unknown coordinate type found: ',D{5}{1}])
            end
            % Skip over a few to the frequency info
            while 1
                a = fgetl(fid);
                if strncmp(strrep(a,' ',''),freqInfoMarker,13), break; end
            end
            fgetl(fid); % Skip underline
            a = fgetl(fid); % Get frequency line
            D = textscan(a,'%f%f%s%s%s%s%s');
            freq = D{2};
            freqUnit = D{3}{1};
            % Skip over until the actual data
            while 1
                a = fgetl(fid);
                if strncmp(strrep(a,' ',''),dataMarker,length(dataMarker)), break; end
            end
            % Read first pol data
            DATA1 = fscanf(fid,'%f%f%f%f',[4, x1N*x2N]).';
            % Skip over until the next data
            while 1
                a = fgetl(fid);
                if strncmp(strrep(a,' ',''),dataMarker,length(dataMarker)), break; end
            end
            % Read second pol data
            DATA2 = fscanf(fid,'%f%f%f%f',[4, x1N*x2N]).';
            
            % Get angles
%             x = deg2rad(DATA1(:,1));
%             y = deg2rad(DATA1(:,2));
            
%             % Sort and format the data
            iPhaseChange = [];
            phaseChange = 0;
            switch gridType
                case 'PhTh'
%                     if x2span > 180 && x1start < 0 && x1stop > 0 % Don't have full sphere, but some data repeats
%                         % Remove the negative part
%                         % TODO: Check if this always works - not sure if
%                         % the negative side is always the continuous part
%                         % Must also still check to only remove the parts
%                         % that repeat
%                         iRemove = find(DATA1(:,1) < 0);
%                         DATA1(iRemove,:) = [];
%                         DATA2(iRemove,:) = [];
%                     end
                    % th=0 is continuous to the negative side, so shift
                    % by 180 degrees
                    iPhaseChange = find(DATA1(:,1) == 0);
                    phaseChange = 180;
                    y = deg2rad(DATA1(:,1));
                    x = deg2rad(DATA1(:,2));
                case 'AzEl'
                    x = deg2rad(DATA1(:,1));
                    y = deg2rad(DATA1(:,2));
            end
            switch coorType
                case {'spherical','Ludwig2AE'}
                    DATA1(iPhaseChange,4) = DATA1(iPhaseChange,4) + phaseChange;
                    DATA2(iPhaseChange,4) = DATA2(iPhaseChange,4) + phaseChange;
                    E1 = lin20(DATA1(:,3)).*exp(1i.*deg2rad(DATA1(:,4)));
                    E2 = lin20(DATA2(:,3)).*exp(1i.*deg2rad(DATA2(:,4)));
            end
            
            % Calculate the power and scaling
            D0 = lin10(directivity_dB);
            directivityDir = [0,0];     % Assume this thing gives directivity at broadside
            iDir = intersect(find(x == directivityDir(1)),find(y == directivityDir(2)));
            if ~isempty(iDir)
                Prad = 4*pi*r.^2.*(abs(E1(iDir(1))).^2 + abs(E2(iDir(1))).^2)./(2.*376.7303.*D0);
            else
                % Just use the maximum for directivity
                Prad = 4*pi*r.^2.*max(abs(E1).^2 + abs(E2).^2)./(2.*376.7303.*D0);
            end
            
            % Create the FarField object
            polType = 'linear'; % Probably always the case?
            radEff = ones(size(freq)); %replace with manual radiation efficiency calculation
            
            FF = FarField(x,y,E1,E2,freq,Prad,radEff,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
            
%             % Set to a standard grid - try to get positive x range and 180
%             % y range
%             FF = FF.setYrange(180);
%             % Insert the missing cut if required
%             % TODO: Remove after setYrange and setXrange have been fixed to
%             % be more general
%             FF = FF.copyAndInsertXcut(-pi,pi);
%             FF = FF.setXrange('pos');
%             FF = FF.setBase;
        end
        
        function FF = readFITS(inputStruct,gridType,coorType,polType,varargin)
            % READFITS Create a FarField object from a FITS image file.
            
            % Expect the data format as: data([Nx,Ny,Nf,Ne,Nc,Np])
            %   Nx - number of x values
            %   Ny - number of y values
            %   Nf - number of frequencies
            %   Ne - number of patterns in the data
            %   Nc - number of field components [1 | 2]
            %   Np - number of complex data parts [1 | 2] (if this is just
            %   1, the other complex component is expected from another
            %   fits file - see inputStruct.type)
            % inputStruct is a struct of important input parameters,
            % containing a subset of:
            %   .pathName1 - the path and name of the first component (required)
            %   .pathName2 - the path and name of the second component (optional)
            %   .type1 - the type of the first component {('real') | 'imag' | 'mag' | 'phase'} (optional)
            %   .type2 - the type of the second component {('real') | 'imag' | 'mag' | 'phase'} (optional)
            %   .scale1 - the scale of the first component {('lin') | 'dB' | ('rad') | 'deg'} (optional)
            %   .scale2 - the scale of the second component {('lin') | 'dB' | ('rad') | 'deg'} (optional)
            %   .scaleFuncGrid - can be a function handle or scaling factor
            %   for the grid. if it is an array with 2 elements, the first
            %   is used for x, and the second for y
            
            % ToDo: This is very hard to generalise - learn as we go and
            % get more variations I guess
            
            % Parsing through the inputs
            parseobj = inputParser;
            parseobj.FunctionName = 'readFITS';
            
            typeValidator_inputStruct = @(x) isa(x,'struct') || isa(x,'char');
            parseobj.addRequired('inputStruct',typeValidator_inputStruct);
            
            parseobj.addRequired('gridType');
            parseobj.addRequired('coorType');
            parseobj.addRequired('polType');
            
            typeValidation_range = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,2]},'readFITS');
            parseobj.addParameter('xRange',[0,360],typeValidation_range);
            parseobj.addParameter('yRange',[0,180],typeValidation_range);
            parseobj.addParameter('fRange',[1e9,2e9],typeValidation_range);
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR0','BOR1'};
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readFITS');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readFITS');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readFITS');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.addParameter('freqUnit','Hz');
            
            % ToDo - add keyboard functionality to read the second
            % component as well as type/scale
            if nargin == 0 || isempty(inputStruct)
                [name,path] = uigetfile('*.fits','Component 1:');
                inputStruct.pathName1 = [path,name];
            end
            parseobj.parse(inputStruct,gridType,coorType,polType,varargin{:})
            
            inputStruct = parseobj.Results.inputStruct;
            gridType = parseobj.Results.gridType;
            coorType = parseobj.Results.coorType;
            polType = parseobj.Results.polType;
            xRange = parseobj.Results.xRange;
            yRange = parseobj.Results.yRange;
            fRange = parseobj.Results.fRange;
            symmetryXZ = parseobj.Results.symmetryXZ;
            symmetryYZ = parseobj.Results.symmetryYZ;
            symmetryXY = parseobj.Results.symmetryXY;
            symmetryBOR = parseobj.Results.symmetryBOR;
            r = parseobj.Results.r;
            orientation = parseobj.Results.orientation;
            earthLocation = parseobj.Results.earthLocation;
            time = parseobj.Results.time;
            freqUnit = parseobj.Results.freqUnit;
            
            eta0 = 3.767303134749689e+02;
            
            pathName1 = inputStruct.pathName1;
            
            % Set some defaults
            [type1,scale1] = deal('real','lin');
            [type2,scale2] = deal('imag','lin');
            [pathName2] = deal([]);
            if isfield(inputStruct,'type1') && ~isempty(inputStruct.type1)
                type1 = inputStruct.type1;
            end
            if isfield(inputStruct,'scale1') && ~isempty(inputStruct.scale1)
                scale1 = inputStruct.scale1;
            end
            if isfield(inputStruct,'pathName2') && ~isempty(inputStruct.pathName2)
                pathName2 = inputStruct.pathName2;
            end
            if isfield(inputStruct,'type2') && ~isempty(inputStruct.type2)
                type2 = inputStruct.type2;
            else
                if ~isempty(pathName2), type2 = 'imag'; end
            end
            if isfield(inputStruct,'scale2') && ~isempty(inputStruct.scale2)
                scale2 = inputStruct.scale2;
            else
                if ~isempty(pathName2), scale2 = 'lin'; end
            end
            if isfield(inputStruct,'scaleFuncGrid') && ~isempty(inputStruct.scaleFuncGrid)
                scaleFuncGrid = inputStruct.scaleFuncGrid;
            else
                scaleFuncGrid = 1;
            end
            if isempty(gridType), gridType = 'PhTh'; end
            if isempty(coorType), coorType = 'spherical'; end
            if isempty(polType), polType = 'linear'; end
            
            % Test the formats - not complete, but should catch most issues
            if strcmp(type1,'real')
                assert(strcmp(scale1,'lin'),'For cartesian field components the scale must be lin');
                if ~isempty(type2)
                    assert(strcmp(type2,'imag'),'Cannot mix cartesian and polar forms for complex field');
                end
            elseif strcmp(type1,'mag')
                assert(strcmp(scale1,{'lin','dB10','dB20'}),'Scale must be lin, dB10 or dB20 for magnitude types');
                if ~isempty(type2)
                    assert(strcmp(type2,'phase'),'Cannot mix cartesian and polar forms for complex field');
                end
            else
                error('type1 must be real or mag')
            end
            if ~isempty(type2)
                if strcmp(type2,'imag')
                    assert(strcmp(scale2,'lin'),'For cartesian field components the scale must be lin');
                elseif strcmp(type2,'phase')
                    assert(strcmp(scale2,{'deg','rad'}),'Scale must be deg or rad for phase types');
                else
                    error('type2 must be imag or phase')
                end
            end
            
            % sort out extensions
            [~,~,ext1] = fileparts(pathName1);
            if isempty(ext1), pathName1 = [pathName1,'.fits']; end
            if ~isempty(pathName2)
                [~,~,ext2] = fileparts(pathName2);
                if isempty(ext2), pathName2 = [pathName2,'.fits']; end
            end
            
            % Read the info and data and process the fields
            info1 = fitsinfo(pathName1);
            [xvect1,yvect1,fvect1] = fitsGrid(info1,xRange,yRange,fRange);
            if isa(scaleFuncGrid,'function_handle')
                [X,Y] = meshgrid(scaleFuncGrid(xvect1),scaleFuncGrid(yvect1));
            else
                [X,Y] = meshgrid(scaleFuncGrid.*xvect1,scaleFuncGrid.*yvect1);
            end
            [Nx,Ny,Nf] = deal(numel(xvect1),numel(yvect1),numel(fvect1));
            data1 = fitsread(pathName1);
            
            % Ne = number of fields; Nc = number of components; Np = number
            % of complex parts
            [~,~,~,Ne,Nc,Np] = size(data1(1,1,1,:,:,:)); % Can be more than one field in here - column 4 represents number of patterns. 
            % E-field naming: EXy; X = component, y = complex part
            E1a = (data1(:,:,:,:,1,:));
            [E1b,E2a,E2b] = deal([]);
            if Nc > 1, E2a = (data1(:,:,:,:,2,:)); end
            % Sort out type and scale
            switch scale1
                case 'dB10'
                    E1a = lin10(E1a);
                    if Nc > 1, E2a = lin10(E2a); end
                case 'dB20'
                    E1a = lin20(E1a);
                    if Nc > 1, E2a = lin20(E2a); end
            end
            % Get the second complex parameter
            if ~isempty(pathName2) && Np == 1
                info2 = fitsinfo(pathName2);
                [xvect2,yvect2,fvect2] = fitsGrid(info2,xRange,yRange,fRange);
                data2 = fitsread(pathName2);
                assert(all(xvect1 == xvect2)&&all(yvect1 == yvect2)&&all(fvect1 == fvect2)&&all(size(data1)==size(data2)),'The two patterns must have identical grids and number of fields')
                E1b = (data2(:,:,:,:,1));
                if Nc > 1, E2b = (data2(:,:,:,:,2)); end
            elseif Np == 2
                E1b = data1(:,:,:,:,1,2);
                if Nc > 1, E2b = data1(:,:,:,:,2,2); end
            end
            E1 = E1a(:,:,:,:,1,1);
            if Nc > 1, E2 = E2a(:,:,:,:,1,1); end
            
            % Handle the complex number types
            switch type2
                case 'imag'
                    if ~isempty(E1b), E1 = E1 + 1i.*E1b; end
                    if Nc > 1
                        if ~isempty(E2b), E2 = E2 + 1i.*E2b; end
                    end
                case 'phase'
                    if strcmp(scale2,'deg')
                        E1b = deg2rad(E1b);
                        if Nc > 1
                            if ~isempty(E2b), E2b = deg2rad(E2b); end
                        end
                    end
                    if ~isempty(E1b), E1 = E1.*exp(1i.*E1b); end
                    if Nc > 1
                        if ~isempty(E2b), E2 = E2.*exp(1i.*E2b); end
                    end
            end
                
            if isempty(E1b), coorType = 'power'; end     % Force this, since only one component provided
            
            freq = fvect1;
            Prad = [];
            radEff = 1;
            FF(1:Ne) = FarField;
            for ee = 1:Ne
                FF(ee) = FarField(X(:),Y(:),reshape(E1(:,:,:,ee,1),Nx*Ny,Nf),reshape(E2(:,:,:,ee,1),Nx*Ny,Nf),freq,Prad,radEff,...
                'coorType',coorType,'polType',polType,'gridType',gridType,'freqUnit',freqUnit,...
                'symmetryXZ',symmetryXZ,'symmetryYZ',symmetryYZ,'symmetryXY',symmetryXY,'symmetryBOR',symmetryBOR,...
                'r',r,'orientation',orientation,'earthLocation',earthLocation,'time',time);
            end
            
            function [xvect,yvect,fvect] = fitsGrid(info,xRange,yRange,fRange)
                % Internal helper function to set up the grids
                % Assume an [x,y,freq] format...
                keywords = info.PrimaryData.Keywords;
                Nx = cell2mat(keywords(strcmp(keywords(:,1),'NAXIS1'),2));
                Ny = cell2mat(keywords(strcmp(keywords(:,1),'NAXIS2'),2));
                Nf = cell2mat(keywords(strcmp(keywords(:,1),'NAXIS3'),2));
                if ismember('CRVAL1',keywords(:,1)) && ismember('CRPIX1',keywords(:,1)) && ismember('CDELT1',keywords(:,1))
                    startValX  = cell2mat(keywords(strcmp(keywords(:,1),'CRVAL1'),2));
                    startPosX = cell2mat(keywords(strcmp(keywords(:,1),'CRPIX1'),2));
                    delX = cell2mat(keywords(strcmp(keywords(:,1),'CDELT1'),2));
                    xvect = (0:delX:(Nx-1)*delX) + startValX - (startPosX-mod(Nx,2))*delX;
                else
                    xvect = linspace(xRange(1),xRange(2),Nx);
                end
                if ismember('CRVAL2',keywords(:,1)) && ismember('CRPIX2',keywords(:,1)) && ismember('CDELT2',keywords(:,1))
                    startValY  = cell2mat(keywords(strcmp(keywords(:,1),'CRVAL2'),2));
                    startPosY = cell2mat(keywords(strcmp(keywords(:,1),'CRPIX2'),2));
                    delY = cell2mat(keywords(strcmp(keywords(:,1),'CDELT2'),2));
                    yvect = (0:delY:(Ny-1)*delY) + startValY - (startPosY-mod(Ny,2))*delY;
                else
                    yvect = linspace(yRange(1),yRange(2),Ny);
                end
                if ismember('CRVAL3',keywords(:,1)) && ismember('CRPIX3',keywords(:,1)) && ismember('CDELT3',keywords(:,1))
                    startValF  = cell2mat(keywords(strcmp(keywords(:,1),'CRVAL3'),2));
                    startPosF = cell2mat(keywords(strcmp(keywords(:,1),'CRPIX3'),2));
                    delF = cell2mat(keywords(strcmp(keywords(:,1),'CDELT3'),2));
                    fvect = (0:delF:(Nf-1)*delF) + startValF - (startPosF-mod(Nf,2))*delF;
                else
                    fvect = linspace(fRange(1),fRange(2),Nf);
                end
            end
        end
        
        function FF = farFieldFromPowerPattern(x,y,U,freq,varargin)
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
            
            expected_fieldPol = {'linearX','linearY','circularLH','circularRH','power'};
            parseobj.addParameter('fieldPol','linearY', @(x) any(validatestring(x,expected_fieldPol)));
            
            expected_gridType = {'PhTh','DirCos','AzEl','ElAz','Horiz','TrueView','ArcSin','Mollweide','RAdec','GalLongLat'};
            parseobj.addParameter('gridType','PhTh', @(x) any(validatestring(x,expected_gridType)));
            
            expected_symPlane = {'none','electric','magnetic'};
            parseobj.addParameter('symmetryXZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryYZ','none', @(x) any(validatestring(x,expected_symPlane)));
            parseobj.addParameter('symmetryXY','none', @(x) any(validatestring(x,expected_symPlane)));
            
            expected_symBOR = {'none','BOR1'};  % Can't do BOR0 symmetry
            parseobj.addParameter('symmetryBOR','none', @(x) any(validatestring(x,expected_symBOR)));
            
            typeValidation_scalar = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','scalar'},'readGRASPcut');
            parseobj.addParameter('r',1,typeValidation_scalar);
            
            typeValidation_orientation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('orientation',[0,0,0],typeValidation_orientation);
            
            typeValidation_earthLocation = @(x) validateattributes(x,{'numeric'},{'real','finite','nonnan','size',[1,3]},'readGRASPcut');
            parseobj.addParameter('earthLocation',[deg2rad(18.86) deg2rad(-33.93) 300],typeValidation_earthLocation);
            
            typeValidation_time = @(x) isa(x,'datetime');
            parseobj.addParameter('time',datetime(2018,7,22,0,0,0),typeValidation_time);
            
            parseobj.parse(x,y,U,freq,varargin{:})
            
            x = parseobj.Results.x;
            y = parseobj.Results.y;
            U = parseobj.Results.P;
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
            
            if isscalar(U)
                U = U./(4.*pi).*ones(size(ph));
            end
            %From power pattern and polarization parameters, generate E1 and E2 accordingly
            if strcmp(symmetryBOR,'none')
                coorType = 'Ludwig3';
                switch fieldPol
                    case 'linearX' % linearly polarised along X-axis
                        polType = 'linear';
                        E1  = sqrt(U./r^2.*2*eta0);
                        E2  = zeros(size(U));
                    case 'linearY' % linearly polarised along Y-axis
                        polType = 'linear';
                        E1  = zeros(size(U));
                        E2  = sqrt(U./r^2.*2*eta0);
                    case 'circularLH'  % Lefthand Circular polarization
                        polType = 'circular';
                        E1  = sqrt(U./r^2.*2*eta0);
                        E2  = zeros(size(U));
                    case 'circularRH'  % Righthand Circular polarization
                        polType = 'circular';
                        E1  = zeros(size(U));
                        E2  = sqrt(U./r^2.*2*eta0);
                    case 'power' % Only real power pattern of interest
                        polType = 'none';
                        coorType = 'power';
                        E1 = sqrt(U./r^2.*2*eta0);
                        E2 = [];
                    otherwise
                        error('fieldPol input string unrecognised')
                end
            else
                assert(strcmp(gridType,'PhTh'),'PhTh grid required for BOR1 definition')
                assert(all(abs(unique(x) - [0;pi/2]) < tol),'Invalid range for BOR1 symmetry (E-plane and H-plane required)')
                coorType = 'spherical';
                iph0 = find(abs(x - 0) < tol);
                iph90 = find(abs(x - pi/2) < tol);
                [E1,E2] = deal(zeros(size(U)));
                switch fieldPol
                    case 'linearX' 
                        polType = 'linear';
                        E1(iph0,:)  = sqrt(U(iph0,:)./r^2.*2*eta0);
                        E2(iph90,:)  = -sqrt(U(iph90,:)./r^2.*2*eta0);
                    case 'linearY'
                        polType = 'linear';
                        E1(iph90,:)  = sqrt(U(iph90,:)./r^2.*2*eta0);
                        E2(iph0,:)  = sqrt(U(iph0,:)./r^2.*2*eta0);
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
    
    methods (Access = private)
        
        %% Base representation methods
        function obj = setBase(obj)
            % SETBASE Set base grid and fields
            
            obj = obj.setBaseGrid;
            obj = obj.setBaseFields;
        end
        
        function obj = setBaseGrid(obj)
            % SETBASEGRID Set base grid
            
            obj.xBase = obj.x;
            obj.yBase = obj.y;
            obj.gridTypeBase = obj.gridType;
            % Set after the gridType hase been set for the check in the
            % dependent property getter. ph and th not stored for astro
            % grid bases
            [obj.phBase,obj.thBase] = obj.getPhTh; % Faster than the 2 dependent getters
        end
        
        function obj = setBaseFields(obj)
            % SETBASEFIELDS Set base fields
            
            obj.E1Base = obj.E1;
            obj.E2Base = obj.E2;
            obj.E3Base = obj.E3;
            obj.coorTypeBase = obj.coorType;
            obj.polTypeBase = obj.polType;
        end
        
        function obj = clearBase(obj)
            % CLEARBASE Clear base grid and fields
            
            obj = obj.clearBaseGrid;
            obj = obj.clearBaseFields;
        end
        
        function obj = clearBaseGrid(obj)
            % CLEARBASEGRID Clear base grid
            
            obj.xBase = [];
            obj.yBase = [];
            obj.gridTypeBase = [];
            obj.phBase = [];
            obj.thBase = [];
        end
        
        function obj = clearBaseFields(obj)
            % CLEARBASEFIELDS Clear the base fields
            
            obj.E1Base = [];
            obj.E2Base = [];
            obj.E3Base = [];
            obj.coorTypeBase = [];
            obj.polTypeBase = [];
        end
        
        function [FF1,FF2,level] = mathSetup(obj1,obj2)
            % MATHSETUP returns 2 equal sampled FarFields for math
            % operations
            
            % level = 0: current grids used
            %         1: base grids used
            %         2: resampled onto new grid
            
            % Check if current grids and types are equal
            currentEqual = [isGridEqual(obj1,obj2),typesAreEqual(obj1,obj2)];
            level = 0;  
            if ~all(currentEqual)   % Try the bases
                obj1 = obj1.reset2Base;
                obj2 = obj2.reset2Base;
                level = level + 1;
                baseEqual = [isGridEqual(obj1,obj2),typesAreEqual(obj1,obj2)];
                if ~all(baseEqual)  % Need to resample
                    % TODO - check for cases where the grid ranges are not
                    % compatible
                    obj2 = transformTypes(obj2, obj1); % Make sure the types are equal
                    level = level + 1;
                    if obj1.freqHz == obj2.freqHz
                        [E1i,E2i,E3i] = deal(zeros(obj1.Nang,obj2.Nf));
                        for ff = 1:obj2.Nf
                            E1i(:,ff) = interpolateGrid(obj2,'E1',obj1.x,obj1.y);
                            if numel(obj2.E2) > 0 && numel(obj1.E1) > 0
                                E2i(:,ff) = interpolateGrid(obj2,'E2',obj1.x,obj1.y);
                            else
                                E2i = [];
                            end
                            if numel(obj2.E3) > 0 && numel(obj1.E3) > 0
                                E3i(:,ff) = interpolateGrid(obj2,'E3',obj1.x,obj1.y);
                            else
                                E3i = [];
                            end
                        end
                        obj2 = FarField(obj1.x,obj1.y,E1i,E2i,obj2.freq,obj2.Prad,obj2.radEff,...
                            'coorType',obj1.coorType,'polType',obj1.polType,'gridType',obj1.gridType,'freqUnit',obj1.freqUnit,...
                            'symmetryXZ',obj1.symmetryXZ,'symmetryYZ',obj1.symmetryYZ,'symmetryXY',obj1.symmetryXY,'symmetryBOR',obj1.symmetryBOR,'E3',E3i,...
                            'r',obj1.r,'orientation',obj1.orientation,'earthLocation',obj1.earthLocation,'time',obj1.time);
                    else
                        % TODO - obj2 frequency resample on obj1.freqHz
                        error('mathSetup not implemented for unequal base frequencies yet')
                    end
                end
            end
            [FF1,FF2] = deal(obj1,obj2);
        end
        
        %% Grid getters
        function [u, v, w] = getDirCos(obj)
            % GETDIRCOS Get DirCos grid.
            
            switch obj.gridType
                case 'DirCos'
                    u = obj.x;
                    v = obj.y;
                    w = sqrt(1 - u.^2 - v.^2);
                otherwise
                    handle2DirCos = str2func([obj.gridType,'2DirCos']);
                    [u,v,w] = handle2DirCos(obj.x,obj.y);
            end
        end
        
        function [Xg, Yg] = getTrueView(obj)
            % GETTRUEVIEW Get TrueView grid.
            
            switch obj.gridType
                case 'TrueView'
                    Xg = obj.x;
                    Yg = obj.y;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [Xg,Yg] = DirCos2TrueView(u,v,w);
            end
        end
        
        function [asinu, asinv] = getArcSin(obj)
            % GETARCSIN Get ArcSin grid.
            
            switch obj.gridType
                case 'ArcSin'
                    asinu = obj.x;
                    asinv = obj.y;
                otherwise
                    [u,v,w] = getDirCos(obj);
                    [asinu,asinv] = DirCos2ArcSin(u,v,w);
            end
        end
        
        function [ph, th] = getPhTh(obj)
            % GETPHTH Get PhTh grid.
            
            if strcmp(obj.gridType,'PhTh')
                ph = obj.x;
                th = obj.y;
            else
                if any(strcmp(obj.gridType,obj.astroGrids))
                    if strcmp(obj.gridType,'Horiz')
                        az = obj.x;
                        alt = obj.y;
                    else
                        handle2Horiz = str2func([obj.gridType,'2Horiz']);
                        [az,alt] = handle2Horiz(obj.x,obj.y,obj.julDate,obj.earthLocation(1:2));
                    end
                    [u,v,w] = Horiz2DirCos(az,alt);
                else
                    [u,v,w] = getDirCos(obj);
                end
                [ph,th] = DirCos2PhTh(u,v,w);
            end
        end
        
        function [az, el] = getAzEl(obj)
            % GETAZEL Get AzEl grid.
            
            if strcmp(obj.gridType,'AzEl')
                az = obj.x;
                el = obj.y;
            else
                if any(strcmp(obj.gridType,obj.astroGrids))
                    if strcmp(obj.gridType,'Horiz')
                        azH = obj.x;
                        alt = obj.y;
                    else
                        handle2Horiz = str2func([obj.gridType,'2Horiz']);
                        [azH,alt] = handle2Horiz(obj.x,obj.y,obj.julDate,obj.earthLocation(1:2));
                    end
                    [u,v,w] = Horiz2DirCos(azH,alt);
                else
                    [u,v,w] = getDirCos(obj);
                end
                [az,el] = DirCos2AzEl(u,v,w);
            end
        end
        
        function [ep, al] = getElAz(obj)
            % GETELAZ Get ElAz grid.
            
            if strcmp(obj.gridType,'ElAz')
                ep = obj.x;
                al = obj.y;
            else
                if any(strcmp(obj.gridType,obj.astroGrids))
                    if strcmp(obj.gridType,'Horiz')
                        azH = obj.x;
                        alt = obj.y;
                    else
                        handle2Horiz = str2func([obj.gridType,'2Horiz']);
                        [azH,alt] = handle2Horiz(obj.x,obj.y,obj.julDate,obj.earthLocation(1:2));
                    end
                    [u,v,w] = Horiz2DirCos(azH,alt);
                else
                    [u,v,w] = getDirCos(obj);
                end
                [ep,al] = DirCos2ElAz(u,v,w);
            end
        end
        
        function [az,alt] = getHoriz(obj)
            % GETHoriz Get Horiz grid.
            
            switch obj.gridType
                case 'Horiz'
                    az = obj.x;
                    alt = obj.y;
                case [obj.localGrids,obj.projectionGrids]
                    [u,v,w] = getDirCos(obj);
                    [az,alt] = DirCos2Horiz(u,v,w);
                otherwise
                    handle2Horiz = str2func([obj.gridType,'2Horiz']);
                    [az,alt] = handle2Horiz(obj.x,obj.y,obj.julDate,obj.earthLocation(1:2));
            end
        end
        
        function [RA, dec] = getRAdec(obj)
            % GETRADEC Get RAdec grid.
            
            switch obj.gridType
                case 'RAdec'
                    RA = obj.x;
                    dec = obj.y;
                case obj.projectionGrids
                    [u,v,w] = getDirCos(obj);
                    [RA,dec] = DirCos2RAdec(u,v,w);
                case obj.localGrids
                    [az,alt] = PhTh2Horiz(obj.ph,obj.th);
                    [RA,dec] = Horiz2RAdec(az,alt,obj.julDate,obj.earthLocation(1:2));
                case 'Horiz'
                    [RA,dec] = Horiz2RAdec(obj.x,obj.y,obj.julDate,obj.earthLocation(1:2));
                case 'GalLongLat'
                    [RA,dec] = GalLongLat2RAdec(obj.x,obj.y);
            end
        end
        
        function [long, lat] = getGalLongLat(obj)
            % GETGALLINGLAT Get GalLongLat grid.
            
            switch obj.gridType
                case 'GalLongLat'
                    long = obj.x;
                    lat = obj.y;
                case obj.projectionGrids
                    [u,v,w] = getDirCos(obj);
                    [long,lat] = DirCos2GalLongLat(u,v,w);
                case obj.localGrids
                    [az,alt] = PhTh2Horiz(obj.ph,obj.th);
                    [long,lat] = Horiz2GalLongLat(az,alt,obj.julDate,obj.earthLocation(1:2));
                case 'Horiz'
                    [long,lat] = Horiz2GalLongLat(obj.x,obj.y,obj.julDate,obj.earthLocation(1:2));
                case 'RAdec'
                    [long,lat] = RAdec2GalLongLat(obj.x,obj.y);
            end

        end
        
        %% Coordinate system getters
        function [Eth, Eph, Er] = getEspherical(obj)
            % GETESPHERICAL Get Espherical coordinates.
            
            tol = 10^(-obj.nSigDig);
            % Get the base Ph and Th - changed to base grid externally
            [Ph,Th] = getPhTh(obj); % Faster than using the dependent variables - just called once
            TH = repmat(Th(:,1),1,obj.Nf);
            PH = repmat(Ph(:,1),1,obj.Nf);
            switch obj.coorType
                case 'spherical'
                    Eth = obj.E1;
                    Eph = obj.E2;
                case 'Ludwig1'
                    Eth = cos(TH).*cos(PH).*obj.E1 + cos(TH).*sin(PH).*obj.E2 - sin(TH).*obj.E3;
                    Eph = -sin(PH).*obj.E1 + cos(PH).*obj.E2;
                case 'Ludwig2AE'
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Del = cos(PH).^2 + cos(TH).^2.*sin(PH).^2;
                    Eth = (cosEl./Del).*(cos(PH).*obj.E1 + cos(TH).*sin(PH).*obj.E2);
                    Eph = (cosEl./Del).*(-cos(TH).*sin(PH).*obj.E1 + cos(PH).*obj.E2);
                    % Pole at th = +-90; ph = +-90
                    polePos = ((abs(abs(wrap2pi(TH))-pi/2) < tol) & (abs(abs(wrap2pi(PH))-pi/2) < tol));
                    Eth(polePos) = obj.E1(polePos);
                    Eph(polePos) = obj.E2(polePos);
                case 'Ludwig2EA'
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Del = cos(TH).^2.*cos(PH).^2 + sin(PH).^2;
                    Eth = (cosAl./Del).*(cos(TH).*cos(PH).*obj.E1 + sin(PH).*obj.E2);
                    Eph = (cosAl./Del).*(-sin(PH).*obj.E1 + cos(TH).*cos(PH).*obj.E2);
                    % Poles at th = +-90; ph = [-180|0|180|360]
                    % Keep continuous in ph-cuts - get limits along constant phi
                    polePos = ((abs(abs(wrap2pi(TH))-pi/2) < tol) & ((abs(abs(wrap2pi(PH))-0) < tol) | (abs(abs(wrap2pi(PH))-pi) < tol) | (abs(abs(wrap2pi(PH))-2*pi) < tol)));
                    Eth(polePos) = -cos(PH(polePos)).*obj.E1(polePos);
                    Eph(polePos) = cos(PH(polePos)).*obj.E2(polePos);
                case 'Ludwig3'
                    Del = 1;
                    Eth = (1./Del).*(cos(PH).*obj.E1 + sin(PH).*obj.E2);
                    Eph = (1./Del).*(-sin(PH).*obj.E1 + cos(PH).*obj.E2);
            end
            if nargout > 2
                Er = zeros(size(Eth));
            end
        end
        
        function [Ex, Ey, Ez] = getELudwig1(obj)
            % GETELUDWIG1 Get ELudwig1 coordinates.
            
            switch obj.coorType
                case 'Ludwig1'
                    Ex = obj.E1;
                    Ey = obj.E2;
                    if nargout > 2, Ez = obj.E3; end
                otherwise
                    [Eth, Eph, ~] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    % Assume farfield so no radial E-field
                    Ex = cos(TH).*cos(PH).*Eth - sin(PH).*Eph;
                    Ey = cos(TH).*sin(PH).*Eth + cos(PH).*Eph;
                    if nargout > 2, Ez = zeros(size(Ex)); end  % Strict definition in the Ludwig paper 
            end
        end
        
        function [Eaz, Eel, E3] = getELudwig2AE(obj)
            % GETELUDWIG2AE Get ELudwig2AE coordinates.
            
            switch obj.coorType
                case 'Ludwig2AE'
                    Eaz = obj.E1;
                    Eel = obj.E2;
                    if nargout > 2, E3 = obj.E3; end
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    cosEl = sqrt(1 - sin(TH).^2.*sin(PH).^2);
                    Eaz = (1./cosEl).*(cos(PH).*Eth - cos(TH).*sin(PH).*Eph);
                    Eel = (1./cosEl).*(cos(TH).*sin(PH).*Eth + cos(PH).*Eph);
                    if nargout > 2, E3 = zeros(size(Eaz)); end
                    % Take limit in poles along constant ph
                    phPoles = deg2rad([90,270].');
                    poleMat = [ones(2,1).*deg2rad(90),phPoles;ones(2,1).*deg2rad(270),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[abs(obj.th),abs(obj.ph)],'rows');
                    iPole = iPole(iPole>0);
                    Eaz(iPole,:) = sin(TH(iPole)).*sin(PH(iPole)).*Eph(iPole,:);
                    Eel(iPole,:) = -sin(TH(iPole)).*sin(PH(iPole)).*Eth(iPole,:);
                    if nargout > 2, E3(iPole,:) = 0; end
            end
        end
        
        function [Eal, Eep, E3] = getELudwig2EA(obj)
            % GETELUDWIG2EA Get ELudwig2EA coordinates.
            
            switch obj.coorType
                case 'Ludwig2EA'
                    Eal = obj.E1;
                    Eep = obj.E2;
                    if nargout > 2, E3 = obj.E3; end
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    TH = repmat(obj.th(:,1),1,obj.Nf);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    cosAl = sqrt(1 - sin(TH).^2.*cos(PH).^2);
                    Eal = (1./cosAl).*(cos(TH).*cos(PH).*Eth - sin(PH).*Eph);
                    Eep = (1./cosAl).*(sin(PH).*Eth + cos(TH).*cos(PH).*Eph);
                    if nargout > 2, E3 = zeros(size(Eal)); end
                    % Take limit in poles along constant ph
                    phPoles = deg2rad([0,180,360].');
                    poleMat = [ones(3,1).*deg2rad(90),phPoles;ones(3,1).*deg2rad(270),phPoles]; % [th=90,ph]
                    [~,iPole] = ismember(poleMat,[abs(obj.th),abs(obj.ph)],'rows');
                    iPole = iPole(iPole>0);
                    Eal(iPole,:) = -sin(TH(iPole)).*cos(PH(iPole)).*Eth(iPole,:);
                    Eep(iPole,:) = -sin(TH(iPole)).*cos(PH(iPole)).*Eph(iPole,:);
                    if nargout > 2, E3(iPole,:) = 0; end
            end
        end
        
        function [Eh, Ev, E3] = getELudwig3(obj)
            % GETELUDWIG3 Get ELudwig3 coordinates.
            
            switch obj.coorType
                case 'Ludwig3'
                    Eh = obj.E1;
                    Ev = obj.E2;
                    if nargout > 2, E3 = obj.E3; end
                otherwise
                    [Eth, Eph] = getEspherical(obj);
                    PH = repmat(obj.ph(:,1),1,obj.Nf);
                    Eh = cos(PH).*Eth - sin(PH).*Eph;
                    Ev = sin(PH).*Eth + cos(PH).*Eph;
                    if nargout > 2, E3 = zeros(size(Eh)); end
            end
        end
        
        
        %% Polarization type getters
        function [E1lin, E2lin, E3lin] = getElin(obj)
            % GETELIN Get linear polarization.
            
            % This is the base version of the polarization transformers.
            % Everything goes to and from here
            coorTypeIn = obj.coorType;
            coorTypeH = str2func(['coor2',coorTypeIn]);
            obj = obj.reset2Base; % This is done just to change the fields - don't care about the grid
            obj = coorTypeH(obj,false);
            switch obj.polType % Should be the same as the transformed object - can use obj or obj1
                case 'linear'
                    E1lin = obj.E1;
                    E2lin = obj.E2;
                case 'circular'
                    Del = 2*1i;
                    E1lin = sqrt(2)./Del.*(1i.*obj.E1 + 1i.*obj.E2);
                    E2lin = sqrt(2)./Del.*(-obj.E1 + obj.E2);
                case 'slant'
                    PSI = ones(size(obj.E1)).*obj.slant; % Slant of input object
                    Del = 1;
                    E1lin = 1./Del.*(cos(PSI).*obj.E1 + sin(PSI).*obj.E2);
                    E2lin = 1./Del.*(-sin(PSI).*obj.E1 + cos(PSI).*obj.E2);
            end
            E3lin = [];
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
                    E3circ = [];
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
                    E3slant = [];
            end
        end
        
        %% Name setters
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
                case 'Horiz'
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
                case 'none' % power only pattern
                    E1name = 'sqrt(2\eta_0W)';
                    E2name = 'none';
                otherwise
                    error(['Unknown polType property: ', obj.polType]);
            end
        end
        
        %% Range changers
        function obj = range180sym(obj)
            % RANGE180SYM sets the range to the base: 'sym', '180'
            
            assert(any(strcmp(obj.gridType,obj.sphereGrids)),['Can only do range transformations on spherical grids, not on a ',obj.gridType,' grid.'])
            
            tol = 10^(-obj.nSigDig);
            if strcmp(obj.xRangeType,'sym') && strcmp(obj.yRangeType,'180')
                % Do nothing - already there
                return;
            else
                % Get original coordinate system
                coorTypeHandle = str2func(['coor2',obj.coorType]);
                obj = obj.rangeChangeCoorTrans;
                
                % Bookkeeping on the original object
                xSpanOrig = max(obj.x) - min(obj.x);
                ySpanOrig = max(obj.y) - min(obj.y);
                NangOrig = obj.Nang;
                
                if strcmp(obj.xRangeType,'pos') && strcmp(obj.yRangeType,'180')
                    obj.x = wrap2pi(obj.x);
                    [val1,val2] = deal(-pi,pi);
                    % Insert possible missing cuts
                    if abs(xSpanOrig - 2*pi) < tol
                        i1 = find(abs(obj.x - val1) < tol);
                        i2 = find(abs(obj.x - val2) < tol);
                        iin = [i1;i2];
                        if numel(iin) > 0
                            xAdd = [ones(size(i1)).*val2;ones(size(i2)).*val1];
                            yAdd = obj.y(iin);
                            obj = insertMissingCuts(obj,iin,xAdd,yAdd);
                        end
                    end
                elseif strcmp(obj.xRangeType,'pos') && strcmp(obj.yRangeType,'360')
                    xshift = @(x) -pi + x;
                    if strcmp(obj.gridType,'PhTh')
                        iShift = find(obj.y > pi + tol);
                        yshift = @(y) 2*pi - y;
                        obj.x(iShift) = xshift(obj.x(iShift));
                        obj.y(iShift) = yshift(obj.y(iShift));
                        iin = find(abs(obj.y - pi) < tol);
                        poleSign = -1;
                    else
                        obj.y = wrap2pi(obj.y);
                        iShift = find(abs(obj.y) > pi/2 + tol);
                        yshift = @(y) sign(y).*(pi - abs(y));
                        obj.x(iShift) = xshift(obj.x(iShift));
                        obj.y(iShift) = yshift(obj.y(iShift));
                        iin = find(abs(abs(obj.y) - pi/2) < tol);
                        poleSign = -1;
                    end
                    % Fix E-field signs
                    iE = repmat(ismember((1:obj.Nang).',iShift(:)),1,obj.Nf);
                    obj.E1(iE) = -obj.E1(iE);
                    obj.E2(iE) = -obj.E2(iE);
                    if abs(ySpanOrig - 2*pi) < tol
                        xAdd = xshift(obj.x(iin));
                        yAdd = yshift(obj.y(iin));
                        % Try to fix the inserted field sign
                        obj = obj.insertDirs(xAdd,yAdd,poleSign.*obj.E1(iin,:),poleSign.*obj.E2(iin,:));
                    end
                elseif strcmp(obj.xRangeType,'sym') && strcmp(obj.yRangeType,'360')
                    xshift = @(x) x - (sign(x) + ~sign(x)).*pi; % Make the sign(x) = 0 -> 1
                    if strcmp(obj.gridType,'PhTh')
                        iShift = find(obj.y < -tol);
                        yshift = @(y) -y;
                    else
                        iShift = find(abs(obj.y) > pi/2+tol);
                        yshift = @(y) sign(y).*(pi-abs(y));
                    end
                    % Fix E-field signs
                    iE = repmat(ismember((1:obj.Nang).',iShift(:)),1,obj.Nf);
                    obj.E1(iE) = -obj.E1(iE);
                    obj.E2(iE) = -obj.E2(iE);
                    obj.x(iShift) = xshift(obj.x(iShift));
                    obj.y(iShift) = yshift(obj.y(iShift));
                    if abs(ySpanOrig - 2*pi) < tol
                        if strcmp(obj.gridType,'PhTh')
                            iin = find(abs(abs(obj.y) - 0) < tol);
                            poleSign = -1;
                        else
                            iin = find(abs(abs(obj.y) - pi/2) < tol);
                            poleSign = -1;
                        end
                        xAdd = xshift(obj.x(iin));
                        yAdd = yshift(obj.y(iin));
                        % Try to fix the inserted field sign
                        obj = obj.insertDirs(xAdd,yAdd,poleSign.*obj.E1(iin,:),poleSign.*obj.E2(iin,:));
                    end
                    % Add the x = 180 cut
                    if abs(xSpanOrig - pi) < tol
                        iin = find(abs(obj.x + pi) < tol);
                        xAdd = obj.x(iin) + 2*pi;
                        yAdd = obj.y(iin);
                        obj = insertMissingCuts(obj,iin,xAdd,yAdd);
                    end
                    
                    
                end
                % Sort
                obj = obj.sortGrid;
                % Reset coordinate Type
                obj = coorTypeHandle(obj,false);
            end
        end
        
        function obj = range180pos(obj)
            % RANGE180POS sets the range to 'pos', '180'
            
            assert(any(strcmp(obj.gridType,obj.sphereGrids)),['Can only do range transformations on spherical grids, not on a ',obj.gridType,' grid.'])
            
            tol = 10^(-obj.nSigDig);
            
            if strcmp(obj.xRangeType,'pos') && strcmp(obj.yRangeType,'180') 
                % Do nothing - already there
                return;
            else
                % Set to standard base
                obj = obj.rangeChangeCoorTrans;
                
                % Set to standard base
                obj = obj.range180sym;
                % Bookkeeping on the original object
                xSpanOrig = max(obj.x) - min(obj.x);
                
                obj.x = wrap22pi(obj.x);
                [val1,val2] = deal(0,2*pi);
                % Insert possible missing cuts
                if abs(xSpanOrig - 2*pi) < tol
                    i1 = find(abs(obj.x - val1) < tol);
                    i2 = find(abs(obj.x - val2) < tol);
                    iin = [i1;i2];
                    if numel(iin) > 0
                        xAdd = [ones(size(i1)).*val2;ones(size(i2)).*val1];
                        yAdd = obj.y(iin);
                        obj = insertMissingCuts(obj,iin,xAdd,yAdd);
                    end
                end
                % Sort
                obj = obj.sortGrid;
            end
        end
        
        function obj = range360sym(obj)
            % RANGE360SYM sets the range to 'sym', '360'
            
            assert(any(strcmp(obj.gridType,obj.sphereGrids)),['Can only do range transformations on spherical grids, not on a ',obj.gridType,' grid.'])

            
            tol = 10^(-obj.nSigDig);
            
            if strcmp(obj.xRangeType,'sym') && strcmp(obj.yRangeType,'360')
                % Do nothing - already there
                return;
            else
                % Get original coordinate system
                coorTypeHandle = str2func(['coor2',obj.coorType]);
                obj = obj.rangeChangeCoorTrans;
                
                % Set to standard base
                obj = obj.range180sym;
                % Bookkeeping on the base object
                xSpanOrig = diff(obj.xRange);
                ySpanOrig = diff(obj.yRange);
                NangOrig = obj.Nang;
                
                iShift = find(abs(obj.x) > pi/2);
                xshift = @(x) -sign(x).*(pi - abs(x));
                if strcmp(obj.gridType,'PhTh')
                    yshift = @(y) -y;
                    % Grab the required y = 0 cut for later insertion to
                    % ensure a smooth pole
                    ixPoleCuts = find(abs(obj.x) <= pi/2+tol);
                    iy0 = intersect(find(abs(obj.y - 0) < tol),ixPoleCuts);
                    [xy0,yy0,E1y0,E2y0] = deal(obj.x(iy0),obj.y(iy0),obj.E1(iy0,:),obj.E2(iy0,:));
                else
                    yshift = @(y) (sign(y) + ~sign(y)).*pi-y; % Make the sign(y) = 0 -> 1
                end
                % Fix E-field signs
                iE = repmat(ismember((1:obj.Nang).',iShift(:)),1,obj.Nf);
                obj.E1(iE) = -obj.E1(iE);
                obj.E2(iE) = -obj.E2(iE);

                obj.x(iShift) = xshift(obj.x(iShift));
                obj.y(iShift) = yshift(obj.y(iShift));
                if abs(xSpanOrig - 2*pi) < tol
                    iin = find(abs(abs(obj.x) - pi/2) < tol);
                    xAdd = xshift(obj.x(iin));
                    yAdd = yshift(obj.y(iin));
                    if strcmp(obj.gridType,'PhTh')
                        % Negative sign for y < 0, positive for the rest
                        ySignChange = 0;
                        E1in = bsxfun(@times,obj.E1(iin,:),1-2.*(yAdd < (ySignChange-tol)));
                        E2in = bsxfun(@times,obj.E2(iin,:),1-2.*(yAdd < (ySignChange-tol)));
                    else
                        % Negative sign for |y| > pi/2, positive for the rest
                        ySignChange = pi/2;
                        E1in = bsxfun(@times,obj.E1(iin,:),1-2.*(abs(yAdd) > (ySignChange-tol)));
                        E2in = bsxfun(@times,obj.E2(iin,:),1-2.*(abs(yAdd) > (ySignChange-tol)));
                    end
                    obj = insertMissingCuts(obj,iin,xAdd,yAdd,E1in,E2in);
                end
                if strcmp(obj.gridType,'PhTh')
                    % Remove the pole, and insert the original
                    iout = find(abs(obj.y) < tol);
                    obj = removeDirs(obj,iout);
                else
                    if abs(ySpanOrig - pi) < tol
                        % Copy the pi cut to -pi
                        iy0 = find(abs(obj.y - pi) < tol);
                    else
                        iy0 = [];
                    end
                    [xy0,yy0,E1y0,E2y0] = deal(obj.x(iy0),obj.y(iy0),obj.E1(iy0,:),obj.E2(iy0,:));
                    yy0 = -yy0;
                end
                obj = insertDirs(obj,xy0,yy0,E1y0,E2y0);
                % Sort
                obj = obj.sortGrid;
%                 % Reset coordinate Type
                obj = coorTypeHandle(obj,false);
            end
        end
        
        function obj = range360pos(obj)
            % RANGE360POS sets the range to 'pos', '360'

            assert(any(strcmp(obj.gridType,obj.sphereGrids)),['Can only do range transformations on spherical grids, not on a ',obj.gridType,' grid.'])

            tol = 10^(-obj.nSigDig);
            
            if strcmp(obj.xRangeType,'pos') && strcmp(obj.yRangeType,'360') 
                % Do nothing - already there
                return;
            else% Get original coordinate system
                coorTypeHandle = str2func(['coor2',obj.coorType]);
                obj = obj.rangeChangeCoorTrans;
                
                % Set to standard base
                obj = obj.range180sym;
                % Bookkeeping on the base object
                xSpanOrig = max(obj.x) - min(obj.x);
                ySpanOrig = diff(obj.yRange);
                NangOrig = obj.Nang;
                
                iShift = find(obj.x < -tol);
                xshift = @(x) pi + x;
                ixPoleCuts = find(obj.x >= -tol);
                if strcmp(obj.gridType,'PhTh')
                    yshift = @(y) 2*pi - y;
                    % Grab the required y = pi cut for later insertion to
                    % ensure a smooth pole
                    iy0 = intersect(find(abs(obj.y - pi) < tol),ixPoleCuts);
                else
                    yshift = @(y) pi - (sign(y) + ~sign(y)).*abs(y); % Make the sign(y) = 0 -> 1
                    % Grab the required y = 0 cut for later insertion to
                    % ensure a smooth pole
                    iy0 = intersect(find(abs(obj.y - 0) < tol),ixPoleCuts);
                end
                % Fix E-field signs
                iE = repmat(ismember((1:obj.Nang).',iShift(:)),1,obj.Nf);
                obj.E1(iE) = -obj.E1(iE);
                obj.E2(iE) = -obj.E2(iE);
                [xy0,yy0,E1y0,E2y0] = deal(obj.x(iy0),obj.y(iy0),obj.E1(iy0,:),obj.E2(iy0,:));
                obj.x(iShift) = xshift(obj.x(iShift));
                obj.y(iShift) = yshift(obj.y(iShift));
                obj.y = wrap22pi(obj.y);    % Works for all cases, but not needed for PhTh since already positive from the base...
                if abs(xSpanOrig - 2*pi) < tol
                    iin = find(abs(obj.x) < tol);
                    xAdd = xshift(obj.x(iin));
                    yAdd = wrap22pi(yshift(obj.y(iin)));
                    if strcmp(obj.gridType,'PhTh')
                        % Negative sign for y > pi, positive for the rest
                        ySignChange = pi;
                    else
                        % Negative sign for y > pi/2, positive for the rest
                        ySignChange = pi/2;
                    end
                    E1in = bsxfun(@times,obj.E1(iin,:),1-2.*(yAdd >= (ySignChange-tol)));
                    E2in = bsxfun(@times,obj.E2(iin,:),1-2.*(yAdd >= (ySignChange-tol)));
                    obj = insertMissingCuts(obj,iin,xAdd,yAdd,E1in,E2in);
                end
                % Remove the 180 pole, and insert the original
                if strcmp(obj.gridType,'PhTh')
                    iout = find(abs(obj.y - pi) < tol);
                    obj = removeDirs(obj,iout);
                else
                    yy0 = yy0 + 2*pi; % Not overwriting a pole, but inserting a new one somewhere else
                end
                if abs(ySpanOrig - pi) < tol
                    obj = insertDirs(obj,xy0,yy0,E1y0,E2y0); % Sorted in this function
                end
                % Sort
                obj = obj.sortGrid;
                % Reset coordinate Type
                obj = coorTypeHandle(obj,false);
            end
        end

        function obj = rangeChangeCoorTrans(obj)
            % RANGECHANGECOORTRANS transforms to the appropriate coordinate
            % type for a given gridType for a range change operation
            
            if ~strcmp(obj.coorType,'power')    % If power, do nothing
                switch obj.gridType
                    case 'PhTh'
                        obj = obj.coor2spherical(false);
                    case 'AzEl'
                        obj = obj.coor2Ludwig2AE(false);
                    case 'ElAz'
                        obj = obj.coor2Ludwig2EA(false);
                    otherwise % Astrogrids
                        obj = obj.coor2power; % Not really sure when this will be used, so play safe
                end
            end
        end
        
        %% Other utilities
        function [xRangeType,yRangeType] = setRangeTypes(obj)
            % SETRANGETYPES Returns current rangeType
            
            % Try to figure out what the current rangeType is.
            % Not much error checking is done - assume somewhat
            % sensible inputs are provided most of the time.
            xRangeType = 'sym';
            if any(strcmp(obj.gridType,obj.sphereGrids))
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
        
        function objNew = insertMissingCuts(obj,iin,xAdd,yAdd,E1Add,E2Add,E3Add)
            % INSERTMISSINGCUTS Inserts missing cuts. Used by
            % set{X|Y}range.
            
            objNew = obj;
            if nargin < 7
                if numel(obj.E3) > 0
                    E3Add = obj.E3(iin,:);
                else
                    E3Add = [];
                end
            end
            if nargin < 6
                if numel(obj.E2) > 0
                    E2Add = obj.E2(iin,:);
                else
                    E2Add = [];
                end
            end
            if nargin < 5
                E1Add = obj.E1(iin,:);
            end
            objNew.x = [objNew.x;xAdd];
            objNew.y = [objNew.y;yAdd];
            objNew.E1 = [objNew.E1;E1Add];
            objNew.E2 = [objNew.E2;E2Add];
            objNew.E3 = [objNew.E3;E3Add];
        end
        
        function obj = removeDirs(obj,iout)
            % REMOVEDIRS removes the directions at indices iout.
            
            % Does not change the base...
            obj.x(iout) = [];
            obj.y(iout) = [];
            obj.E1(iout,:) = [];
            if numel(obj.E2) > 0
                obj.E2(iout,:) = [];
            end
            if numel(obj.E3) > 0
                obj.E3(iout,:) = [];
            end
        end
        
        function obj = insertDirs(obj,xin,yin,E1in,E2in,E3in)
            % INSERTDIRS inserts the provided directions in the object
            
            % Does not change the base
            if nargin < 6
                E3in = [];
            end
            if nargin < 5
                E2in = [];
            end
            obj.x = [obj.x;xin];
            obj.y = [obj.y;yin];
            obj.E1 = [obj.E1;E1in];
            if numel(obj.E2) > 0 && numel(E2in) > 0 
                obj.E2 = [obj.E2;E2in];
            end
            if numel(obj.E3) > 0 && numel(E3in) > 0
                obj.E3 = [obj.E3;E3in];
            end
%             obj = obj.sortGrid(obj.nSigDig-1);
        end
        
        function [stepx,stepy] = gridStep(obj)
            obj = obj.roundGrid;
            NxPlot = obj.Nx;
            NyPlot = obj.Ny;
            if NxPlot*NyPlot == obj.Nang
                stepx = diff(obj.xRange)/(NxPlot-1);
                stepy = diff(obj.yRange)/(NyPlot-1);
            else
                N = round(sqrt(obj.Nang)/2)*2+1;
                stepx = diff(obj.xRange)/(N-1);
                stepy = diff(obj.yRange)/(N-1);
            end
        end
    end
    
end