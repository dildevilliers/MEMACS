classdef SphereField < EMField
    %SPHEREFIELD Electromagnetic field sampled on a spherical surface.
    %
    % Canonical representation:
    %
    %   ph      : Np x 1 azimuth angles [deg]
    %   th      : Np x 1 polar angles   [deg]
    %
    %   Eph     : Np x Nf
    %   Eth     : Np x Nf
    %   Er      : Np x Nf, optional
    %
    %   r       : scalar sphere radius [m]
    %   Prad    : 1 x Nf radiated power [W], optional
    %
    % Angles are always stored in degrees.
    %
    % Field samples are stored point-by-point. Structured Ph-Th grid
    % information is detected automatically and stored in gridInfo.
    %
    % No alternative coordinate or polarization representations are
    % stored in the object.
    %
    % Eph, Eth and Er are far-field pattern quantities.
    %
    % The physical electric field at radius r is
    %
    %   E(r) = Epattern * exp(-1j*k*r)/r
    %
    % where k depends on frequency.


    properties (SetAccess = protected)

        % Azimuth angle phi [deg], Np x 1.
        ph (:,1) double = []

        % Polar angle theta [deg], Np x 1.
        th (:,1) double = []

        % Spherical electric-field components [V/m], Np x Nf.
        Eph double = []
        Eth double = []

        % Optional radial electric-field component [V/m], Np x Nf.
        Er double = []

        % Sphere radius [m].
        r (1,1) double = 1

        % Radiated power [W], 1 x Nf.
        %
        % Empty means that radiated power is unknown / not supplied.
        Prad (1,:) double = []

        % Radiation efficiency, 1 x Nf.
        %
        % Dimensionless, normally between 0 and 1.
        % Defaults to unity when not supplied.
        etaRad (1,:) double = 1

    end


    properties (SetAccess = private)

        % Information describing the structured Ph-Th sampling grid.
        %
        % Fields:
        %   isStructured
        %   phVec
        %   thVec
        %   Nph
        %   Nth
        %   gridIndex
        %   hasPhiPeriodicEndpoint
        %   coversFullPhi
        %   coversFullTheta
        %   isFullSphere
        %   hasRedundantPhi
        %   hasNorthPole
        %   hasSouthPole
        %   northPoleIndex
        %   southPoleIndex
        %
        % gridIndex maps stored samples onto an [Nth x Nph] grid:
        %
        %   A(gridIndex) = sampleValues
        %
        % where rows correspond to theta and columns to phi.
        gridInfo struct = struct()

        symmetry struct = struct("YZ","none","XZ","none","XY","none")
        % symmetrySide meaning:
        % +1   stored data lie on positive side
        % -1   stored data lie on negative side
        % 0   no symmetry on this plane
        symmetrySide struct = struct("YZ",0,"XZ",0,"XY",0)
    end


    properties (Dependent, SetAccess = private)

        Np
        hasEr

        % Grid information
        isStructured
        Nph
        Nth
        coversFullPhi
        coversFullTheta
        isFullSphere
        hasPhiPeriodicEndpoint
        hasRedundantPhi
        hasNorthPole
        hasSouthPole

        % Symmetry convenience
        hasSymmetry
        symmetryFactor
        numberOfSymmetries
        hasFullSphereSymmetryCoverage
        representsFullSphere
    end


    methods

        function obj = SphereField(ph, th, Eph, Eth, freqHz, options)
            %SPHEREFIELD Construct a spherical electromagnetic field.
            %
            % Required inputs:
            %
            %   ph      azimuth coordinates [deg]
            %   th      polar coordinates [deg]
            %   Eph     phi-directed electric field
            %   Eth     theta-directed electric field
            %   freqHz  frequency/frequencies [Hz]
            %
            % ph and th may have any shape, but must contain the same
            % number of elements. They are converted internally to column
            % vectors.
            %
            % Eph and Eth must contain Np x Nf samples. For a
            % single-frequency field, arbitrary vector shape is accepted.
            %
            % Name-value options:
            %
            %   Er
            %   r
            %   Prad
            %   Metadata
            %   Provenance

            arguments
                ph double {mustBeFinite} = []
                th double {mustBeFinite} = []
                Eph double = []
                Eth double = []
                freqHz double {mustBeVector, mustBeFinite, mustBePositive} = 1e9  % Dummy

                options.Er double = []

                options.r (1,1) double ...
                    {mustBeFinite, mustBePositive} = 1

                options.Prad double = []
                options.etaRad double = []

                options.Metadata struct = struct()
                options.Provenance struct = struct()

                options.Symmetry struct = struct()
            end

            % ---------------------------------------------------------
            % Default synthetic field
            % ---------------------------------------------------------

            if nargin == 0

                phVec = 0:10:350;
                thVec = 0:5:180;

                [PH,TH] = meshgrid(phVec,thVec);

                ph = PH(:);
                th = TH(:);

                % Simple smooth synthetic pattern (dipole-like).
                Eth = sind(th);
                Eph = zeros(size(Eth));

                freqHz = 1e9;

                options.r = 1;

                options.Metadata = struct( ...
                    'name', "Default SphereField", ...
                    'description', "Synthetic test pattern");

                options.Provenance = struct( ...
                    'source', "SphereField default constructor");
            elseif nargin < 5

                error( ...
                    'SphereField:NotEnoughInputs', ...
                    ['SphereField requires either no inputs, or the five ' ...
                    'canonical inputs ph, th, Eph, Eth and freqHz.']);

            end

            % Initialise superclass.
            obj@EMField( ...
                freqHz, ...
                Metadata = options.Metadata, ...
                Provenance = options.Provenance);


            % ---------------------------------------------------------
            % Coordinates
            % ---------------------------------------------------------

            if numel(ph) ~= numel(th)
                error( ...
                    'SphereField:CoordinateSizeMismatch', ...
                    'ph and th must contain the same number of points.');
            end

            obj.ph = ph(:);
            obj.th = th(:);

            % ---------------------------------------------------------
            % Detect structured grid
            % ---------------------------------------------------------

            obj.gridInfo = obj.detectStructuredGrid();

            % ---------------------------------------------------------
            % Field components
            % ---------------------------------------------------------

            Np_ = numel(obj.ph);

            obj.Eph = obj.normaliseFieldArray( ...
                Eph, Np_, obj.Nf, 'Eph');

            obj.Eth = obj.normaliseFieldArray( ...
                Eth, Np_, obj.Nf, 'Eth');


            if ~isempty(options.Er)

                obj.Er = obj.normaliseFieldArray( ...
                    options.Er, Np_, obj.Nf, 'Er');
            end

            % Symmetry
            obj.symmetry=obj.validateSymmetry(options.Symmetry);
            obj.symmetrySide=obj.validateSymmetryDomain();
            obj.validateSymmetryBoundary();

            % ---------------------------------------------------------
            % Radius
            % ---------------------------------------------------------

            obj.r = options.r;

            % ---------------------------------------------------------
            % Radiated power
            % ---------------------------------------------------------

            if ~isempty(options.Prad)
                if numel(options.Prad) ~= obj.Nf
                    error('SphereField:PradSizeMismatch', ['Prad must contain one value per frequency. ' ...
                        'Expected %d values, received %d.'],obj.Nf, numel(options.Prad));
                end
                obj.Prad = reshape(options.Prad, 1, []);
            end

            if isempty(options.etaRad)
                obj.etaRad = ones(1,obj.Nf);
            else
                if numel(options.etaRad) ~= obj.Nf
                    error('SphereField:EtaRadSizeMismatch',['etaRad must contain one value per frequency. ' ...
                        'Expected %d values, received %d.'],obj.Nf, numel(options.etaRad));
                end
                etaRad = reshape(options.etaRad,1,[]);
                if any(~isfinite(etaRad)) || any(etaRad < 0) || any(etaRad > 1)
                    warning('SphereField:InvalidEtaRad','etaRad must contain finite values between 0 and 1.');
                end
                obj.etaRad = etaRad;
            end

            % % ---------------------------------------------------------
            % % Detect structured grid
            % % ---------------------------------------------------------
            % 
            % obj.gridInfo = obj.detectStructuredGrid();

        end


        %% Dependent properties

        function value = get.Np(obj)
            value = numel(obj.ph);
        end

        function value = get.hasEr(obj)
            value = ~isempty(obj.Er);
        end

        function value = get.isStructured(obj)
            value = obj.gridInfo.isStructured;
        end

        function value = get.Nph(obj)
            if obj.isStructured
                value = obj.gridInfo.Nph;
            else
                value = [];
            end
        end

        function value = get.Nth(obj)
            if obj.isStructured
                value = obj.gridInfo.Nth;
            else
                value = [];
            end
        end

        function value = get.coversFullPhi(obj)
            value = obj.gridInfo.coversFullPhi;
        end

        function value = get.coversFullTheta(obj)
            value = obj.gridInfo.coversFullTheta;
        end

        function value = get.isFullSphere(obj)
            value = obj.gridInfo.isFullSphere;
        end

        function value = get.hasPhiPeriodicEndpoint(obj)
            value = obj.gridInfo.hasPhiPeriodicEndpoint;
        end

        function value = get.hasRedundantPhi(obj)
            value = obj.gridInfo.hasRedundantPhi;
        end

        function value = get.hasNorthPole (obj)
            value = obj.gridInfo.hasNorthPole ;
        end

        function value = get.hasSouthPole  (obj)
            value = obj.gridInfo.hasSouthPole  ;
        end

        function tf = get.hasSymmetry(obj)
            tf = any([obj.symmetry.YZ,obj.symmetry.XZ,obj.symmetry.XY]~="none");
        end

        function n = get.numberOfSymmetries(obj)
            n = sum([obj.symmetry.YZ,obj.symmetry.XZ,obj.symmetry.XY]~="none");
        end

        function f = get.symmetryFactor(obj)
            f = 2^obj.numberOfSymmetries;
        end

        function tf = get.hasFullSphereSymmetryCoverage(obj)
            tf = obj.checkFullSphereSymmetryCoverage();
        end

        function tf = get.representsFullSphere(obj)
            tf = obj.isFullSphere || obj.hasFullSphereSymmetryCoverage;
        end

        %% Actual public methods
        function [Eph, Eth, Er] = getEfield(obj)
            %GETEFIELD Physical electric field at radius r [V/m].
            %
            % Returns the spherical components after applying the far-field
            % propagation factor exp(-1j*k*r)/r.

            FFfact = exp(-1i .* obj.wavenumber .* obj.r) ./ obj.r;

            Eph = obj.Eph .* reshape(FFfact,1,[]);
            Eth = obj.Eth .* reshape(FFfact,1,[]);

            if obj.hasEr
                Er = obj.Er .* reshape(FFfact,1,[]);
            else
                Er = [];
            end
        end

        function W = getW(obj)
            %GETW Time-average radial power density [W/m^2].
            %
            % Uses the physical electric field at radius r.

            W = obj.getWFromPattern(obj.Eph,obj.Eth);
        end

        function U = getU(obj)
            %GETU Radiation intensity [W/sr].

            U = obj.getUFromPattern(obj.Eph,obj.Eth);
        end

        function D = getDirectivity(obj, options)
            %GETDIRECTIVITY Linear directivity.
            %
            % Name-value option:
            %
            %   PowerSource
            %       "auto"       - use stored Prad if available, otherwise
            %                      calculate it from the field
            %       "stored"     - require stored Prad
            %       "integrated" - calculate Prad from the field

            arguments
                obj
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
            end

            switch options.PowerSource
                case "stored"
                    if isempty(obj.Prad), error('SphereField:NoStoredPrad','No stored radiated power is available.'); end
                    P = obj.Prad;
                case "integrated"
                    P = obj.integratePower();
                case "auto"
                    if isempty(obj.Prad), P = obj.integratePower(); else, P = obj.Prad; end
            end

            D = obj.getDirectivityFromPattern(obj.Eph,obj.Eth,P);
        end

        function G = getGain(obj, options)
            %GETGAIN Linear gain.
            %
            %   G = etaRad .* D
            %   PowerSource
            %       "auto"       - use stored Prad if available, otherwise
            %                      calculate it from the field
            %       "stored"     - require stored Prad
            %       "integrated" - calculate Prad from the field

            arguments
                obj
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
            end

            D = obj.getDirectivity(PowerSource=options.PowerSource);
            G = D.*reshape(obj.etaRad,1,[]);
        end

        function [rHat, thHat, phHat] = getCartesianDirections(obj)
            %GETCARTESIANDIRECTIONS Cartesian spherical unit vectors.
            %
            % Outputs are 3 x Np:
            %
            %   rHat  = radial unit vector
            %   thHat = increasing-theta unit vector
            %   phHat = increasing-phi unit vector

            sth = sind(obj.th).';
            cth = cosd(obj.th).';

            sph = sind(obj.ph).';
            cph = cosd(obj.ph).';

            rHat = [sth.*cph; sth.*sph; cth];
            thHat = [cth.*cph; cth.*sph; -sth];
            phHat = [-sph; cph; zeros(1,obj.Np)];

        end

        function ECart = getCartesianE(obj)
            %GETCARTESIANE Electric field in Cartesian coordinates.
            %
            % Output:
            %   ECart : 3 x Np x Nf
            %
            % Applies the canonical outgoing-wave factor exp(-jkr)/r.

            ECart = obj.getCartesianPattern();
            FFfact = reshape(exp(-1i.*obj.wavenumber.*obj.r)./obj.r,[1 1 obj.Nf]);
            ECart = ECart.*FFfact;
        end

        function ECart = getCartesianPattern(obj)
            %GETCARTESIANPATTERN Canonical far-field pattern in Cartesian components.
            %
            % Output:
            %   ECart : 3 x Np x Nf
            %
            % Uses the stored pattern quantities Eth, Eph and Er [V], without the
            % exp(-jkr)/r propagation factor.

            ECart = obj.sphericalPatternToCartesian(obj.ph,obj.th,obj.Eph,obj.Eth,obj.Er);
        end

        function [isConsistent, details] = checkPoleConsistency(obj, options)
            %CHECKPOLECONSISTENCY Check repeated samples at the spherical poles.
            %
            % All samples representing the same pole must correspond to the same
            % Cartesian electric field, although their Eth/Eph components may
            % differ because the spherical basis rotates with phi.
            %
            % Outputs:
            %
            %   isConsistent
            %       True if all existing pole samples are mutually consistent.
            %
            %   details
            %       Structure containing:
            %           northConsistent
            %           southConsistent
            %           northMaxError
            %           southMaxError
            %
            % Max errors are returned separately for each frequency.

            arguments
                obj
                options.RelativeTolerance (1,1) double {mustBeNonnegative} = 1e-10
                options.AbsoluteTolerance (1,1) double {mustBeNonnegative} = 1e-12
            end

            E = obj.getCartesianE();

            details = struct();

            [details.northConsistent, details.northMaxError] = checkPole(obj.gridInfo.northPoleIndex);

            [details.southConsistent, details.southMaxError] = checkPole(obj.gridInfo.southPoleIndex);

            isConsistent = details.northConsistent && details.southConsistent;

            function [tf, maxError] = checkPole(idx)

                % No pole, or only one stored sample, is trivially consistent.
                if numel(idx) <= 1
                    tf = true;
                    maxError = zeros(1,obj.Nf);
                    return
                end

                Epole = E(:,idx,:);

                % Use first stored pole sample as reference.
                Eref = Epole(:,1,:);
                dE = Epole - Eref;
                err = sqrt(sum(abs(dE).^2,1));
                mag = sqrt(sum(abs(Epole).^2,1));
                maxError = reshape(max(err,[],2), 1, []);
                refScale = reshape(max(mag,[],2), 1, []);
                tol = max(options.AbsoluteTolerance, options.RelativeTolerance.*refScale);
                tf = all(maxError <= tol);
            end

        end

        function G = getGridView(obj,viewType,options)
            %GETGRIDVIEW Return structured views of the canonical Ph-Th grid.
            %
            % Views:
            %   "stored"        - exact stored samples
            %   "singlePeriod"  - one complete phi period
            %   "fullSphere"    - one of four equivalent full-sphere layouts
            %   "interpolation" - sphere extended beyond 4*pi for interpolation
            %
            % Full-sphere layouts:
            %
            %   PhiRange      ThetaRange     nominal range
            %   ---------------------------------------------------------
            %   "positive"    "180"          phi 0..360,   th 0..180
            %   "symmetric"   "180"          phi -180..180, th 0..180
            %   "positive"    "360"          phi 0..180,   th -180..180
            %   "symmetric"   "360"          phi -90..90,  th -180..180
            %
            % No interpolation is performed.

            arguments
                obj
                viewType (1,1) string {mustBeMember(viewType,["stored","singlePeriod","fullSphere","interpolation"])} = "stored"
                options.PhiRange (1,1) string {mustBeMember(options.PhiRange,["positive","symmetric"])} = "positive"
                options.ThetaRange (1,1) string {mustBeMember(options.ThetaRange,["180","360"])} = "180"
                options.PaddingSamples (1,1) double {mustBeInteger,mustBeNonnegative} = 2
            end

            switch viewType
                case "stored"
                    G = obj.getStoredGridView();
                case "singlePeriod"
                    G = obj.getSinglePeriodGridView();
                case "fullSphere"
                    G = obj.getFullSphereGridView(options.PhiRange,options.ThetaRange);
                case "interpolation"
                    G = obj.getInterpolationGridView(options.PaddingSamples);
            end

            if G.isStructured
                G.Nth=numel(G.thVec);
                G.Nph=numel(G.phVec);
            else
                G.Nth=[];
                G.Nph=[];
            end
            if G.isStructured
                assert(size(G.Eth,1)==G.Nth && size(G.Eth,2)==G.Nph,'SphereField:InvalidGridView','Grid dimensions are inconsistent with field-array dimensions.');
            end
        end

        function PradInt = integratePower(obj)
            %INTEGRATEPOWER Integrate total radiated power over the represented sphere.
            %
            % The stored field must either explicitly cover the full sphere or form
            % a complete symmetry sector from which the full sphere is represented.
            %
            % Symmetry sectors are integrated over the stored domain and multiplied
            % by symmetryFactor.

            if ~obj.representsFullSphere
                error('SphereField:PowerIntegrationRequiresFullSphere',...
                    'Radiated power integration requires full-sphere coverage, explicitly or through symmetry.');
            end

            hasAzSymmetry = obj.symmetry.YZ ~= "none" || obj.symmetry.XZ ~= "none";

            if obj.isFullSphere || ~hasAzSymmetry
                G = obj.getGridView("singlePeriod");
            else
                G = obj.getGridView("stored");
            end

            Eph_ = reshape(G.Eph,[],obj.Nf);
            Eth_ = reshape(G.Eth,[],obj.Nf);
            U = obj.getUFromPattern(Eph_,Eth_);
            U = reshape(U,G.Nth,G.Nph,obj.Nf);

            phRad = deg2rad(G.phVec);
            thRad = deg2rad(G.thVec);

            sinTh = reshape(sin(thRad),G.Nth,1,1);
            integrand = U.*sinTh;

            intPhi = trapz(phRad,integrand,2);
            PradInt = squeeze(trapz(thRad,intPhi,1)).';

            PradInt = obj.symmetryFactor.*PradInt;
        end

        function [E1,E2,E3,info] = getLinearBasis(obj, options)
            %GETLINEARBASIS Return canonical field in a requested linear basis.
            %
            %   [E1,E2,E3,info] = obj.getLinearBasis(Basis = basis)
            %
            % Available bases:
            %
            %   "spherical"   -> Eth, Eph, optional Er
            %   "Ludwig1"     -> Ex, Ey, Ez
            %   "Ludwig2AE"   -> Eaz, Eel
            %   "Ludwig2EA"   -> Eal, Eep
            %   "Ludwig3"     -> Eh, Ev
            %
            % The returned quantities are far-field pattern quantities, i.e.
            % the canonical exp(-jkr)/r propagation factor is NOT applied.
            %
            % No interpolation or modification of SphereField is performed.

            arguments
                obj
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "spherical"
            end

            [E1,E2,E3,info]=obj.getLinearBasisFromData(obj.ph,obj.th,obj.Eph,obj.Eth,obj.Er,options.Basis);

        end

        function [E1,E2,E3,info] = getField(obj,options)
            %GETFIELD Return field components in requested basis and polarization.
            %
            % Basis:
            %   "spherical"   -> Eth, Eph, optional Er
            %   "Ludwig1"     -> Ex, Ey, Ez
            %   "Ludwig2AE"   -> Eaz, Eel
            %   "Ludwig2EA"   -> Eal, Eep
            %   "Ludwig3"     -> Eh, Ev
            %
            % Polarization:
            %   "linear"   - linear components of requested Basis
            %   "circular" - IEEE RHCP/LHCP components
            %   "slant"    - linear basis rotated by SlantAngle
            %
            % SlantAngle is in degrees and measured from component 1 toward
            % component 2.

            arguments
                obj
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "spherical"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
            end

            [E1lin,E2lin,E3lin,basisInfo]=obj.getLinearBasis(Basis=options.Basis);
            [E1,E2,E3,info]=obj.applyPolarization(E1lin,E2lin,E3lin,basisInfo,options.Polarization,options.SlantAngle);
        end

        function C = getCoordinates(obj,options)
            %GETCOORDINATES Return sample locations in the requested coordinates.
            %
            % No interpolation or modification of SphereField is performed.
            %
            % Coordinates:
            %   "PhTh"     : phi, theta [deg]
            %   "AzEl"     : azimuth, elevation [deg]
            %   "ElAz"     : epsilon, alpha [deg]
            %   "DirCos"   : u, v (with w also returned)
            %   "TrueView" : theta*cos(phi), theta*sin(phi) [deg]
            %   "ArcSin"   : asin(u), asin(v) [deg]

            arguments
                obj
                options.Coordinates (1,1) string {mustBeMember(options.Coordinates,["PhTh","AzEl","ElAz","DirCos","TrueView","ArcSin"])} = "PhTh"
            end

            C = obj.getCoordinatesFromPhTh(obj.ph,obj.th,options.Coordinates);
            C.hasStructuredTopology = obj.isStructured;
        end

        function D = getData(obj,options)
            %GETDATA Extract data in a requested representation.
            %
            % FieldType:
            %   "pattern" - stored far-field pattern quantities rE [V]
            %   "field"   - physical electric field at radius r [V/m]
            %
            % Quantity:
            %   "complex", "magnitude", "phase", "real", "imag"
            %   "W", "U", "directivity", "gain"
            %
            % No interpolation, resampling or modification is performed.

            arguments
                obj
                options.GridView (1,1) string {mustBeMember(options.GridView,["stored","singlePeriod","fullSphere","interpolation"])} = "stored"
                options.PhiRange (1,1) string {mustBeMember(options.PhiRange,["positive","symmetric"])} = "positive"
                options.ThetaRange (1,1) string {mustBeMember(options.ThetaRange,["180","360"])} = "180"
                options.PaddingSamples (1,1) double {mustBeInteger,mustBeNonnegative} = 2
                options.Coordinates (1,1) string {mustBeMember(options.Coordinates,["PhTh","AzEl","ElAz","DirCos","TrueView","ArcSin"])} = "PhTh"
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "spherical"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
                options.FieldType (1,1) string {mustBeMember(options.FieldType,["pattern","field"])} = "pattern"
                options.Quantity (1,1) string {mustBeMember(options.Quantity,["complex","magnitude","phase","real","imag","W","U","directivity","gain"])} = "complex"
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
                options.FrequencyIndex double {mustBeVector,mustBeFinite} = 1
                % options.Component (1,1) string = "all"
                options.ComponentIndex (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.Scale (1,1) string {mustBeMember(options.Scale,["linear","dB10","dB20"])} = "linear"
                options.Normalize (1,1) string {mustBeMember(options.Normalize,["none","max","specified"])} = "none"
                options.ReferenceValue (1,1) double = NaN
            end

            fi = obj.validateFrequencyIndex(options.FrequencyIndex);

            G = obj.getGridView(options.GridView,PhiRange=options.PhiRange,ThetaRange=options.ThetaRange,PaddingSamples=options.PaddingSamples);
            if G.isStructured
                [PH,TH]=meshgrid(G.phVec,G.thVec);
                ph_ = PH(:); th_ = TH(:);
                Eph_ = reshape(G.Eph,[],obj.Nf);
                Eth_ = reshape(G.Eth,[],obj.Nf);
                if isempty(G.Er), Er_ = []; else, Er_ = reshape(G.Er,[],obj.Nf); end
                C = obj.getCoordinatesFromPhTh(ph_,th_,options.Coordinates);
                C.hasStructuredTopology=true;
            else
                ph_ = obj.ph; th_ = obj.th;
                Eph_ = obj.Eph; Eth_ = obj.Eth; Er_ = obj.Er;
                C = obj.getCoordinates(Coordinates=options.Coordinates);
            end

            D = struct();
            D.gridView=options.GridView;
            D.coordinates = options.Coordinates;
            D.x = C.x; D.y = C.y;
            D.xName = C.xName; D.yName = C.yName;
            D.xUnit = C.xUnit; D.yUnit = C.yUnit;
            D.hasStructuredTopology = C.hasStructuredTopology;
            D.frequencyIndex = fi;
            D.freqHz = obj.freqHz(fi);
            D.quantity = options.Quantity;
            D.fieldType = options.FieldType;
            D.basis = options.Basis;
            D.polarization = options.Polarization;
            D.componentNames = string.empty(1,0);
            D.valueUnit = "";
            if any(options.Quantity==["directivity","gain"]), D.powerSource=options.PowerSource; end
            if G.isStructured
                D.Nth = G.Nth;
                D.Nph = G.Nph;
            else
                D.Nth = [];
                D.Nph = [];
            end
            D.isStructured = G.isStructured;
            D.gridView = options.GridView;

            if isfield(C,"z") && ~isempty(C.z)
                D.z = C.z; D.zName = C.zName; D.zUnit = C.zUnit;
            else
                D.z = [];
            end

            fieldQuantities = ["complex","magnitude","phase","real","imag"];

            if any(options.Quantity == fieldQuantities)
                if G.isStructured
                    [E1lin,E2lin,E3lin,info] = obj.getLinearBasisFromData(ph_,th_,Eph_,Eth_,Er_,options.Basis);
                    [E1,E2,E3,info] = obj.applyPolarization(E1lin,E2lin,E3lin,info,options.Polarization,options.SlantAngle);
                else
                    [E1,E2,E3,info] = obj.getField(Basis=options.Basis,Polarization=options.Polarization,SlantAngle=options.SlantAngle);
                end
                E = cat(3,E1(:,fi),E2(:,fi));
                if ~isempty(E3), E = cat(3,E,E3(:,fi)); end

                if options.FieldType == "field"
                    FFfact = exp(-1i*obj.wavenumber(fi)*obj.r)/obj.r;
                    E = E.*reshape(FFfact,1,[],1);
                    D.valueUnit = "V/m";
                else
                    D.valueUnit = "V";
                end

                D.componentNames = info.componentNames;
                if isfield(info,"convention"), D.polarizationConvention = info.convention; end
                if isfield(info,"slantAngle"), D.slantAngle = info.slantAngle; end

                switch options.Quantity
                    case "complex"
                        D.values = E;
                    case "magnitude"
                        D.values = abs(E);
                    case "phase"
                        D.values = angle(E)*180/pi;
                        D.valueUnit = "deg";
                    case "real"
                        D.values = real(E);
                    case "imag"
                        D.values = imag(E);
                end
            else
                P = [];
                if any(options.Quantity == ["directivity","gain"])
                    switch options.PowerSource
                        case "stored"
                            if isempty(obj.Prad), error('SphereField:NoStoredPrad','No stored radiated power is available.'); end
                            P = obj.Prad;
                        case "integrated"
                            P = obj.integratePower();
                        case "auto"
                            if isempty(obj.Prad), P = obj.integratePower(); else, P = obj.Prad; end
                    end
                end
                switch options.Quantity
                    case "W"
                        X=obj.getWFromPattern(Eph_,Eth_);
                        D.values=X(:,fi);
                        D.componentNames="W";
                        D.valueUnit="W/m^2";

                    case "U"
                        X=obj.getUFromPattern(Eph_,Eth_);
                        D.values=X(:,fi);
                        D.componentNames="U";
                        D.valueUnit="W/sr";

                    case "directivity"
                        X=obj.getDirectivityFromPattern(Eph_,Eth_,P);
                        D.values=X(:,fi);
                        D.componentNames="D";
                        D.valueUnit="";

                    case "gain"
                        X=obj.getGainFromPattern(Eph_,Eth_,P);
                        D.values=X(:,fi);
                        D.componentNames="G";
                        D.valueUnit="";
                end
            end

            % Normalization
            switch options.Normalize
                case "none"
                    ref = ones(1,numel(fi));

                case "max"
                    % One common reference over all directions and components,
                    % independently for each selected frequency.
                    if ndims(D.values) > 2
                        ref = max(max(abs(D.values),[],1),[],3);
                    else
                        ref = max(abs(D.values),[],1);
                    end
                    ref = reshape(ref,1,[]);

                    nz = ref ~= 0;
                    if ndims(D.values) > 2
                        D.values(:,nz,:) = D.values(:,nz,:)./reshape(ref(nz),1,[],1);
                    else
                        D.values(:,nz) = D.values(:,nz)./ref(nz);
                    end

                case "specified"
                    ref = options.ReferenceValue;

                    if ~isfinite(real(ref)) || ~isfinite(imag(ref)) || ref == 0
                        error('SphereField:InvalidReferenceValue',...
                            'ReferenceValue must be finite and nonzero.');
                    end

                    if options.Quantity == "complex"
                        D.values = D.values/ref;
                    else
                        D.values = D.values/abs(ref);
                    end
            end

            D.normalization = options.Normalize;
            D.referenceValue = ref;

            % Component selection
            Nc = numel(D.componentNames);

            if options.ComponentIndex > Nc
                compList = strings(1,Nc);
                for ii = 1:Nc
                    compList(ii) = sprintf('%d=%s',ii,D.componentNames(ii));
                end
                error('SphereField:InvalidComponentIndex',...
                    'ComponentIndex %d is invalid. Use 0 for all components, or one of: %s.',...
                    options.ComponentIndex,strjoin(compList,", "));
            end

            if options.ComponentIndex > 0
                if Nc > 1, D.values = D.values(:,:,options.ComponentIndex); end
                D.componentNames = D.componentNames(options.ComponentIndex);
            end

            % Scaling
            switch options.Scale
                case "linear"

                case "dB20"
                    if options.Quantity ~= "magnitude"
                        error('SphereField:InvalidScale',...
                            'dB20 scaling is only valid for magnitude field data.');
                    end
                    D.values = 20*log10(D.values);
                    D.valueUnit = "dB";

                case "dB10"
                    if ~any(options.Quantity == ["W","U","directivity","gain"])
                        error('SphereField:InvalidScale',...
                            'dB10 scaling is only valid for power-like quantities.');
                    end
                    D.values = 10*log10(D.values);
                    D.valueUnit = "dB";
            end

            D.scale = options.Scale;
        end

        function objNew = resample(obj,phNew,thNew,options)
            %RESAMPLE Resample a structured SphereField onto a new Ph-Th tensor grid.
            %
            %   objNew=obj.resample(phNew,thNew)
            %   objNew=obj.resample(phNew,thNew,Method="makima")
            %
            % phNew and thNew are axis vectors defining the new tensor grid.
            % Interpolation is performed on Cartesian far-field pattern components
            % to avoid spherical-basis discontinuities at phi seams and theta poles.
            %
            % Available methods:
            %   "nearest", "linear", "cubic", "spline", "makima"
            %
            % Spherical continuation through phi seams and poles is used where the
            % source-grid topology permits it. No ordinary extrapolation is performed.
            %
            % PaddingSamples controls the number of exact spherical-continuation
            % samples added around supported boundaries. If omitted, a suitable value
            % is selected automatically from the interpolation method.
            %
            % The returned SphereField contains canonical pattern quantities rE and
            % retains r, frequency, Prad, etaRad, metadata and provenance information.

            arguments
                obj
                phNew double {mustBeVector,mustBeReal,mustBeFinite}
                thNew double {mustBeVector,mustBeReal,mustBeFinite}
                options.Method (1,1) string {mustBeMember(options.Method,["nearest","linear","cubic","spline","makima"])} = "linear"
                options.PaddingSamples double = []
            end

            if ~obj.isStructured, error('SphereField:ResampleRequiresStructuredGrid','resample currently requires structured source data.'); end

            phNew=phNew(:).'; thNew=thNew(:).';
            if isempty(phNew) || isempty(thNew), error('SphereField:EmptyResampleGrid','phNew and thNew must be non-empty.'); end
            if numel(unique(phNew))~=numel(phNew) || numel(unique(thNew))~=numel(thNew), error('SphereField:DuplicateResampleCoordinates','phNew and thNew must contain unique values.'); end
            if any(thNew<0 | thNew>180), error('SphereField:InvalidResampleTheta','Canonical resampled theta values must lie in [0,180] deg.'); end

            if isempty(options.PaddingSamples)
                switch options.Method
                    case {"nearest","linear"}, pad=1;
                    case {"cubic","makima"}, pad=2;
                    case "spline", pad=3;
                end
            else
                pad=options.PaddingSamples;
                if ~isscalar(pad) || ~isfinite(pad) || pad<0 || fix(pad)~=pad, error('SphereField:InvalidPaddingSamples','PaddingSamples must be a nonnegative integer.'); end
            end

            G=obj.getGridView("interpolation",PaddingSamples=pad);

            % Full-phi fields are periodic, so map arbitrary phi representations
            % onto the source period before interpolation.
            phEval=phNew;
            if obj.coversFullPhi
                G0=obj.getGridView("singlePeriod");
                ph0=G0.phVec(1);
                phEval=mod(phNew-ph0,360)+ph0;
            end

            tol=1e-10;
            if min(phEval)<min(G.phVec)-tol || max(phEval)>max(G.phVec)+tol
                error('SphereField:ResampleOutsidePhiSupport','Requested phi values extend outside the available interpolation support.');
            end
            if min(thNew)<min(G.thVec)-tol || max(thNew)>max(G.thVec)+tol
                error('SphereField:ResampleOutsideThetaSupport','Requested theta values extend outside the available interpolation support.');
            end

            if options.Method=="cubic"
                if ~obj.isUniformVector(G.phVec) || ~obj.isUniformVector(G.thVec)
                    error('SphereField:CubicRequiresUniformGrid','Cubic interpolation requires uniformly spaced phi and theta grids.');
                end
            end

            % Convert extended spherical pattern components to Cartesian pattern.
            [PH,TH]=meshgrid(G.phVec,G.thVec);
            sth=sind(TH); cth=cosd(TH); sph=sind(PH); cph=cosd(PH);

            Ex=cth.*cph.*G.Eth-sph.*G.Eph;
            Ey=cth.*sph.*G.Eth+cph.*G.Eph;
            Ez=-sth.*G.Eth;

            if ~isempty(G.Er)
                Ex=Ex+sth.*cph.*G.Er;
                Ey=Ey+sth.*sph.*G.Er;
                Ez=Ez+cth.*G.Er;
            end

            % Interpolate Cartesian components.
            ExNew=obj.interpolateComplexGrid(G.thVec,G.phVec,Ex,thNew,phEval,options.Method);
            EyNew=obj.interpolateComplexGrid(G.thVec,G.phVec,Ey,thNew,phEval,options.Method);
            EzNew=obj.interpolateComplexGrid(G.thVec,G.phVec,Ez,thNew,phEval,options.Method);

            % Project back onto the spherical basis at the requested coordinates.
            [PHnew,THnew]=meshgrid(phNew,thNew);
            sth=sind(THnew); cth=cosd(THnew); sph=sind(PHnew); cph=cosd(PHnew);

            EthNew=cth.*cph.*ExNew+cth.*sph.*EyNew-sth.*EzNew;
            EphNew=-sph.*ExNew+cph.*EyNew;
            if obj.hasEr, ErNew=sth.*cph.*ExNew+sth.*sph.*EyNew+cth.*EzNew; else, ErNew=[]; end

            EphNew=reshape(EphNew,[],obj.Nf);
            EthNew=reshape(EthNew,[],obj.Nf);
            if ~isempty(ErNew), ErNew=reshape(ErNew,[],obj.Nf); end

            prov=obj.provenance;
            prov.operation="resample";
            prov.interpolationMethod=options.Method;
            prov.paddingSamples=pad;

            if isempty(ErNew)
                objNew=SphereField(PHnew(:),THnew(:),EphNew,EthNew,obj.freqHz,r=obj.r,Prad=obj.Prad,etaRad=obj.etaRad,Metadata=obj.metadata,Provenance=prov);
            else
                objNew=SphereField(PHnew(:),THnew(:),EphNew,EthNew,obj.freqHz,Er=ErNew,r=obj.r,Prad=obj.Prad,etaRad=obj.etaRad,Metadata=obj.metadata,Provenance=prov);
            end
        end

        function objNew = resampleLike(obj,target,options)
            %RESAMPLELIKE Resample this SphereField onto the structured grid of another SphereField.
            %
            %   objNew=obj.resampleLike(target)
            %   objNew=obj.resampleLike(target,Method=method)
            %   objNew=obj.resampleLike(target,Method=method,PaddingSamples=n)
            %
            % Uses the stored Ph-Th grid of target as the new sampling grid and calls
            % resample internally. The electromagnetic data in target are not used.
            %
            % target must currently contain structured spherical data.
            %
            % Available interpolation methods:
            %   "nearest", "linear", "cubic", "spline", "makima"
            %
            % PaddingSamples controls the number of exact spherical-continuation samples
            % added around supported phi seams and poles before interpolation. If omitted,
            % the default for the selected interpolation method is used.
            %
            % The returned SphereField retains the frequency, radius, Prad, etaRad,
            % metadata and provenance of the source object, not those of target.
            %
            % Example:
            %   SF2r=SF2.resampleLike(SF1,Method="makima");
            %
            % This is useful when two fields must be placed on exactly the same angular
            % sampling grid for comparison or further processing.

            arguments
                obj
                target SphereField
                options.Method (1,1) string {mustBeMember(options.Method,["nearest","linear","cubic","spline","makima"])} = "linear"
                options.PaddingSamples double = []
            end

            if ~target.isStructured, error('SphereField:ResampleLikeRequiresStructuredTarget','resampleLike currently requires a structured target SphereField.'); end

            G=target.getGridView("stored");
            if isempty(options.PaddingSamples)
                objNew=obj.resample(G.phVec,G.thVec,Method=options.Method);
            else
                objNew=obj.resample(G.phVec,G.thVec,Method=options.Method,PaddingSamples=options.PaddingSamples);
            end
        end

        function objNew = expandSymmetry(obj)
            %EXPANDSYMMETRY Explicitly expand all declared field symmetries.
            %
            %   objNew = obj.expandSymmetry()
            %
            % Electric and magnetic symmetry across the YZ, XZ and XY Cartesian
            % planes are expanded by reflecting both the propagation directions and
            % Cartesian field vectors.
            %
            % Expansion is applied to whatever angular samples are stored. Full-sphere
            % coverage is not required. Consequently, the returned SphereField may be
            % structured or unstructured and may or may not cover the full sphere.
            %
            % The source SphereField is unchanged. All symmetry flags in the returned
            % SphereField are cleared because the symmetry has been made explicit.

            if ~obj.hasSymmetry
                objNew = obj;
                return
            end

            ph_ = obj.ph;
            th_ = obj.th;
            Eph_ = obj.Eph;
            Eth_ = obj.Eth;
            Er_ = obj.Er;

            planes = ["YZ","XZ","XY"];

            for ii = 1:numel(planes)
                p = planes(ii);
                type = obj.symmetry.(p);

                if type == "none", continue; end

                [ph2,th2,Eph2,Eth2,Er2] = obj.reflectSymmetryData(ph_,th_,Eph_,Eth_,Er_,p,type);

                ph_ = [ph_; ph2];
                th_ = [th_; th2];
                Eph_ = [Eph_; Eph2];
                Eth_ = [Eth_; Eth2];
                if ~isempty(Er_), Er_ = [Er_; Er2]; end
            end

            [ph_,th_,Eph_,Eth_,Er_] = obj.mergeSymmetryData(ph_,th_,Eph_,Eth_,Er_);

            prov = obj.provenance;
            prov.operation = "expandSymmetry";
            prov.originalSymmetry = obj.symmetry;

            objNew = SphereField(ph_,th_,Eph_,Eth_,obj.freqHz,r=obj.r,Er=Er_,Prad=obj.Prad,...
                etaRad=obj.etaRad,Metadata=obj.metadata,Provenance=prov);
        end

        % Plotting
        function h = plotGrid(obj,viewType,options)
            %PLOTGRID Plot the SphereField sampling grid in the phi-theta plane.
            %
            %   plotGrid(obj)
            %   plotGrid(obj,viewType)
            %   plotGrid(obj,viewType,Name=Value)
            %   plotGrid(obj,G,Name=Value)
            %
            % Displays the angular sampling points associated with a SphereField grid.
            % This is primarily a diagnostic method for visually checking grid topology,
            % periodic extensions, full-sphere representations, and interpolation padding.
            %
            % viewType may be:
            %   "stored"        - exact stored grid
            %   "singlePeriod"  - one complete phi period
            %   "fullSphere"    - one of the equivalent full-sphere representations
            %   "interpolation" - extended grid used to support interpolation across
            %                     phi seams and theta poles where physically possible
            %
            % Alternatively, G may be a grid struct previously returned by getGridView.
            %
            % Name-value options:
            %   PhiRange       = "positive" | "symmetric"
            %   ThetaRange     = "180" | "360"
            %   PaddingSamples = interpolation-view padding
            %   ShowLines      = show tensor-grid lines instead of points
            %   Marker         = marker symbol, e.g. ".", "o", "x", "+"
            %   Color          = MATLAB color specification
            %   MarkerSize     = marker size
            %   Axes           = axes into which the grid is plotted
            %
            % Examples:
            %   obj.plotGrid()
            %   obj.plotGrid("fullSphere",PhiRange="symmetric",ThetaRange="360")
            %
            %   G=obj.getGridView("interpolation",PaddingSamples=3);
            %   obj.plotGrid(G,Marker="x",Color="r")
            %
            %   ax=axes; hold(ax,"on")
            %   obj.plotGrid("stored",Axes=ax,Marker="o",Color="b")
            %   obj.plotGrid("interpolation",Axes=ax,Marker="x",Color="r")
            %
            % No interpolation, resampling, or modification of SphereField is performed.

            arguments
                obj
                viewType = "stored"
                options.PhiRange (1,1) string {mustBeMember(options.PhiRange,["positive","symmetric"])} = "positive"
                options.ThetaRange (1,1) string {mustBeMember(options.ThetaRange,["180","360"])} = "180"
                options.PaddingSamples (1,1) double {mustBeInteger,mustBeNonnegative} = 2
                options.ShowLines (1,1) logical = false
                options.Marker (1,1) string = "."
                options.Color = []
                options.MarkerSize (1,1) double {mustBePositive} = 8
                options.Axes = []
            end

            if isstruct(viewType)
                G=viewType;
            else
                viewType=string(viewType);
                mustBeMember(viewType,["stored","singlePeriod","fullSphere","interpolation"]);
                G=obj.getGridView(viewType,PhiRange=options.PhiRange,ThetaRange=options.ThetaRange,PaddingSamples=options.PaddingSamples);
            end

            [PH,TH]=meshgrid(G.phVec,G.thVec);

            if isempty(options.Axes)
                figure;
                ax=axes;
            else
                ax=options.Axes;
            end

            holdState=ishold(ax);
            hold(ax,"on");

            plotArgs={"Marker",options.Marker,"MarkerSize",options.MarkerSize};
            if ~isempty(options.Color), plotArgs=[plotArgs,{"Color",options.Color}]; end

            if options.ShowLines
                h(1)=plot(ax,PH,TH,"-",plotArgs{:});
                h(2)=plot(ax,PH.',TH.',"-",plotArgs{:});
            else
                h=plot(ax,PH(:),TH(:),"LineStyle","none",plotArgs{:});
            end

            xlabel(ax,"\phi (deg)");
            ylabel(ax,"\theta (deg)");
            grid(ax,"on");
            axis(ax,"tight");

            if ~holdState, hold(ax,"off"); end
        end

        function h = plot(obj,options)
            %PLOT Plot SphereField data on a structured angular grid.
            %
            %   obj.plot()
            %   obj.plot(Name=Value)
            %
            % GridView controls the representation of the spherical sampling grid:
            %   "stored"       - stored grid
            %   "singlePeriod" - one complete phi period
            %   "fullSphere"   - equivalent full-sphere representation
            %
            % For GridView="fullSphere":
            %   PhiRange   = "positive" | "symmetric"
            %   ThetaRange = "180" | "360"
            %
            % Symmetry:
            %   "stored" - plot only the explicitly stored samples (default)
            %   "expand" - explicitly expand all declared electric/magnetic symmetries
            %              before obtaining and plotting the requested data
            %
            % Symmetry expansion does not modify the original SphereField.
            %
            % Examples:
            %   obj.plot(Quantity="magnitude",Component="Eth")
            %   obj.plot(GridView="fullSphere",PhiRange="symmetric",ThetaRange="180",...
            %       Quantity="magnitude",Component="Eth")
            %   obj.plot(GridView="fullSphere",PhiRange="symmetric",ThetaRange="360",...
            %       Quantity="real",Component="Eth")
            %   obj.plot(Quantity="directivity",Scale="dB10")
            %   obj.plot(Component="Eth")
            %   obj.plot(Component="Eth",Symmetry="expand")
            %   obj.plot(Component="Eth",Symmetry="expand",GridView="fullSphere",...
            %       PhiRange="symmetric",ThetaRange="360")

            arguments
                obj
                options.Coordinates (1,1) string {mustBeMember(options.Coordinates,["PhTh","AzEl","ElAz","DirCos","TrueView","ArcSin"])} = "PhTh"
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "spherical"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
                options.FieldType (1,1) string {mustBeMember(options.FieldType,["pattern","field"])} = "pattern"
                options.Quantity (1,1) string {mustBeMember(options.Quantity,["magnitude","phase","real","imag","W","U","directivity","gain"])} = "directivity"
                % options.Component (1,1) string = "all"
                options.ComponentIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.Scale (1,1) string {mustBeMember(options.Scale,["linear","dB10","dB20"])} = "linear"
                options.FrequencyIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.GridView (1,1) string {mustBeMember(options.GridView,["stored","singlePeriod","fullSphere"])} = "stored"
                options.PhiRange (1,1) string {mustBeMember(options.PhiRange,["positive","symmetric"])} = "positive"
                options.ThetaRange (1,1) string {mustBeMember(options.ThetaRange,["180","360"])} = "180"
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
                options.Symmetry (1,1) string {mustBeMember(options.Symmetry,["stored","expand"])} = "stored"
                options.Normalize (1,1) string {mustBeMember(options.Normalize,["none","max","specified"])} = "none"
                options.ReferenceValue (1,1) double = NaN
                options.Axes = []
                options.ShowColorbar (1,1) logical = true
                options.CLim double = []
                options.Title string = ""
            end

            if options.Symmetry == "expand"
                objPlot = obj.expandSymmetry();
            else
                objPlot = obj;
            end

            if ~objPlot.isStructured, error('SphereField:PlotRequiresStructuredGrid','plot currently requires structured data.'); end
            D = objPlot.getData(Coordinates=options.Coordinates,Basis=options.Basis,...
                Polarization=options.Polarization,SlantAngle=options.SlantAngle,...
                FieldType=options.FieldType,Quantity=options.Quantity,...
                ComponentIndex=options.ComponentIndex,Scale=options.Scale,...
                FrequencyIndex=options.FrequencyIndex,GridView=options.GridView,...
                PhiRange=options.PhiRange,ThetaRange=options.ThetaRange,...
                PowerSource=options.PowerSource,Normalize=options.Normalize,...
                ReferenceValue=options.ReferenceValue);
            
            Z = reshape(D.values,D.Nth,D.Nph);
            X = reshape(D.x,D.Nth,D.Nph);
            Y = reshape(D.y,D.Nth,D.Nph);

            if isempty(options.Axes)
                figure;
                ax = axes;
            else
                ax = options.Axes;
            end

            h = surf(ax,X,Y,Z,'EdgeColor','interp','FaceColor','interp');
            view(ax,2);
            axis(ax,'tight');
            xlabel(ax,D.xName);
            ylabel(ax,D.yName);

            if options.ShowColorbar, colorbar(ax); end
            if ~isempty(options.CLim), clim(ax,options.CLim); end

            if options.Title~=""
                title(ax,options.Title);
            else
                title(ax,objPlot.makePlotTitle(D));
            end
        end

        function [h,C] = plotCut(obj,options)
            %PLOTCUT Plot a full great-circle plane cut through a SphereField.
            %
            %   h = obj.plotCut()
            %   h = obj.plotCut(Name=Value)
            %   [h,C] = obj.plotCut(Name=Value)
            %
            % Phi specifies the azimuth of the cut plane. The returned cut angle
            % runs from -180 to 180 degrees:
            %
            %   angle >= 0 : phi = Phi,       theta =  angle
            %   angle <  0 : phi = Phi + 180, theta = -angle
            %
            % The field is interpolated in phi when required.
            %
            % Examples:
            %   obj.plotCut(Phi=0)
            %   obj.plotCut(Phi=45,Basis="Ludwig3",ComponentIndex=2)
            %   obj.plotCut(Phi=0,Quantity="directivity",Scale="dB10")
            %   obj.plotCut(Phi=0,Symmetry="expand")

            arguments
                obj
                options.Phi (1,1) double {mustBeFinite} = 0
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "spherical"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
                options.FieldType (1,1) string {mustBeMember(options.FieldType,["pattern","field"])} = "pattern"
                options.Quantity (1,1) string {mustBeMember(options.Quantity,["magnitude","phase","real","imag","W","U","directivity","gain"])} = "directivity"
                options.ComponentIndex (1,1) double {mustBeInteger,mustBeNonnegative} = 1
                options.Scale (1,1) string {mustBeMember(options.Scale,["linear","dB10","dB20"])} = "linear"
                options.FrequencyIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.Method (1,1) string {mustBeMember(options.Method,["nearest","linear","cubic","spline","makima"])} = "linear"
                options.Symmetry (1,1) string {mustBeMember(options.Symmetry,["stored","expand"])} = "stored"
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
                options.Normalize (1,1) string {mustBeMember(options.Normalize,["none","max","specified"])} = "none"
                options.ReferenceValue (1,1) double = NaN
                options.Axes = []
                options.LineSpec string = ""
                options.Title string = ""
            end

            if options.Symmetry == "expand"
                objPlot = obj.expandSymmetry();
            else
                objPlot = obj;
            end

            C = objPlot.getCutData(options.Phi,Basis=options.Basis,...
                Polarization=options.Polarization,SlantAngle=options.SlantAngle,...
                FieldType=options.FieldType,Quantity=options.Quantity,...
                ComponentIndex=options.ComponentIndex,Scale=options.Scale,...
                Normalize=options.Normalize,ReferenceValue=options.ReferenceValue,...
                FrequencyIndex=options.FrequencyIndex,Method=options.Method,...
                PowerSource=options.PowerSource);

            if isempty(options.Axes)
                figure;
                ax = axes;
            else
                ax = options.Axes;
            end

            Y = reshape(C.values,numel(C.angle),[]);
            if options.LineSpec == ""
                h = plot(ax,C.angle,Y);
            else
                h = plot(ax,C.angle,Y,options.LineSpec);
            end

            grid(ax,'on');
            xlabel(ax,'Cut angle (deg)');
            ylabel(ax,objPlot.makeCutYLabel(C));

            a = C.angle(isfinite(C.angle));
            if numel(a) > 1, xlim(ax,[min(a) max(a)]); end

            if options.Title ~= ""
                title(ax,options.Title);
            else
                title(ax,sprintf('\\phi = %.6g^\\circ cut, %s, %.6g GHz',...
                    C.phi,strjoin(C.componentNames,", "),C.freqHz/1e9));
            end
        end

        function h = plotPrincipalCuts(obj,options)
            %PLOTPRINCIPALCUTS Plot principal great-circle cuts.
            %
            %   obj.plotPrincipalCuts()
            %   obj.plotPrincipalCuts(Name=Value)
            %
            % By default, plots the phi = 0, 45 and 90 degree planes.
            %
            % For field quantities, all available field components are plotted.
            % A single normalization reference is used for all cuts and components.
            %
            % Examples:
            %   obj.plotPrincipalCuts()
            %   obj.plotPrincipalCuts(Basis="Ludwig3")
            %   obj.plotPrincipalCuts(Basis="Ludwig3",Quantity="magnitude",Scale="dB20")
            %   obj.plotPrincipalCuts(Quantity="directivity",Scale="dB10")
            %   obj.plotPrincipalCuts(Phi=[0 90])

            arguments
                obj
                options.Phi (1,:) double {mustBeFinite} = [0 45 90]
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "Ludwig3"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
                options.FieldType (1,1) string {mustBeMember(options.FieldType,["pattern","field"])} = "pattern"
                options.Quantity (1,1) string {mustBeMember(options.Quantity,["magnitude","phase","real","imag","W","U","directivity","gain"])} = "magnitude"
                options.ComponentIndex (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.Scale (1,1) string {mustBeMember(options.Scale,["linear","dB10","dB20"])} = "dB20"
                options.FrequencyIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.Method (1,1) string {mustBeMember(options.Method,["nearest","linear","cubic","spline","makima"])} = "linear"
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
                options.Symmetry (1,1) string {mustBeMember(options.Symmetry,["stored","expand"])} = "stored"
                options.Normalize (1,1) string {mustBeMember(options.Normalize,["none","max","specified"])} = "max"
                options.ReferenceValue (1,1) double = NaN
                options.Axes = []
                options.ShowLegend (1,1) logical = true
                options.Title string = ""
            end

            if options.Symmetry == "expand"
                objPlot = obj.expandSymmetry();
            else
                objPlot = obj;
            end

            if ~objPlot.isStructured
                error('SphereField:PlotPrincipalCutsRequiresStructuredGrid',...
                    'plotPrincipalCuts currently requires structured data.');
            end

            % Resolve the full-field normalization once. All cuts then use exactly
            % the same normalization reference.
            normalize = options.Normalize;
            referenceValue = options.ReferenceValue;

            if normalize == "max"
                Dref = objPlot.getData(Basis=options.Basis,...
                    Polarization=options.Polarization,...
                    SlantAngle=options.SlantAngle,...
                    FieldType=options.FieldType,...
                    Quantity=options.Quantity,...
                    ComponentIndex=0,...
                    Scale="linear",...
                    FrequencyIndex=options.FrequencyIndex,...
                    PowerSource=options.PowerSource,...
                    Normalize="max");

                referenceValue = Dref.referenceValue;
                normalize = "specified";
            end

            Ncut = numel(options.Phi);
            C = cell(1,Ncut);

            for ii = 1:Ncut
                C{ii} = objPlot.getCutData(options.Phi(ii),...
                    Basis=options.Basis,...
                    Polarization=options.Polarization,...
                    SlantAngle=options.SlantAngle,...
                    FieldType=options.FieldType,...
                    Quantity=options.Quantity,...
                    ComponentIndex=options.ComponentIndex,...
                    Scale=options.Scale,...
                    Normalize=normalize,...
                    ReferenceValue=referenceValue,...
                    FrequencyIndex=options.FrequencyIndex,...
                    Method=options.Method,...
                    PowerSource=options.PowerSource);
            end

            if isempty(options.Axes)
                figure;
                ax = axes;
            else
                ax = options.Axes;
            end

            holdState = ishold(ax);
            hold(ax,'on');

            Nc = size(C{1}.values,3);
            h = gobjects(Ncut,Nc);

            % Colour identifies the cut; line style identifies the component.
            lineStyles = ["-","--",":","-."];
            colors = ax.ColorOrder;
            Ncolor = size(colors,1);

            for ii = 1:Ncut
                col = colors(mod(ii - 1,Ncolor) + 1,:);
                Y = reshape(C{ii}.values,numel(C{ii}.angle),[]);

                for jj = 1:Nc
                    ls = lineStyles(mod(jj - 1,numel(lineStyles)) + 1);

                    h(ii,jj) = plot(ax,C{ii}.angle,Y(:,jj),...
                        'Color',col,...
                        'LineStyle',ls,...
                        'DisplayName',sprintf('\\phi = %.6g^\\circ, %s',...
                        C{ii}.phi,C{ii}.componentNames(jj)));

                    % Components of the same cut should have the same colour.
                    if jj > 1
                        h(ii,jj).Color = h(ii,1).Color;
                    end
                end
            end

            if ~holdState, hold(ax,'off'); end

            grid(ax,'on');
            xlabel(ax,'Cut angle (deg)');

            % All cuts should normally have the same angular extent.
            a = C{1}.angle(isfinite(C{1}.angle));
            if numel(a) > 1, xlim(ax,[min(a) max(a)]); end

            if C{1}.valueUnit ~= ""
                ylabel(ax,sprintf('%s (%s)',options.Quantity,C{1}.valueUnit));
            else
                ylabel(ax,options.Quantity);
            end

            if options.ShowLegend
                legend(ax,'show','Location','best');
            end

            if options.Title ~= ""
                title(ax,options.Title);
            else
                title(ax,sprintf('Principal cuts, %.6g GHz',C{1}.freqHz/1e9));
            end
        end

        function h = plot3D(obj,options)
            %PLOT3D Plot SphereField data as a 3-D radiation-pattern surface.
            %
            %   obj.plot3D()
            %   obj.plot3D(Name=Value)
            %
            % The requested field quantity determines the colour of the surface.
            % The radial coordinate is normalized for visualization.
            %
            % For linear data:
            %   R = abs(value)/max(abs(value))
            %
            % For dB data, Floor specifies the radial dynamic range relative to
            % the maximum plotted value. For example, Floor=-40 maps:
            %   peak       -> R = 1
            %   peak-20 dB -> R = 0.5
            %   peak-40 dB -> R = 0
            %
            % Examples:
            %   obj.plot3D()
            %   obj.plot3D(Basis="Ludwig3",ComponentIndex=1)
            %   obj.plot3D(Basis="Ludwig3",ComponentIndex=2,...
            %       Quantity="magnitude",Scale="dB20",Normalize="max")
            %   obj.plot3D(Quantity="directivity",Scale="dB10",Normalize="max")
            %   obj.plot3D(Symmetry="expand")

            arguments
                obj
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "Ludwig3"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
                options.FieldType (1,1) string {mustBeMember(options.FieldType,["pattern","field"])} = "pattern"
                options.Quantity (1,1) string {mustBeMember(options.Quantity,["magnitude","W","U","directivity","gain"])} = "magnitude"
                options.ComponentIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.Scale (1,1) string {mustBeMember(options.Scale,["linear","dB10","dB20"])} = "dB20"
                options.FrequencyIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.GridView (1,1) string {mustBeMember(options.GridView,["stored","singlePeriod","fullSphere"])} = "fullSphere"
                options.PhiRange (1,1) string {mustBeMember(options.PhiRange,["positive","symmetric"])} = "positive"
                options.ThetaRange (1,1) string {mustBeMember(options.ThetaRange,["180","360"])} = "180"
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
                options.Symmetry (1,1) string {mustBeMember(options.Symmetry,["stored","expand"])} = "stored"
                options.Normalize (1,1) string {mustBeMember(options.Normalize,["none","max","specified"])} = "max"
                options.ReferenceValue (1,1) double = NaN
                options.Floor (1,1) double {mustBeFinite} = -40
                options.ShowReferenceSphere (1,1) logical = true
                options.ShowAxes (1,1) logical = true
                options.Axes = []
                options.ShowColorbar (1,1) logical = true
                options.CLim double = []
                options.Title string = ""
            end

            if options.Floor >= 0
                error('SphereField:InvalidPlot3DFloor',...
                    'Floor must be negative and specifies the dB range below the peak.');
            end

            if options.Symmetry == "expand"
                objPlot = obj.expandSymmetry();
            else
                objPlot = obj;
            end

            if ~objPlot.isStructured
                error('SphereField:Plot3DRequiresStructuredGrid',...
                    'plot3D currently requires structured data.');
            end

            D = objPlot.getData(Basis=options.Basis,...
                Polarization=options.Polarization,...
                SlantAngle=options.SlantAngle,...
                FieldType=options.FieldType,...
                Quantity=options.Quantity,...
                ComponentIndex=options.ComponentIndex,...
                Scale=options.Scale,...
                FrequencyIndex=options.FrequencyIndex,...
                GridView=options.GridView,...
                PhiRange=options.PhiRange,...
                ThetaRange=options.ThetaRange,...
                PowerSource=options.PowerSource,...
                Normalize=options.Normalize,...
                ReferenceValue=options.ReferenceValue);

            if ~D.isStructured
                error('SphereField:Plot3DRequiresStructuredGrid',...
                    'The requested grid representation is not structured.');
            end

            % Obtain angular grid corresponding exactly to the requested data view.
            G = objPlot.getGridView(options.GridView,...
                PhiRange=options.PhiRange,...
                ThetaRange=options.ThetaRange);

            [PH,TH] = meshgrid(G.phVec,G.thVec);

            V = reshape(D.values,D.Nth,D.Nph);

            % Convert requested quantity to a normalized plotting radius.
            switch options.Scale
                case "linear"
                    A = abs(V);
                    vmax = max(A,[],'all');

                    if vmax > 0
                        R = A/vmax;
                    else
                        R = zeros(size(A));
                    end

                case {"dB10","dB20"}
                    vmax = max(V,[],'all');
                    vfloor = vmax + options.Floor;

                    Vrad = max(V,vfloor);
                    R = (Vrad - vfloor)/(vmax - vfloor);
            end

            % Spherical pattern surface -> Cartesian plotting coordinates.
            X = R.*sind(TH).*cosd(PH);
            Y = R.*sind(TH).*sind(PH);
            Z = R.*cosd(TH);

            if isempty(options.Axes)
                figure;
                ax = axes;
            else
                ax = options.Axes;
            end

            holdState = ishold(ax);
            hold(ax,'on');

            axis(ax,'equal');
            axis(ax,[-1 1 -1 1 -1 1]*1.15);
            axis(ax,'off');
            view(ax,3);

            a = linspace(0,360,361);

            if options.ShowReferenceSphere
                plot3(ax,cosd(a),sind(a),zeros(size(a)),...
                    'k:','LineWidth',0.5);
                plot3(ax,cosd(a),zeros(size(a)),sind(a),...
                    'k:','LineWidth',0.5);
                plot3(ax,zeros(size(a)),cosd(a),sind(a),...
                    'k:','LineWidth',0.5);
            end

            if options.ShowAxes
                L = 1.1;

                plot3(ax,[-L L],[0 0],[0 0],...
                    'k-','LineWidth',0.5);
                plot3(ax,[0 0],[-L L],[0 0],...
                    'k-','LineWidth',0.5);
                plot3(ax,[0 0],[0 0],[-L L],...
                    'k-','LineWidth',0.5);

                text(ax,L,0,0,'  x');
                text(ax,0,L,0,'  y');
                text(ax,0,0,L,'  z');
            end

            h = surf(ax,X,Y,Z,V,...
                'EdgeColor','none',...
                'FaceColor','interp');

            if ~holdState, hold(ax,'off'); end

            if options.ShowColorbar
                colorbar(ax);
            end

            if ~isempty(options.CLim)
                clim(ax,options.CLim);
            elseif any(options.Scale == ["dB10","dB20"])
                clim(ax,[vmax + options.Floor vmax]);
            end

            if options.Title ~= ""
                title(ax,options.Title);
            else
                title(ax,sprintf('%s, %s, %.6g GHz',...
                    D.quantity,D.componentNames(1),D.freqHz/1e9));
            end
        end
    end


    methods (Access = private)

        function E = normaliseFieldArray(~, E, Np_, Nf_, name)
            %NORMALISEFIELDARRAY Convert field component to Np x Nf.

            if isempty(E)
                error('SphereField:EmptyFieldComponent','%s cannot be empty.', name);
            end

            % Single-frequency case: any vector containing Np samples is
            % acceptable.
            if Nf_ == 1 && isvector(E) && numel(E) == Np_
                E = E(:);
                return
            end

            % General multi-frequency form.
            if ~ismatrix(E) || size(E,1) ~= Np_ || size(E,2) ~= Nf_
                error('SphereField:FieldSizeMismatch', '%s must have size Np x Nf = %d x %d. Received %s.', name, Np_, Nf_, mat2str(size(E)));
            end

        end

        function S = validateSymmetry(~,S)
            S0 = struct("YZ","none","XZ","none","XY","none");
            names = fieldnames(S);
            for ii = 1:numel(names)
                name = names{ii};
                if ~isfield(S0,name)
                    error('SphereField:InvalidSymmetryPlane','Unknown symmetry plane "%s". Valid planes are YZ, XZ and XY.',name);
                end
                value = string(S.(name));
                if ~isscalar(value) || ~any(value == ["none","electric","magnetic"])
                    error('SphereField:InvalidSymmetryType','Symmetry for plane %s must be "none", "electric" or "magnetic".',name);
                end
                S0.(name) = value;
            end
            S = S0;
        end

        function validateSymmetryBoundary(obj)
            if ~obj.hasSymmetry, return; end

            [rHat,~,~] = obj.getCartesianDirections();
            E = obj.getCartesianPattern();

            q = {rHat(1,:).',rHat(2,:).',rHat(3,:).'};
            planes = ["YZ","XZ","XY"];
            normalComponent = [1 2 3];

            geomTol = 1e-10;
            relTol = 1e-8;

            for ii = 1:3
                p = planes(ii);
                type = obj.symmetry.(p);
                if type == "none", continue; end

                idx = abs(q{ii})<=geomTol;
                if ~any(idx), continue; end

                Ep = E(:,idx,:);
                scale = max(abs(Ep),[],'all');
                if scale==0, continue; end

                n = normalComponent(ii);

                if type=="magnetic"
                    err = max(abs(Ep(n,:,:)),[],'all');
                else
                    tangential = setdiff(1:3,n);
                    err = max(abs(Ep(tangential,:,:)),[],'all');
                end

                if err>relTol*scale
                    error('SphereField:SymmetryBoundaryViolation',...
                        'Field does not satisfy the declared %s symmetry on the %s plane.',type,p);
                end
            end
        end

        function side = getSymmetrySide(obj,plane)
            x = sind(obj.th).*cosd(obj.ph);
            y = sind(obj.th).*sind(obj.ph);
            z = cosd(obj.th);

            switch plane
                case "YZ", q = x;
                case "XZ", q = y;
                case "XY", q = z;
            end

            tol = 1e-10;

            hasPos = any(q > tol);
            hasNeg = any(q <- tol);

            if hasPos && hasNeg
                error('SphereField:InvalidSymmetryDomain','Stored samples occur on both sides of the %s symmetry plane.',plane);
            elseif hasPos
                side = 1;
            elseif hasNeg
                side = -1;
            else
                error('SphereField:InvalidSymmetryDomain','All stored samples lie on the %s symmetry plane; the symmetry side cannot be determined.',plane);
            end
        end

        function side = validateSymmetryDomain(obj)
            side = struct("YZ",0,"XZ",0,"XY",0);

            planes = ["YZ","XZ","XY"];
            for ii = 1:numel(planes)
                p = planes(ii);
                if obj.symmetry.(p) ~= "none"
                    side.(p) = obj.getSymmetrySide(p);
                end
            end
        end

        function tf = checkFullSphereSymmetryCoverage(obj)
            %CHECKFULLSPHERESYMMETRYCOVERAGE True if stored symmetry sector is
            % sufficient to reconstruct the complete sphere.

            tf = false;
            if ~obj.hasSymmetry || ~obj.isStructured, return; end

            ph_ = obj.gridInfo.phVec;
            th_ = obj.gridInfo.thVec;
            tol = 1e-8;

            % ----- Theta coverage -----
            if obj.symmetry.XY == "none"
                thRequired = [0 180];
            elseif obj.symmetrySide.XY > 0
                thRequired = [0 90];
            else
                thRequired = [90 180];
            end

            if abs(min(th_) - thRequired(1)) > tol || abs(max(th_) - thRequired(2)) > tol, return; end

            % ----- Phi coverage -----
            hasYZ = obj.symmetry.YZ ~= "none";
            hasXZ = obj.symmetry.XZ ~= "none";
            nAz = hasYZ + hasXZ;

            if nAz == 0
                tf = obj.coversFullPhi;
                return
            end

            % Unwrap so representations such as 270:10:450 are accepted.
            phu = rad2deg(unwrap(deg2rad(ph_)));
            [phMin,iMin] = min(phu);
            [phMax,iMax] = max(phu);
            span = phMax-phMin;

            requiredSpan = 180/nAz;
            if abs(span-requiredSpan)>tol, return; end

            % Sector boundaries must coincide with the active symmetry plane(s).
            p1 = ph_(iMin);
            p2 = ph_(iMax);

            if hasYZ && ~hasXZ
                tf = abs(cosd(p1)) <= tol && abs(cosd(p2)) <= tol;

            elseif hasXZ && ~hasYZ
                tf = abs(sind(p1)) <= tol && abs(sind(p2)) <= tol;

            else
                p1YZ = abs(cosd(p1)) <= tol;
                p1XZ = abs(sind(p1)) <= tol;
                p2YZ = abs(cosd(p2)) <= tol;
                p2XZ = abs(sind(p2)) <= tol;

                tf = (p1YZ && p2XZ) || (p1XZ && p2YZ);
            end
        end

        function info = detectStructuredGrid(obj)
            %DETECTSTRUCTUREDGRID Analyse the canonical Ph-Th sampling grid.
            %
            % The stored samples are never modified or reordered.
            %
            % 0 deg and 360 deg phi samples are retained separately in phVec,
            % but are recognised as periodic representations of the same
            % physical azimuth.

            phVec = unique(obj.ph, 'sorted');
            thVec = unique(obj.th, 'sorted');

            Nph_ = numel(phVec);
            Nth_ = numel(thVec);

            info = struct();

            info.isStructured = false;

            info.phVec = [];
            info.thVec = [];

            info.Nph = [];
            info.Nth = [];

            info.gridIndex = [];

            info.hasPhiPeriodicEndpoint = false;
            info.hasRedundantPhi = false;

            info.coversFullPhi = false;
            info.coversFullTheta = false;
            info.isFullSphere = false;


            % Structured tensor-product grid
            if Nph_*Nth_ ~= obj.Np
                return
            end

            [tfPh, iPh] = ismember(obj.ph, phVec);
            [tfTh, iTh] = ismember(obj.th, thVec);

            if ~all(tfPh) || ~all(tfTh)
                return
            end

            gridIndex = sub2ind([Nth_, Nph_], iTh, iPh);

            if numel(unique(gridIndex)) ~= obj.Np
                return
            end

            info.isStructured = true;

            info.phVec = reshape(phVec,1,[]);
            info.thVec = reshape(thVec,1,[]);

            info.Nph = Nph_;
            info.Nth = Nth_;

            info.gridIndex = gridIndex;

            % Theta coverage
            angleTol = 1e-8;

            info.northPoleIndex = find(abs(obj.th) <= angleTol);
            info.southPoleIndex = find(abs(obj.th - 180) <= angleTol);

            info.hasNorthPole = ~isempty(info.northPoleIndex);
            info.hasSouthPole = ~isempty(info.southPoleIndex);

            info.coversFullTheta = abs(thVec(1)) <= angleTol && abs(thVec(end) - 180) <= angleTol;

            % Phi periodicity and coverage
            if Nph_ >= 2

                phSpan = phVec(end) - phVec(1);

                % Explicit duplicate periodic endpoint, e.g.
                %
                %   0 ... 360
                %  -180 ... 180
                %
                info.hasPhiPeriodicEndpoint = abs(phSpan - 360) <= angleTol;
                info.hasRedundantPhi = phSpan > 360 + angleTol;

                if phSpan >= 360 - angleTol
                    info.coversFullPhi = true;
                else

                    % Interior angular gaps.
                    dph = diff(phVec);

                    % Periodic gap from final sample back to first sample.
                    seamGap = 360 - phSpan;

                    % A complete periodic grid should have a seam gap no larger
                    % than the largest ordinary grid gap.
                    %
                    % Examples:
                    %
                    %   0:10:350   -> seamGap = 10  -> full
                    %   0:10:180   -> seamGap = 180 -> incomplete
                    %
                    maxInteriorGap = max(dph);
                    gapTol = max(angleTol, 1e-8*maxInteriorGap);
                    info.coversFullPhi = seamGap <= maxInteriorGap + gapTol;
                end

            end

            % Complete sphere
            info.isFullSphere = info.coversFullPhi && info.coversFullTheta;

        end

        function G = getStoredGridView(obj)

            G = struct();

            G.viewType = "stored";
            G.isStructured = obj.isStructured;
            G.freqHz = obj.freqHz;

            if obj.isStructured
                G.phVec = obj.gridInfo.phVec;
                G.thVec = obj.gridInfo.thVec;

                G.Eph = obj.fieldToGrid(obj.Eph);
                G.Eth = obj.fieldToGrid(obj.Eth);

                if obj.hasEr
                    G.Er = obj.fieldToGrid(obj.Er);
                else
                    G.Er = [];
                end

                G.hasPhiPeriodicEndpoint = obj.hasPhiPeriodicEndpoint;

            else
                G.ph = obj.ph;
                G.th = obj.th;

                G.Eph = obj.Eph;
                G.Eth = obj.Eth;
                G.Er = obj.Er;

                G.hasPhiPeriodicEndpoint = false;
            end
        end

        function G = getSinglePeriodGridView(obj)

            if ~obj.isStructured
                error('SphereField:SinglePeriodRequiresStructuredGrid', 'singlePeriod is only defined for structured Ph-Th grids.');
            end

            if ~obj.coversFullPhi
                error('SphereField:SinglePeriodRequiresFullPhi', 'singlePeriod requires complete phi coverage.');
            end

            G0 = obj.getStoredGridView();

            phVec = G0.phVec;

            % Use the first stored phi as the start of the period.
            ph0 = phVec(1);
            phEnd = ph0 + 360;

            angleTol = 1e-8;

            % Keep all stored cuts belonging to the first complete revolution.
            keep = phVec >= ph0 - angleTol & phVec <= phEnd + angleTol;

            phOut = phVec(keep);

            Eph_ = G0.Eph(:,keep,:);
            Eth_ = G0.Eth(:,keep,:);

            if obj.hasEr
                Er_ = G0.Er(:,keep,:);
            else
                Er_ = [];
            end

            % If the periodic endpoint is not explicitly present, append it from
            % the first phi cut.
            if abs(phOut(end) - phEnd) > angleTol

                phOut(end+1) = phEnd;

                Eph_(:,end+1,:) = Eph_(:,1,:);
                Eth_(:,end+1,:) = Eth_(:,1,:);

                if obj.hasEr
                    Er_(:,end+1,:) = Er_(:,1,:);
                end

            end

            G = struct();

            G.viewType = "singlePeriod";
            G.isStructured = true;

            G.phVec = phOut;
            G.thVec = G0.thVec;

            G.Eph = Eph_;
            G.Eth = Eth_;
            G.Er = Er_;

            G.freqHz = obj.freqHz;

            G.hasPhiPeriodicEndpoint = true;
        end

        function G = getFullSphereGridView(obj,phiRange,thetaRange)
            if ~obj.isStructured, error('SphereField:FullSphereRequiresStructuredGrid','Full-sphere views require a structured grid.'); end
            if ~obj.isFullSphere, error('SphereField:FullSphereRequiresFullSphere','Full-sphere views require full-sphere data.'); end

            G0 = obj.getSinglePeriodGridView();

            if thetaRange == "180"
                if phiRange == "positive", phStart = 0; else, phStart = -180; end
                G = obj.reindexPhiPeriod(G0,phStart,360);
            else
                if phiRange == "positive", phStart = 0; else, phStart = -90; end
                G = obj.makeFullThetaGrid(G0,phStart);
            end

            G.viewType = "fullSphere";
            G.phiRange = phiRange;
            G.thetaRange = thetaRange;
        end

        function G = getInterpolationGridView(obj,pad)
            %GETINTERPOLATIONGRIDVIEW Return a structured grid extended wherever
            % exact spherical continuation is possible.
            %
            % Possible extensions are treated independently:
            %   phi seam    : requires complete phi coverage
            %   north pole  : requires theta=0 and complete phi coverage
            %   south pole  : requires theta=180 and complete phi coverage
            %
            % No interpolation or extrapolation is performed. All added samples
            % are exact copies of physically equivalent spherical directions.

            if ~obj.isStructured, error('SphereField:InterpolationViewRequiresStructuredGrid','Interpolation views require a structured grid.'); end

            canWrapPhi=obj.coversFullPhi;
            canExtendNorth=obj.hasNorthPole && obj.coversFullPhi;
            canExtendSouth=obj.hasSouthPole && obj.coversFullPhi;

            % For complete phi coverage, first reduce the data to one clean period.
            % For partial phi coverage, retain the stored grid exactly.
            if canWrapPhi
                G=obj.getSinglePeriodGridView();
            else
                G=obj.getStoredGridView();
            end

            G.viewType="interpolation";
            G.paddingSamples=pad;
            G.interpolationSupport=struct("phi",canWrapPhi,"northPole",canExtendNorth,"southPole",canExtendSouth);
            G.extendedPhi=false;
            G.extendedNorthPole=false;
            G.extendedSouthPole=false;

            if pad==0, return; end

            % Pole continuation must be done before extending the phi seam.
            if canExtendNorth || canExtendSouth
                G=obj.extendThetaForInterpolation(G,pad,canExtendNorth,canExtendSouth);
            end
            if canWrapPhi
                G=obj.extendPhiForInterpolation(G,pad);
            end
        end

        function G = extendThetaForInterpolation(obj,G,pad,extendNorth,extendSouth)
            Nth_=numel(G.thVec);
            if pad>Nth_-1, error('SphereField:TooMuchThetaPadding','PaddingSamples is too large for the theta grid.'); end

            % Same physical direction across a pole occurs at phi+180 deg.
            i180=obj.findPhiColumns(G.phVec,G.phVec+180);
            Eth180=G.Eth(:,i180,:);
            Eph180=G.Eph(:,i180,:);
            if isempty(G.Er), Er180=[]; else, Er180=G.Er(:,i180,:); end

            EthLo=[]; EphLo=[]; ErLo=[]; thLo=[];
            EthHi=[]; EphHi=[]; ErHi=[]; thHi=[];

            if extendNorth
                if pad>Nth_-1, error('SphereField:TooMuchNorthPolePadding','Not enough theta samples to extend through the north pole.'); end
                ii=(pad+1):-1:2;
                thLo=-G.thVec(ii);
                EthLo=-Eth180(ii,:,:);
                EphLo=-Eph180(ii,:,:);
                if ~isempty(G.Er), ErLo=Er180(ii,:,:); end
                G.extendedNorthPole=true;
            end

            if extendSouth
                if pad>Nth_-1, error('SphereField:TooMuchSouthPolePadding','Not enough theta samples to extend through the south pole.'); end
                ii=(Nth_-1):-1:(Nth_-pad);
                thHi=360-G.thVec(ii);
                EthHi=-Eth180(ii,:,:);
                EphHi=-Eph180(ii,:,:);
                if ~isempty(G.Er), ErHi=Er180(ii,:,:); end
                G.extendedSouthPole=true;
            end

            G.thVec=[thLo,G.thVec,thHi];
            G.Eth=cat(1,EthLo,G.Eth,EthHi);
            G.Eph=cat(1,EphLo,G.Eph,EphHi);
            if ~isempty(G.Er), G.Er=cat(1,ErLo,G.Er,ErHi); end
        end

        function G = extendPhiForInterpolation(~,G,pad)
            NphUnique=numel(G.phVec)-double(G.hasPhiPeriodicEndpoint);
            if pad>NphUnique-1, error('SphereField:TooMuchPhiPadding','PaddingSamples is too large for the phi grid.'); end

            iLo=(NphUnique-pad+1):NphUnique;
            iHi=2:(pad+1);

            phLo=G.phVec(iLo)-360;
            phHi=G.phVec(iHi)+360;

            G.phVec=[phLo,G.phVec,phHi];
            G.Eth=cat(2,G.Eth(:,iLo,:),G.Eth,G.Eth(:,iHi,:));
            G.Eph=cat(2,G.Eph(:,iLo,:),G.Eph,G.Eph(:,iHi,:));
            if ~isempty(G.Er), G.Er=cat(2,G.Er(:,iLo,:),G.Er,G.Er(:,iHi,:)); end

            G.hasPhiPeriodicEndpoint=false;
            G.extendedPhi=true;
        end

        function G = reindexPhiPeriod(~,G0,phStart,span)
            phSrc = G0.phVec;
            if G0.hasPhiPeriodicEndpoint, phSrc = phSrc(1:end-1); end

            phNew = mod(phSrc-phStart,span)+phStart;
            [phNew,ord] = sort(phNew);

            G = G0;
            G.phVec = phNew;
            G.Eth = G0.Eth(:,ord,:);
            G.Eph = G0.Eph(:,ord,:);
            if ~isempty(G0.Er), G.Er = G0.Er(:,ord,:); end
            G.hasPhiPeriodicEndpoint = false;

            tol = 1e-8;
            if abs(phNew(1)-phStart)<=tol
                G.phVec(end+1) = phStart+span;
                G.Eth(:,end+1,:) = G.Eth(:,1,:);
                G.Eph(:,end+1,:) = G.Eph(:,1,:);
                if ~isempty(G.Er), G.Er(:,end+1,:) = G.Er(:,1,:); end
                G.hasPhiPeriodicEndpoint = true;
            end
        end

        function G = makeFullThetaGrid(obj,G0,phStart)
            phSrc = G0.phVec;
            if G0.hasPhiPeriodicEndpoint, phSrc = phSrc(1:end-1); end

            phHalf = mod(phSrc-phStart,180)+phStart;
            phHalf = unique(phHalf,'sorted');

            i0 = obj.findPhiColumns(G0.phVec,phHalf);
            i180 = obj.findPhiColumns(G0.phVec,phHalf+180);

            th_ = G0.thVec;
            iNeg = numel(th_):-1:2;
            thNew = [-th_(iNeg),th_];

            G = G0;
            G.phVec = phHalf;
            G.thVec = thNew;

            G.Eth = cat(1,-G0.Eth(iNeg,i180,:),G0.Eth(:,i0,:));
            G.Eph = cat(1,-G0.Eph(iNeg,i180,:),G0.Eph(:,i0,:));
            if isempty(G0.Er), G.Er = []; else, G.Er = cat(1,G0.Er(iNeg,i180,:),G0.Er(:,i0,:)); end

            G.hasPhiPeriodicEndpoint = false;
            if abs(phHalf(1)-phStart)<=1e-8
                G.phVec(end+1) = phStart+180;
                G.Eth(:,end+1,:) = G.Eth(:,1,:);
                G.Eph(:,end+1,:) = G.Eph(:,1,:);
                if ~isempty(G.Er), G.Er(:,end+1,:) = G.Er(:,1,:); end
                G.hasPhiPeriodicEndpoint = true;
            end
        end

        function idx = findPhiColumns(~,phVec,phTarget)
            targetSize = size(phTarget);
            phVec = phVec(:);
            phTarget = phTarget(:);
            idx = zeros(size(phTarget));
            tol = 1e-8;

            for ii = 1:numel(phTarget)
                d = mod(phVec-phTarget(ii)+180,360)-180;
                jj = find(abs(d)<=tol,1);
                if isempty(jj), error('SphereField:PhiSamplingNotCompatible','Required phi sample %.12g deg is not present in the source grid.',phTarget(ii)); end
                idx(ii) = jj;
            end
            idx = reshape(idx,targetSize);
        end

        function [E1,E2,E3,info] = getLinearBasisFromData(~,ph,th,Eph,Eth,Er,basis)
            
            E3 = [];
            info = struct();
            info.basis = basis;

            switch basis
                case "spherical"
                    E1=Eth; E2=Eph; E3=Er;
                    if isempty(Er), info.componentNames=["Eth","Eph"]; else, info.componentNames=["Eth","Eph","Er"]; end
                    info.isTransverse=isempty(Er);

                case "Ludwig1"
                    cth=cosd(th); sth=sind(th);
                    cph=cosd(ph); sph=sind(ph);
                    if isempty(Er), Er_=0; else, Er_=Er; end

                    E1=cth.*cph.*Eth-sph.*Eph+sth.*cph.*Er_;
                    E2=cth.*sph.*Eth+cph.*Eph+sth.*sph.*Er_;
                    E3=-sth.*Eth+cth.*Er_;

                    info.componentNames = ["Ex","Ey","Ez"];
                    info.isTransverse = false;

                case "Ludwig2AE"
                    cth = cosd(th);
                    sth = sind(th);

                    cph = cosd(ph);
                    sph = sind(ph);

                    cosEl = sqrt(max( 0, 1 - sth.^2.*sph.^2));

                    singular = cosEl < 1e-12;

                    % Regular points
                    E1 = zeros(size(Eth),'like',Eth);
                    E2 = zeros(size(Eth),'like',Eth);

                    regular = ~singular;

                    E1(regular,:) = (cph(regular).*Eth(regular,:) - cth(regular).*sph(regular).*Eph(regular,:))./cosEl(regular);
                    E2(regular,:) = (cth(regular).*sph(regular).*Eth(regular,:) + cph(regular).*Eph(regular,:))./cosEl(regular);

                    % Coordinate singularities:
                    % theta = 90 deg, phi = +/-90 deg.
                    %
                    % Use the constant-phi limiting value.
                    if any(singular)
                        s = sth(singular).*sph(singular);
                        E1(singular,:) = s.*Eph(singular,:);
                        E2(singular,:) = -s.*Eth(singular,:);
                    end

                    info.componentNames = ["Eaz","Eel"];
                    info.isTransverse = true;

                case "Ludwig2EA"
                    cth = cosd(th);
                    sth = sind(th);

                    cph = cosd(ph);
                    sph = sind(ph);

                    cosAl = sqrt(max(0,1 - sth.^2.*cph.^2));

                    singular = cosAl < 1e-12;

                    E1 = zeros(size(Eth),'like',Eth);
                    E2 = zeros(size(Eth),'like',Eth);

                    regular = ~singular;

                    E1(regular,:) = (cth(regular).*cph(regular).*Eth(regular,:) - sph(regular).*Eph(regular,:))./cosAl(regular);
                    E2(regular,:) = (sph(regular).*Eth(regular,:) + cth(regular).*cph(regular).*Eph(regular,:))./cosAl(regular);

                    % Coordinate singularities:
                    % theta = 90 deg, phi = 0/180/360 deg.
                    %
                    % Use the constant-phi limiting value.
                    if any(singular)
                        s = sth(singular).*cph(singular);
                        E1(singular,:) = -s.*Eth(singular,:);
                        E2(singular,:) = -s.*Eph(singular,:);
                    end

                    info.componentNames = ["Eal","Eep"];
                    info.isTransverse = true;

                case "Ludwig3"
                    cph = cosd(ph);
                    sph = sind(ph);

                    E1 = cph.*Eth - sph.*Eph;
                    E2 = sph.*Eth + cph.*Eph;

                    info.componentNames = ["Eh","Ev"];
                    info.isTransverse = true;
            end
        end

        function [E1,E2,E3,info] = applyPolarization(~,E1lin,E2lin,E3lin,basisInfo,polarization,slantAngle)
            E3=[];
            info=basisInfo;
            info.polarization=polarization;

            switch polarization
                case "linear"
                    E1=E1lin; E2=E2lin; E3=E3lin;

                case "circular"
                    if basisInfo.basis=="Ludwig1", error('SphereField:PolarizationRequiresTransverseBasis','Circular polarization requires a two-component transverse basis.'); end
                    E1=(E1lin+1i*E2lin)/sqrt(2);
                    E2=(E1lin-1i*E2lin)/sqrt(2);
                    info.componentNames=["RHCP","LHCP"];
                    info.convention="IEEE";

                case "slant"
                    if basisInfo.basis=="Ludwig1", error('SphereField:PolarizationRequiresTransverseBasis','Slant polarization requires a two-component transverse basis.'); end
                    E1=cosd(slantAngle).*E1lin+sind(slantAngle).*E2lin;
                    E2=-sind(slantAngle).*E1lin+cosd(slantAngle).*E2lin;
                    info.componentNames=["Eslant","Eorth"];
                    info.slantAngle=slantAngle;
            end
        end

        function C = getCoordinatesFromPhTh(~,ph,th,coordinates)
            u = sind(th).*cosd(ph);
            v = sind(th).*sind(ph);
            w = cosd(th);

            C = struct();
            C.coordinates = coordinates;
            C.z = [];

            switch coordinates
                case "PhTh"
                    C.x = ph; C.y = th;
                    C.xName = "phi"; C.yName = "theta";
                    C.xUnit = "deg"; C.yUnit = "deg";

                case "DirCos"
                    C.x = u; C.y = v; C.z = w;
                    C.xName = "u"; C.yName = "v"; C.zName = "w";
                    C.xUnit = ""; C.yUnit = ""; C.zUnit = "";

                case "AzEl"
                    C.x = atan2d(u,w); C.y = asind(max(-1,min(1,v)));
                    C.xName = "az"; C.yName = "el";
                    C.xUnit = "deg"; C.yUnit = "deg";

                case "ElAz"
                    C.x = atan2d(v,w); C.y = asind(max(-1,min(1,u)));
                    C.xName = "epsilon"; C.yName = "alpha";
                    C.xUnit = "deg"; C.yUnit = "deg";

                case "TrueView"
                    C.x = th.*cosd(ph); C.y = th.*sind(ph);
                    C.xName = "Xg"; C.yName = "Yg";
                    C.xUnit = "deg"; C.yUnit = "deg";

                case "ArcSin"
                    C.x = asind(max(-1,min(1,u))); C.y = asind(max(-1,min(1,v)));
                    C.xName = "asin(u)"; C.yName = "asin(v)";
                    C.xUnit = "deg"; C.yUnit = "deg";
            end
        end

        function ECart = sphericalPatternToCartesian(~,ph,th,Eph,Eth,Er)
            %SPHERICALPATTERNTOCARTESIAN Convert spherical pattern components to Cartesian.
            %
            % Inputs:
            %   ph, th : Np x 1 angular coordinates [deg]
            %   Eph    : Np x Nf phi-directed pattern component
            %   Eth    : Np x Nf theta-directed pattern component
            %   Er     : Np x Nf radial pattern component, or []
            %
            % Output:
            %   ECart  : 3 x Np x Nf Cartesian pattern field

            Np_ = numel(ph);
            Nf = size(Eth,2);

            sth = sind(th).'; cth = cosd(th).';
            sph = sind(ph).'; cph = cosd(ph).';

            rHat = [sth.*cph; sth.*sph; cth];
            thHat = [cth.*cph; cth.*sph; -sth];
            phHat = [-sph; cph; zeros(1,Np_)];

            ECart = reshape(thHat,[3 Np_ 1]).*reshape(Eth,[1 Np_ Nf]) + reshape(phHat,[3 Np_ 1]).*reshape(Eph,[1 Np_ Nf]);

            if ~isempty(Er)
                ECart = ECart + reshape(rHat,[3 Np_ 1]).*reshape(Er,[1 Np_ Nf]);
            end
        end

        function [Eph,Eth,Er] = cartesianPatternToSpherical(~,ph,th,ECart,hasEr)
            %CARTESIANPATTERNTOSPHERICAL Convert Cartesian pattern field to spherical components.

            Np_ = numel(ph);
            Nf = size(ECart,3);

            sth = sind(th).'; cth = cosd(th).';
            sph = sind(ph).'; cph = cosd(ph).';

            rHat = [sth.*cph; sth.*sph; cth];
            thHat = [cth.*cph; cth.*sph; -sth];
            phHat = [-sph; cph; zeros(1,Np_)];

            Eth = reshape(sum(ECart.*reshape(thHat,[3 Np_ 1]),1),Np_,Nf);
            Eph = reshape(sum(ECart.*reshape(phHat,[3 Np_ 1]),1),Np_,Nf);

            if hasEr
                Er = reshape(sum(ECart.*reshape(rHat,[3 Np_ 1]),1),Np_,Nf);
            else
                Er = [];
            end
        end

        function [ph2,th2,Eph2,Eth2,Er2] = reflectSymmetryData(obj,ph,th,Eph,Eth,Er,plane,type)
            %REFLECTSYMMETRYDATA Reflect field samples across a symmetry plane.
            %
            % Magnetic symmetry:
            %   E2 = R*E
            %
            % Electric symmetry:
            %   E2 = -R*E

            switch plane
                case "YZ"
                    R = diag([-1 1 1]);
                    ph2 = mod(180 - ph,360);
                    th2 = th;

                case "XZ"
                    R = diag([1 -1 1]);
                    ph2 = mod(-ph,360);
                    th2 = th;

                case "XY"
                    R = diag([1 1 -1]);
                    ph2 = mod(ph,360);
                    th2 = 180 - th;

                otherwise
                    error('SphereField:InvalidSymmetryPlane','Unknown symmetry plane "%s".',plane);
            end

            switch type
                case "magnetic"
                    s = 1;
                case "electric"
                    s = -1;
                otherwise
                    error('SphereField:InvalidSymmetryType','Symmetry type must be "electric" or "magnetic".');
            end

            ECart = obj.sphericalPatternToCartesian(ph,th,Eph,Eth,Er);

            for ff = 1:size(ECart,3)
                ECart(:,:,ff) = s*R*ECart(:,:,ff);
            end

            [Eph2,Eth2,Er2] = obj.cartesianPatternToSpherical(ph2,th2,ECart,~isempty(Er));
        end

        function [ph,th,Eph,Eth,Er] = mergeSymmetryData(~,ph,th,Eph,Eth,Er)
            %MERGESYMMETRYDATA Remove duplicate angular samples after symmetry expansion.
            %
            % Samples are considered duplicates only when both phi and theta are
            % equal within tolerance. Geometrically equivalent pole samples with
            % different phi values are retained.

            tol = 1e-9;

            ph = mod(ph,360);
            ph(abs(ph) < tol | abs(ph - 360) < tol) = 0;
            th(abs(th) < tol) = 0;
            th(abs(th - 180) < tol) = 180;

            key = round([ph th]/tol)*tol;
            [~,ia] = unique(key,'rows','stable');

            ph = ph(ia);
            th = th(ia);
            Eph = Eph(ia,:);
            Eth = Eth(ia,:);
            if ~isempty(Er), Er = Er(ia,:); end

            % Keep a predictable theta-major ordering.
            [~,idx] = sortrows([th ph],[1 2]);
            ph = ph(idx);
            th = th(idx);
            Eph = Eph(idx,:);
            Eth = Eth(idx,:);
            if ~isempty(Er), Er = Er(idx,:); end
        end
        
        function W = getWFromPattern(obj,Eph,Eth)
            %GETWFROMPATTERN Calculate physical radial power density from pattern data.
            %
            % Eph and Eth are canonical far-field pattern quantities rE [V], with
            % size Np x Nf. W is the physical time-average radial power density at
            % obj.r [W/m^2].

            FFfact = exp(-1i*obj.wavenumber*obj.r)/obj.r;
            EphF = Eph.*reshape(FFfact,1,[]);
            EthF = Eth.*reshape(FFfact,1,[]);
            W = (abs(EthF).^2+abs(EphF).^2)/(2*obj.eta0);
        end

        function U = getUFromPattern(obj,Eph,Eth)
            U = obj.r^2.*obj.getWFromPattern(Eph,Eth);
        end

        function D = getDirectivityFromPattern(obj,Eph,Eth,P)
            U = obj.getUFromPattern(Eph,Eth);
            D = 4*pi.*U./reshape(P,1,[]);
        end

        function G = getGainFromPattern(obj,Eph,Eth,P)
            D = obj.getDirectivityFromPattern(Eph,Eth,P);
            G = D.*reshape(obj.etaRad,1,[]);
        end

        function C = getCutData(obj,phi,options)
            %GETCUTDATA Extract data along a great-circle phi-plane cut.
            %
            % The signed cut angle is defined by:
            %   angle >= 0 : phi = phi0,       theta =  angle
            %   angle <  0 : phi = phi0 + 180, theta = -angle

            arguments
                obj
                phi (1,1) double {mustBeFinite}
                options.Basis (1,1) string {mustBeMember(options.Basis,["spherical","Ludwig1","Ludwig2AE","Ludwig2EA","Ludwig3"])} = "spherical"
                options.Polarization (1,1) string {mustBeMember(options.Polarization,["linear","circular","slant"])} = "linear"
                options.SlantAngle (1,1) double {mustBeFinite} = 45
                options.FieldType (1,1) string {mustBeMember(options.FieldType,["pattern","field"])} = "pattern"
                options.Quantity (1,1) string {mustBeMember(options.Quantity,["complex","magnitude","phase","real","imag","W","U","directivity","gain"])} = "magnitude"
                options.ComponentIndex (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.Scale (1,1) string {mustBeMember(options.Scale,["linear","dB10","dB20"])} = "linear"
                options.Normalize (1,1) string {mustBeMember(options.Normalize,["none","max","specified"])} = "none"
                options.ReferenceValue (1,1) double = NaN
                options.FrequencyIndex (1,1) double {mustBeInteger,mustBePositive} = 1
                options.Method (1,1) string {mustBeMember(options.Method,["nearest","linear","cubic","spline","makima"])} = "linear"
                options.PowerSource (1,1) string {mustBeMember(options.PowerSource,["auto","stored","integrated"])} = "auto"
            end

            if ~obj.isStructured
                error('SphereField:CutRequiresStructuredGrid','Cut extraction currently requires structured data.');
            end

            % Resolve full-field normalization before reducing the data to a cut.
            normalize = options.Normalize;
            referenceValue = options.ReferenceValue;

            if normalize == "max"
                Dref = obj.getData(Basis=options.Basis,Polarization=options.Polarization,...
                    SlantAngle=options.SlantAngle,FieldType=options.FieldType,...
                    Quantity=options.Quantity,ComponentIndex=0,...
                    FrequencyIndex=options.FrequencyIndex,PowerSource=options.PowerSource,...
                    Normalize="max",Scale="linear");

                referenceValue = Dref.referenceValue;
                normalize = "specified";
            end

            G = obj.getGridView("stored");
            thCut = G.thVec;

            ph0 = mod(phi,360);
            phPair = [ph0 mod(ph0 + 180,360)];

            SFcut = obj.resample(phPair,thCut,Method=options.Method);

            D = SFcut.getData(Basis=options.Basis,Polarization=options.Polarization,...
                SlantAngle=options.SlantAngle,FieldType=options.FieldType,...
                Quantity=options.Quantity,ComponentIndex=options.ComponentIndex,...
                Scale=options.Scale,FrequencyIndex=options.FrequencyIndex,...
                PowerSource=options.PowerSource,Normalize=normalize,...
                ReferenceValue=referenceValue);

            Z = reshape(D.values,D.Nth,D.Nph,[]);

            angle = [-fliplr(thCut(2:end)) thCut].';
            values = cat(1,flip(Z(2:end,2,:),1),Z(:,1,:));

            C = struct();
            C.angle = angle;
            C.values = values;
            C.phi = ph0;
            C.componentNames = D.componentNames;
            C.quantity = D.quantity;
            C.fieldType = D.fieldType;
            C.valueUnit = D.valueUnit;
            C.freqHz = D.freqHz;
            C.frequencyIndex = D.frequencyIndex;
            C.scale = D.scale;
            C.normalization = options.Normalize;
            C.referenceValue = D.referenceValue;
            C.basis = D.basis;
            C.polarization = D.polarization;
        end

        function A = fieldToGrid(obj, E)
            %FIELDTOGRID Reshape Np x Nf field samples to Nth x Nph x Nf.

            A = zeros( ...
                obj.Nth, ...
                obj.Nph, ...
                obj.Nf, ...
                'like', E);

            for ff = 1:obj.Nf
                A_ = zeros(obj.Nth,obj.Nph,'like',E);
                A_(obj.gridInfo.gridIndex) = E(:,ff);
                A(:,:,ff) = A_;
            end
        end

        function Vnew = interpolateComplexGrid(~,th,ph,V,thNew,phNew,method)
            Nf=size(V,3);
            Vnew=zeros(numel(thNew),numel(phNew),Nf,'like',V);

            for ff=1:Nf
                Fr=griddedInterpolant({th,ph},real(V(:,:,ff)),char(method),'none');
                Fi=griddedInterpolant({th,ph},imag(V(:,:,ff)),char(method),'none');
                Vnew(:,:,ff)=Fr({thNew,phNew})+1i*Fi({thNew,phNew});
            end
        end

        function tf = isUniformVector(~,x)
            if numel(x)<3, tf=true; return; end
            dx=diff(x);
            tol=1e-10*max(1,max(abs(dx)));
            tf=max(abs(dx-dx(1)))<=tol;
        end

        function str = makePlotTitle(~,D)
            if isscalar(D.componentNames)
                component=D.componentNames;
            else
                component="";
            end

            if component~=""
                str=sprintf('%s, %s, %.6g GHz',component,D.quantity,D.freqHz/1e9);
            else
                str=sprintf('%s, %.6g GHz',D.quantity,D.freqHz/1e9);
            end
        end

        function str = makeCutYLabel(~,D)
            name = D.componentNames;

            if D.valueUnit == ""
                str = sprintf('%s, %s',name,D.quantity);
            else
                str = sprintf('%s, %s (%s)',name,D.quantity,D.valueUnit);
            end
        end
    end

    methods (Static)

        function SF = readCSTffs(pathName,options)
            %READCSTFFS Create a SphereField from a CST .ffs file.
            %
            %   SF = SphereField.readCSTffs(pathName)
            %   SF = SphereField.readCSTffs(pathName,Name=Value)
            %
            % Reads CST far-field source (.ffs) files containing one or more
            % frequencies. The stored field components are interpreted as the
            % spherical far-field pattern quantities rE:
            %
            %   Eth [V]
            %   Eph [V]
            %
            % with angles stored in degrees.
            %
            % CST radiated and accepted powers are used to populate Prad and
            % radiation efficiency.
            %
            % If pathName is empty, a file-selection dialog is opened.

            arguments
                pathName (1,1) string = ""
                options.r (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1
            end

            if pathName == ""
                [name,path] = uigetfile('*.ffs');
                if isequal(name,0)
                    SF = SphereField.empty;
                    return
                end
                pathName = fullfile(path,name);
            end

            [~,~,ext] = fileparts(pathName);
            if ext == ""
                pathName = pathName + ".ffs";
            elseif ~strcmpi(ext,'.ffs')
                error('SphereField:InvalidCSTffsExtension',...
                    'Expected a CST .ffs file.');
            end

            fid = fopen(pathName,'r');
            if fid == -1
                error('SphereField:UnableToOpenFile',...
                    'Unable to open CST far-field file %s.',pathName);
            end
            cleanup = onCleanup(@() fclose(fid));

            % CST markers
            freqMarker = '// #Frequencies';
            posMarker = '// Position';
            zAxisMarker = '// zAxis';
            xAxisMarker = '// xAxis';
            powerFreqMarker = '// Radiated/Accepted/Stimulated Power , Frequency';
            NphNthMarker = '// >> Total #phi samples, total #theta samples';
            fieldMarker = '// >> Phi, Theta, Re(E_Theta), Im(E_Theta), Re(E_Phi), Im(E_Phi):';

            Nf = [];
            Nph = [];
            Nth = [];
            freq = [];
            Prad = [];
            Pacc = [];
            Pstim = [];
            pos = [];
            zAxis = [];
            xAxis = [];

            fCount = 0;

            while ~feof(fid)
                a = fgetl(fid);
                if ~ischar(a), break; end

                if strcmp(a,freqMarker)
                    Nf = fscanf(fid,'%i',1);

                elseif strcmp(a,posMarker)
                    pos = fscanf(fid,'%f%f%f',[3,1]);

                elseif strcmp(a,zAxisMarker)
                    zAxis = fscanf(fid,'%f%f%f',[3,1]);

                elseif strcmp(a,xAxisMarker)
                    xAxis = fscanf(fid,'%f%f%f',[3,1]);

                elseif strncmp(a,powerFreqMarker,length(powerFreqMarker))
                    if isempty(Nf)
                        error('SphereField:InvalidCSTffs',...
                            'Frequency count must appear before the power/frequency block.');
                    end

                    PF = fscanf(fid,'%f',[4,Nf]);

                    Prad = PF(1,:);
                    Pacc = PF(2,:);
                    Pstim = PF(3,:);
                    freq = PF(4,:);

                elseif strncmp(a,NphNthMarker,length(NphNthMarker))
                    NphNth = fscanf(fid,'%i %i',[2,1]);
                    Nph = NphNth(1);
                    Nth = NphNth(2);

                    if isempty(Nf)
                        error('SphereField:InvalidCSTffs',...
                            'Frequency count must appear before the angular-grid definition.');
                    end

                    if fCount == 0
                        Np = Nph*Nth;
                        ph = zeros(Np,1);
                        th = zeros(Np,1);
                        Eth = zeros(Np,Nf);
                        Eph = zeros(Np,Nf);
                    end

                elseif strncmp(a,fieldMarker,length(fieldMarker))
                    if isempty(Nph) || isempty(Nth)
                        error('SphereField:InvalidCSTffs',...
                            'Angular-grid dimensions are missing before the field data.');
                    end

                    fCount = fCount + 1;

                    fData = fscanf(fid,'%f',[6,Nph*Nth]).';

                    if size(fData,1) ~= Nph*Nth
                        error('SphereField:InvalidCSTffs',...
                            'Unexpected number of field samples in frequency block %d.',fCount);
                    end

                    % Angular sampling should be identical at all frequencies.
                    if fCount == 1
                        ph = fData(:,1);
                        th = fData(:,2);
                    else
                        tol = 1e-10;
                        if any(abs(ph - fData(:,1)) > tol) || any(abs(th - fData(:,2)) > tol)
                            error('SphereField:CSTffsFrequencyGridMismatch',...
                                'Angular sampling differs between CST frequency blocks.');
                        end
                    end

                    Eth(:,fCount) = complex(fData(:,3),fData(:,4));
                    Eph(:,fCount) = complex(fData(:,5),fData(:,6));

                    if fCount == Nf
                        break
                    end
                end
            end

            if isempty(Nf) || fCount ~= Nf
                error('SphereField:InvalidCSTffs',...
                    'The file contained %d field blocks, but %d frequencies were expected.',...
                    fCount,Nf);
            end

            if isempty(freq)
                error('SphereField:InvalidCSTffs',...
                    'No CST frequency information was found.');
            end

            etaRad = Prad./Pacc;

            metadata = struct();
            metadata.source = "CST";
            metadata.fileFormat = "ffs";
            metadata.fileName = pathName;
            metadata.position = pos;
            metadata.xAxis = xAxis;
            metadata.zAxis = zAxis;
            metadata.acceptedPower = Pacc;
            metadata.stimulatedPower = Pstim;

            provenance = struct();
            provenance.source = "CST";
            provenance.reader = "readCSTffs";
            provenance.file = pathName;

            SF = SphereField(ph,th,Eph,Eth,freq,...
                Prad=Prad,etaRad=etaRad,r=options.r,...
                Metadata=metadata,Provenance=provenance);
        end

        function C = getPolarizationConvention()
            %GETPOLARIZATIONCONVENTION Return canonical polarization convention.

            C = struct();
            C.standard = "IEEE";
            C.timeConvention = "exp(+j*omega*t)";
            C.outgoingWave = "exp(-j*k*r)/r";
            C.circularComponentOrder = ["RHCP","LHCP"];
            C.description = ...
                "IEEE antenna polarization convention.";
        end

    end

end