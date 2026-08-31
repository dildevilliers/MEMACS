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

                options.Metadata struct = struct()
                options.Provenance struct = struct()
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


            % ---------------------------------------------------------
            % Radius
            % ---------------------------------------------------------

            obj.r = options.r;


            % ---------------------------------------------------------
            % Radiated power
            % ---------------------------------------------------------

            if ~isempty(options.Prad)

                if numel(options.Prad) ~= obj.Nf
                    error( ...
                        'SphereField:PradSizeMismatch', ...
                        ['Prad must contain one value per frequency. ' ...
                        'Expected %d values, received %d.'], ...
                        obj.Nf, numel(options.Prad));
                end

                obj.Prad = reshape(options.Prad, 1, []);

            end


            % ---------------------------------------------------------
            % Detect structured grid
            % ---------------------------------------------------------

            obj.gridInfo = obj.detectStructuredGrid();

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

            [Eph_, Eth_, ~] = obj.getEfield();

            W = (abs(Eth_).^2 + abs(Eph_).^2)./(2*obj.eta0);
        end

        function U = getU(obj)
            %GETU Radiation intensity [W/sr].

            W = obj.getW();
            U = obj.r^2.*W;
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

                options.PowerSource (1,1) string ...
                    {mustBeMember(options.PowerSource, ...
                    ["auto","stored","integrated"])} = "auto"
            end

            switch options.PowerSource
                case "stored"

                    if isempty(obj.Prad)
                        error('SphereField:NoStoredPrad', ...
                            'No stored Prad is available.');
                    end

                    P = obj.Prad;
                case "integrated"

                    P = obj.integratePower();
                case "auto"

                    if ~isempty(obj.Prad)
                        P = obj.Prad;
                    else
                        P = obj.integratePower();
                    end
            end
            U = obj.getU();
            D = 4*pi.*U./reshape(P,1,[]);
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

            rHat = [ ...
                sth.*cph
                sth.*sph
                cth];

            thHat = [ ...
                cth.*cph
                cth.*sph
                -sth];

            phHat = [ ...
                -sph
                cph
                zeros(1,obj.Np)];

        end


        function ECart = getCartesianE(obj)
            %GETCARTESIANE Electric field in Cartesian coordinates.
            %
            % Output:
            %
            %   ECart : 3 x Np x Nf
            %
            % No interpolation or modification of the SphereField is performed.

            [rHat, thHat, phHat] = obj.getCartesianDirections();
            [Eph_, Eth_, Er_] = obj.getEfield();

            Eth_ = reshape(Eth_, [1 obj.Np obj.Nf]);
            Eph_ = reshape(Eph_, [1 obj.Np obj.Nf]);

            ECart = ...
                reshape(thHat,[3 obj.Np 1]).*Eth_ + ...
                reshape(phHat,[3 obj.Np 1]).*Eph_;

            if ~isempty(Er_)
                Er_ = reshape(Er_, [1 obj.Np obj.Nf]);

                ECart = ECart + ...
                    reshape(rHat,[3 obj.Np 1]).*Er_;
            end

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

                options.RelativeTolerance (1,1) double ...
                    {mustBeNonnegative} = 1e-10

                options.AbsoluteTolerance (1,1) double ...
                    {mustBeNonnegative} = 1e-12
            end

            E = obj.getCartesianE();

            details = struct();

            [details.northConsistent, details.northMaxError] = ...
                checkPole(obj.gridInfo.northPoleIndex);

            [details.southConsistent, details.southMaxError] = ...
                checkPole(obj.gridInfo.southPoleIndex);

            isConsistent = ...
                details.northConsistent && ...
                details.southConsistent;


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

                tol = max( ...
                    options.AbsoluteTolerance, ...
                    options.RelativeTolerance .* refScale);

                tf = all(maxError <= tol);

            end

        end

        function G = getGridView(obj, viewType)
            %GETGRIDVIEW Return a numerical view of the canonical Ph-Th samples.
            %
            % G = obj.getGridView("stored")
            % G = obj.getGridView("singlePeriod")
            %
            % "stored"
            %   Returns the canonical samples in their stored form.
            %   Structured grids are returned as [Nth x Nph x Nf] arrays.
            %   Unstructured grids remain in point form.
            %
            % "singlePeriod"
            %   Returns one complete phi period for a structured grid with full
            %   phi coverage. A repeated endpoint is retained. Extra revolutions
            %   are removed. If the stored grid omits the periodic endpoint, it is
            %   appended from the equivalent first phi cut.
            %
            % No interpolation is performed.

            arguments
                obj
                viewType (1,1) string {mustBeMember(viewType, ...
                    ["stored","singlePeriod"])} = "stored"
            end

            switch viewType

                case "stored"
                    G = obj.getStoredGridView();

                case "singlePeriod"
                    G = obj.getSinglePeriodGridView();

            end

        end

        function PradInt = integratePower(obj)
            %INTEGRATEPOWER Integrate radiated power over the full sphere.
            %
            %   PradInt = obj.integratePower()
            %
            % Returns:
            %   PradInt : 1 x Nf radiated power [W]
            %
            % The field must be sampled on a structured grid covering the full
            % sphere. Integration is performed over one complete phi period.
            %
            % The time-average radial power density is assumed to be
            %
            %   S_r = (|Eth|^2 + |Eph|^2)/(2*eta0)
            %
            % corresponding to complex peak-amplitude electric-field phasors.
            %
            % Er does not contribute to the radiated far-field power.

            if ~obj.isStructured
                error( ...
                    'SphereField:PowerIntegrationRequiresStructuredGrid', ...
                    ['Power integration currently requires a structured ' ...
                    'Ph-Th grid.']);
            end

            if ~obj.isFullSphere
                error( ...
                    'SphereField:PowerIntegrationRequiresFullSphere', ...
                    ['Power integration requires complete spherical ' ...
                    'coverage.']);
            end

            % Exactly one complete phi period, including the periodic endpoint.
            G = obj.getGridView("singlePeriod");

            % Angular coordinates in radians.
            th_ = deg2rad(G.thVec);
            ph_ = deg2rad(G.phVec);

            % Transverse field magnitude squared - equivalent to getU (but
            % without the correct grid)
            E2 = abs(G.Eth).^2 + abs(G.Eph).^2;

            % Spherical Jacobian.
            sinTh = reshape(sin(th_), [], 1, 1);

            integrand = E2.*sinTh;

            % Integrate first over theta, then phi.
            intTh = trapz(th_, integrand, 1);
            intOmega = trapz(ph_, intTh, 2);

            intOmega = reshape(intOmega, 1, []);

            % Convert angular integral of |E|^2 to radiated power.
            PradInt = intOmega./(2*obj.eta0);

        end

    end


    methods (Access = private)

        function E = normaliseFieldArray(obj, E, Np_, Nf_, name)
            %NORMALISEFIELDARRAY Convert field component to Np x Nf.

            if isempty(E)
                error( ...
                    'SphereField:EmptyFieldComponent', ...
                    '%s cannot be empty.', ...
                    name);
            end


            % Single-frequency case: any vector containing Np samples is
            % acceptable.
            if Nf_ == 1 && isvector(E) && numel(E) == Np_

                E = E(:);
                return

            end


            % General multi-frequency form.
            if ~ismatrix(E) || ...
                    size(E,1) ~= Np_ || ...
                    size(E,2) ~= Nf_

                error( ...
                    'SphereField:FieldSizeMismatch', ...
                    ['%s must have size Np x Nf = %d x %d. ' ...
                    'Received %s.'], ...
                    name, Np_, Nf_, mat2str(size(E)));
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


            %% Structured tensor-product grid

            if Nph_ * Nth_ ~= obj.Np
                return
            end

            [tfPh, iPh] = ismember(obj.ph, phVec);
            [tfTh, iTh] = ismember(obj.th, thVec);

            if ~all(tfPh) || ~all(tfTh)
                return
            end

            gridIndex = sub2ind( ...
                [Nth_, Nph_], ...
                iTh, ...
                iPh);

            if numel(unique(gridIndex)) ~= obj.Np
                return
            end


            info.isStructured = true;

            info.phVec = reshape(phVec,1,[]);
            info.thVec = reshape(thVec,1,[]);

            info.Nph = Nph_;
            info.Nth = Nth_;

            info.gridIndex = gridIndex;


            %% Theta coverage

            angleTol = 1e-8;

            info.northPoleIndex = find(abs(obj.th) <= angleTol);
            info.southPoleIndex = find(abs(obj.th - 180) <= angleTol);

            info.hasNorthPole = ~isempty(info.northPoleIndex);
            info.hasSouthPole = ~isempty(info.southPoleIndex);

            info.coversFullTheta = ...
                abs(thVec(1)) <= angleTol && ...
                abs(thVec(end) - 180) <= angleTol;


            %% Phi periodicity and coverage

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

                    gapTol = max( ...
                        angleTol, ...
                        1e-8 * maxInteriorGap);

                    info.coversFullPhi = ...
                        seamGap <= maxInteriorGap + gapTol;

                end

            end

            
            %% Complete sphere

            info.isFullSphere = ...
                info.coversFullPhi && ...
                info.coversFullTheta;

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

                G.hasPhiPeriodicEndpoint = ...
                    obj.hasPhiPeriodicEndpoint;

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
                error( ...
                    'SphereField:SinglePeriodRequiresStructuredGrid', ...
                    ['singlePeriod is only defined for structured ' ...
                    'Ph-Th grids.']);
            end

            if ~obj.coversFullPhi
                error( ...
                    'SphereField:SinglePeriodRequiresFullPhi', ...
                    ['singlePeriod requires complete phi coverage.']);
            end

            G0 = obj.getStoredGridView();

            phVec = G0.phVec;

            % Use the first stored phi as the start of the period.
            ph0 = phVec(1);
            phEnd = ph0 + 360;

            angleTol = 1e-8;

            % Keep all stored cuts belonging to the first complete revolution.
            keep = ...
                phVec >= ph0 - angleTol & ...
                phVec <= phEnd + angleTol;

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


    end

end