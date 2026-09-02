classdef (Abstract) EMField
    %EMFIELD Abstract base class for electromagnetic field data.
    %
    % This class contains only infrastructure that is independent of the
    % geometrical representation of the field.
    %
    % Canonical frequency storage:
    %   freqHz : 1 x Nf vector [Hz]
    %
    % Derived quantities:
    %   Nf
    %   wavelength
    %   wavenumber
    %
    % Geometry, field components, interpolation, plotting and export are
    % implemented in subclasses such as SphereField and PlaneField.

    properties (SetAccess = protected)

        % Frequency samples [Hz].
        % Always stored internally as a row vector.
        freqHz (1,:) double = []

        % General descriptive metadata.
        metadata struct = struct()

        % Information describing the origin of the data.
        %
        % Example fields might include:
        %   solver
        %   filename
        %   format
        %   reader
        %   dateRead
        %   notes
        provenance struct = struct()

    end


    properties (Dependent, SetAccess = private)

        % Number of frequency samples.
        Nf
        % Free-space wavelength [m], 1 x Nf.
        wavelength
        % Free-space wavenumber [rad/m], 1 x Nf.
        wavenumber

    end


    properties (Constant)

        % Speed of light in vacuum [m/s].
        c0 = 299792458

        % Vacuum permeability [H/m].
        mu0 = 1.25663706212e-6

        % Vacuum permittivity [F/m].
        eps0 = 8.8541878128e-12

        % Free-space wave impedance [ohm].
        eta0 = 376.730313668

    end


    methods

        function obj = EMField(freqHz, options)
            %EMFIELD Constructor for common EM field infrastructure.
            %
            % Derived classes should call, for example:
            %
            %   obj@EMField(freqHz)
            %
            % or
            %
            %   obj@EMField(freqHz, ...
            %       Metadata = metadata, ...
            %       Provenance = provenance)

            arguments
                freqHz double {mustBeVector,mustBeFinite,mustBePositive}

                options.Metadata struct = struct()
                options.Provenance struct = struct()
            end

            obj.freqHz = reshape(freqHz, 1, []);
            obj.metadata = options.Metadata;
            obj.provenance = options.Provenance;
        end


        function value = get.Nf(obj)
            value = numel(obj.freqHz);
        end

        function value = get.wavelength(obj)
            value = obj.c0 ./ obj.freqHz;
        end

        function value = get.wavenumber(obj)
            value = 2*pi*obj.freqHz ./ obj.c0;
        end

        function idx = frequencyIndex(obj, frequencyHz, options)
            %FREQUENCYINDEX Find the index of a stored frequency.
            %
            % idx = obj.frequencyIndex(frequencyHz)
            %
            % Returns the index of the nearest stored frequency, provided
            % it agrees with the requested value within the specified
            % tolerance.
            %
            % Name-value options:
            %
            %   RelativeTolerance
            %       Relative tolerance. Default = 1e-9.
            %
            %   AbsoluteTolerance
            %       Absolute tolerance [Hz]. Default = 1e-6.

            arguments
                obj
                frequencyHz double {mustBeVector, mustBeFinite}

                options.RelativeTolerance (1,1) double ...
                    {mustBeNonnegative} = 1e-9

                options.AbsoluteTolerance (1,1) double ...
                    {mustBeNonnegative} = 1e-6
            end

            frequencyHz = reshape(frequencyHz, 1, []);

            idx = zeros(size(frequencyHz));

            for ii = 1:numel(frequencyHz)

                f = frequencyHz(ii);

                [df, iMin] = min(abs(obj.freqHz - f));

                tol = max( ...
                    options.AbsoluteTolerance, ...
                    options.RelativeTolerance * abs(f));

                if df > tol
                    error( ...
                        'EMField:FrequencyNotFound', ...
                        ['Frequency %.16g Hz was not found. ' ...
                         'Nearest available frequency is %.16g Hz.'], ...
                        f, obj.freqHz(iMin));
                end

                idx(ii) = iMin;
            end
        end

        function idx = nearestFrequencyIndex(obj, frequencyHz)
            %NEARESTFREQUENCYINDEX Return nearest stored frequency index.
            %
            % Unlike frequencyIndex(), no tolerance test is performed.

            arguments
                obj
                frequencyHz double {mustBeVector, mustBeFinite}
            end

            frequencyHz = reshape(frequencyHz, 1, []);

            idx = zeros(size(frequencyHz));

            for ii = 1:numel(frequencyHz)
                [~, idx(ii)] = min(abs(obj.freqHz - frequencyHz(ii)));
            end
        end

        function freqHz = getFrequency(obj, idx)
            %GETFREQUENCY Return stored frequency/frequencies by index.

            idx = obj.validateFrequencyIndex(idx);
            freqHz = obj.freqHz(idx);
        end

    end


    methods (Access = protected)

        function idx = validateFrequencyIndex(obj, idx)
            %VALIDATEFREQUENCYINDEX Validate frequency indices.
            %
            % Intended for use by subclasses.

            arguments
                obj
                idx double {mustBeVector, mustBeFinite}
            end

            if any(idx ~= fix(idx))
                error( ...
                    'EMField:InvalidFrequencyIndex', ...
                    'Frequency indices must be integers.');
            end

            if any(idx < 1 | idx > obj.Nf)
                error( ...
                    'EMField:FrequencyIndexOutOfRange', ...
                    'Frequency index must lie between 1 and %d.', ...
                    obj.Nf);
            end

            idx = reshape(idx, 1, []);
        end

        function obj = setFrequency(obj, freqHz)
            %SETFREQUENCY Replace the canonical frequency vector.
            %
            % Protected because changing frequency independently of the
            % associated field data would generally make the object
            % inconsistent.
            %
            % Since EMField uses value semantics, the updated object must
            % be returned.

            arguments
                obj
                freqHz double {mustBeVector, ...
                               mustBeFinite, ...
                               mustBePositive}
            end

            obj.freqHz = reshape(freqHz, 1, []);
        end

    end

    methods (Static)
        function freqHz = convertFrequencyToHz(freq,freqUnit)
            freqHz = freq*EMField.frequencyUnitScale(freqUnit);
        end

        function freq = convertFrequencyFromHz(freqHz,freqUnit)
            freq = freqHz/EMField.frequencyUnitScale(freqUnit);
        end
    end

    methods (Static,Access=private)
        function scale = frequencyUnitScale(freqUnit)
            arguments
                freqUnit (1,1) string
            end

            switch lower(strtrim(freqUnit))
                case "hz"
                    scale = 1;
                case "khz"
                    scale = 1e3;
                case "mhz"
                    scale = 1e6;
                case "ghz"
                    scale = 1e9;
                case "thz"
                    scale = 1e12;
                otherwise
                    error('EMField:UnsupportedFrequencyUnit',...
                        'Unsupported frequency unit "%s".',freqUnit);
            end
        end
    end

end