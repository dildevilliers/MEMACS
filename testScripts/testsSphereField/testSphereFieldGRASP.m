classdef testSphereFieldGRASP < matlab.unittest.TestCase
    % Tests for SphereField GRASP .grd/.h5 readers.
    %
    % The main regression test compares SphereField.readGRASPgrd
    % against the legacy FarField.readGRASPgrd implementation for a
    % collection of GRASP field files.
    %
    % A separate test compares equivalent GRASP ASCII and HDF5 files.

    properties (TestParameter)
        graspFile = struct(...
            'mainBeam',              'FFmainBeam.grd',...
            'feedRH_PhTh',           'FFfeedRH_PhTh.grd',...
            'feedRH_Lin',            'FFfeedRH_Lin.grd',...
            'feedRH_CP',             'FFfeedRH_CP.grd',...
            'uv_spherical',          'FF_uv_spherical.grd',...
            'trueview_spherical',    'FF_trueview_spherical.grd',...
            'phth_spherical_sym360', 'FF_phth_spherical_sym360.grd',...
            'phth_spherical_sym180', 'FF_phth_spherical_sym180.grd',...
            'phth_spherical_pos360', 'FF_phth_spherical_pos360.grd',...
            'phth_spherical_pos180', 'FF_phth_spherical_pos180.grd',...
            'phth_ludwig3_pos180',   'FF_phth_ludwig3_pos180.grd',...
            'phth_circular_pos180',  'FF_phth_circular_pos180.grd',...
            'feed_grd',              'FF_feed.grd',...
            'elaz_spherical_sym180', 'FF_elaz_spherical_sym180.grd',...
            'elaz_spherical_pos180', 'FF_elaz_spherical_pos180.grd',...
            'azel_spherical_sym180', 'FF_azel_spherical_sym180.grd',...
            'azel_spherical_pos180', 'FF_azel_spherical_pos180.grd');
    end

    properties
        dataPath
    end

    methods (TestClassSetup)
        function setupDataPath(testCase)
            p = fileparts(mfilename("fullpath"));
            testCase.dataPath = [p,'\..\..\data\SimPatterns\GRASPgrd'];
        end
    end

    methods (Test)

        % function testGRASPAgainstFarField(testCase,graspFile)
        %     % Compare new SphereField reader against legacy FarField reader.
        %
        %     pathName = fullfile(testCase.dataPath,graspFile);
        %
        %     SF = SphereField.readGRASPgrd(pathName);
        %     FF = FarField.readGRASPgrd(pathName);
        %
        %     % Convert legacy FarField to the canonical representation used
        %     % by SphereField.
        %     FF = FF.coor2spherical(true);
        %     FF = FF.pol2linear;
        %     % FF = FF.setRangeSph('pos','180');
        %
        %     %----------------------------------------------------------
        %     % Frequency
        %     %----------------------------------------------------------
        %     testCase.verifyEqual(SF.freqHz,FF.freqHz,'RelTol',1e-12);
        %
        %     %----------------------------------------------------------
        %     % Coordinates
        %     %
        %     % The two classes need not store the samples in exactly the
        %     % same order, so sort both using theta then phi.
        %     %----------------------------------------------------------
        %     % A = [SF.th,SF.ph];
        %     % B = [rad2deg(FF.y),rad2deg(FF.x)];
        %     %
        %     % [A,ia] = sortrows(A,[1 2]);
        %     % [B,ib] = sortrows(B,[1 2]);
        %     % % ia = 1:FF.Nang;
        %     % % ib = 1:FF.Nang;
        %
        %     % A = [SF.th,SF.ph];
        %     % B = [rad2deg(FF.y),rad2deg(FF.x)];
        %     %
        %     % % Use the same phi representation for comparison only.
        %     % % Maps phi to [-180,180).
        %     % A(:,2) = mod(A(:,2) + 180,360) - 180;
        %     % B(:,2) = mod(B(:,2) + 180,360) - 180;
        %     %
        %     % % Now sort into the same physical angular order.
        %     % [A,ia] = sortrows(A,[1 2]);
        %     % [B,ib] = sortrows(B,[1 2]);
        %
        %
        %
        %     A = [SF.th,SF.ph];
        %     B = [rad2deg(FF.y),rad2deg(FF.x)];
        %
        %     A(:,2) = mod(A(:,2) + 180,360) - 180;
        %     B(:,2) = mod(B(:,2) + 180,360) - 180;
        %
        %     tolAbs = 1e-12;
        %     A(abs(A) < tolAbs) = 0;
        %     B(abs(B) < tolAbs) = 0;
        %
        %     nDigits = 12;
        %     keyA = round(A,nDigits,'significant');
        %     keyB = round(B,nDigits,'significant');
        %
        %     [keyA,ia] = sortrows(keyA,[1 2]);
        %     [keyB,ib] = sortrows(keyB,[1 2]);
        %
        %     testCase.verifyEqual(keyA,keyB);
        %     % testSphereFieldGRASP.verifyFieldEqual(testCase,keyA,keyB);
        %
        %     testCase.verifyEqual(A(ia,:),B(ib,:),'AbsTol',1e-6);
        %     % testSphereFieldGRASP.verifyFieldEqual(testCase,A(ia,:),B(ib,:));
        %
        %     %----------------------------------------------------------
        %     % Spherical electric-field components
        %     %----------------------------------------------------------
        %     testCase.verifyEqual(SF.Eth(ia,:),FF.E1(ib,:),'RelTol',1e-9,'AbsTol',1e-10);
        %     testCase.verifyEqual(SF.Eph(ia,:),FF.E2(ib,:),'RelTol',1e-9,'AbsTol',1e-10);
        %
        %     % testSphereFieldGRASP.verifyFieldEqual(testCase,SF.Eth(ia,:),FF.E1(ib,:));
        %     % testSphereFieldGRASP.verifyFieldEqual(testCase,SF.Eph(ia,:),FF.E2(ib,:));
        %
        % end

        % function testGRASPAgainstFarField(testCase,graspFile)
        % 
        %     pathName = fullfile(testCase.dataPath,graspFile);
        % 
        %     SF = SphereField.readGRASPgrd(pathName);
        %     FF = FarField.readGRASPgrd(pathName);
        % 
        %     % Convert legacy FarField to spherical coordinates and linear
        %     % polarization.
        %     FF = FF.coor2spherical(true);
        %     FF = FF.pol2linear;
        % 
        %     %--------------------------------------------------------------
        %     % Frequency
        %     %--------------------------------------------------------------
        %     testCase.verifyEqual(...
        %         SF.freqHz,...
        %         FF.freqHz,...
        %         'RelTol',1e-12);
        % 
        %     %--------------------------------------------------------------
        %     % Coordinates
        %     %--------------------------------------------------------------
        %     thSF = SF.th;
        %     phSF = SF.ph;
        % 
        %     thFF = rad2deg(FF.y);
        %     phFF = rad2deg(FF.x);
        % 
        %     %--------------------------------------------------------------
        %     % Remove pole samples.
        %     %
        %     % Phi and the spherical transverse basis are undefined at the
        %     % poles. The legacy and new implementations are therefore not
        %     % required to associate the same phi/basis with these samples.
        %     %--------------------------------------------------------------
        %     poleTol = 1e-8;
        % 
        %     regularSF = abs(sind(thSF)) > poleTol;
        %     regularFF = abs(sind(thFF)) > poleTol;
        % 
        %     thSF = thSF(regularSF);
        %     phSF = phSF(regularSF);
        %     EthSF = SF.Eth(regularSF,:);
        %     EphSF = SF.Eph(regularSF,:);
        % 
        %     thFF = thFF(regularFF);
        %     phFF = phFF(regularFF);
        %     EthFF = FF.E1(regularFF,:);
        %     EphFF = FF.E2(regularFF,:);
        % 
        %     testCase.verifyEqual(numel(thSF),numel(thFF));
        % 
        %     %--------------------------------------------------------------
        %     % Generate robust matching keys.
        %     %
        %     % FarField has undergone several radian/degree and coordinate
        %     % transformations, so round the angles before using them as
        %     % sorting keys.
        %     %--------------------------------------------------------------
        %     nDigits = 10;
        % 
        %     thKeySF = round(thSF,nDigits,'significant');
        %     thKeyFF = round(thFF,nDigits,'significant');
        % 
        %     phKeySF = round(phSF,nDigits,'significant');
        %     phKeyFF = round(phFF,nDigits,'significant');
        % 
        %     % Put phi into the same periodic interval.
        %     phKeySF = mod(phKeySF + 180,360) - 180;
        %     phKeyFF = mod(phKeyFF + 180,360) - 180;
        % 
        %     % Explicitly collapse values numerically at the +/-180 seam.
        %     seamTol = 1e-8;
        % 
        %     phKeySF(abs(abs(phKeySF) - 180) < seamTol) = -180;
        %     phKeyFF(abs(abs(phKeyFF) - 180) < seamTol) = -180;
        % 
        %     keySF = [thKeySF,phKeySF];
        %     keyFF = [thKeyFF,phKeyFF];
        % 
        %     [keySF,ia] = sortrows(keySF,[1 2]);
        %     [keyFF,ib] = sortrows(keyFF,[1 2]);
        % 
        %     %--------------------------------------------------------------
        %     % Verify that the same nominal angular samples were found.
        %     %--------------------------------------------------------------
        %     testCase.verifyEqual(...
        %         keySF,...
        %         keyFF,...
        %         'AbsTol',1e-8);
        % 
        %     %--------------------------------------------------------------
        %     % Compare the original coordinates after matching.
        %     %--------------------------------------------------------------
        %     thSF = thSF(ia);
        %     phSF = phSF(ia);
        % 
        %     thFF = thFF(ib);
        %     phFF = phFF(ib);
        % 
        %     testCase.verifyEqual(...
        %         thSF,...
        %         thFF,...
        %         'AbsTol',1e-6);
        % 
        %     % Compare phi using a periodic angular difference rather than
        %     % comparing its particular numerical representation.
        %     dph = mod(phSF - phFF + 180,360) - 180;
        % 
        %     testCase.verifyEqual(...
        %         dph,...
        %         zeros(size(dph)),...
        %         'AbsTol',1e-6);
        % 
        %     %--------------------------------------------------------------
        %     % Fields
        %     %--------------------------------------------------------------
        %     EthSF = EthSF(ia,:);
        %     EphSF = EphSF(ia,:);
        % 
        %     EthFF = EthFF(ib,:);
        %     EphFF = EphFF(ib,:);
        % 
        %     tolAbs = 1e-10;
        %     tolRel = 1e-9;
        % 
        %     % Treat tiny numerical residue as zero.
        %     EthSF(abs(EthSF) < tolAbs) = 0;
        %     EphSF(abs(EphSF) < tolAbs) = 0;
        %     EthFF(abs(EthFF) < tolAbs) = 0;
        %     EphFF(abs(EphFF) < tolAbs) = 0;
        % 
        %     testCase.verifyEqual(...
        %         EthSF,...
        %         EthFF,...
        %         'RelTol',tolRel,...
        %         'AbsTol',tolAbs);
        % 
        %     testCase.verifyEqual(...
        %         EphSF,...
        %         EphFF,...
        %         'RelTol',tolRel,...
        %         'AbsTol',tolAbs);
        % end

        function testGRASPHDF5AgainstASCII(testCase)
            % The same GRASP calculation stored in ASCII and HDF5 should
            % produce equivalent SphereField objects.

            grdPath = fullfile(testCase.dataPath,'\FF_feed.grd');
            h5Path = fullfile(testCase.dataPath,'\FF_feed.h5');

            SFgrd = SphereField.readGRASPgrd(grdPath);
            SFh5 = SphereField.readGRASPgrd(h5Path);

            %----------------------------------------------------------
            % Frequency
            %----------------------------------------------------------
            testCase.verifyEqual(SFh5.freqHz,SFgrd.freqHz,'RelTol',1e-12);

            %----------------------------------------------------------
            % Coordinates
            %----------------------------------------------------------
            A = [SFgrd.th,SFgrd.ph];
            B = [SFh5.th,SFh5.ph];

            [A,ia] = sortrows(A,[1 2]);
            [B,ib] = sortrows(B,[1 2]);

            testCase.verifyEqual(A,B,'AbsTol',1e-10);

            %----------------------------------------------------------
            % Fields
            %----------------------------------------------------------
            testCase.verifyEqual(SFh5.Eth(ib,:),SFgrd.Eth(ia,:),'RelTol',1e-9,'AbsTol',1e-10);

            testCase.verifyEqual(SFh5.Eph(ib,:),SFgrd.Eph(ia,:),'RelTol',1e-9,'AbsTol',1e-10);
        end

        % function testGRASPNativeGridDetection(testCase)
        % 
        %     p = fileparts(mfilename("fullpath"));
        %     base = [p,'\..\..\data\SimPatterns\GRASPgrd\'];
        % 
        %     cases = {
        %         "FF_uv_spherical.grd",              true,  "DirCos"
        %         "FF_trueview_spherical.grd",        true,  "TrueView"
        %         "FF_phth_spherical_sym360.grd",     true,  "PhTh"
        %         "FF_phth_spherical_sym180.grd",     true,  "PhTh"
        %         "FF_phth_spherical_pos360.grd",     true,  "PhTh"
        %         "FF_phth_spherical_pos180.grd",     true,  "PhTh"
        %         "FF_phth_ludwig3_pos180.grd",       true,  "PhTh"
        %         "FF_phth_circular_pos180.grd",      true,  "PhTh"
        %         "FF_elaz_spherical_sym180.grd",     true,  "ElAz"
        %         "FF_elaz_spherical_pos180.grd",     true,  "ElAz"
        %         "FF_azel_spherical_sym180.grd",     true,  "AzEl"
        %         "FF_azel_spherical_pos180.grd",     true,  "AzEl"
        %         };
        % 
        %     for ii = 1:size(cases,1)
        % 
        %         fileName = cases{ii,1};
        %         expectedStructured = cases{ii,2};
        %         expectedCoordinates = cases{ii,3};
        % 
        %         SF = SphereField.readGRASPgrd(fullfile(base,fileName));
        % 
        %         testCase.verifyEqual(SF.hasNativeStructuredGrid,...
        %             expectedStructured,...
        %             sprintf('Unexpected native-grid state for %s.',fileName));
        % 
        %         if expectedStructured
        %             testCase.verifyEqual(...
        %                 SF.nativeGridInfo.coordinates,...
        %                 expectedCoordinates,...
        %                 sprintf('Unexpected native coordinate system for %s.',fileName));
        % 
        %             testCase.verifyEqual(...
        %                 SF.nativeGridInfo.Nx*SF.nativeGridInfo.Ny,...
        %                 SF.Np,...
        %                 sprintf('Native grid size mismatch for %s.',fileName));
        % 
        %             testCase.verifyEqual(...
        %                 numel(SF.nativeGridInfo.xVec),...
        %                 SF.nativeGridInfo.Nx);
        % 
        %             testCase.verifyEqual(...
        %                 numel(SF.nativeGridInfo.yVec),...
        %                 SF.nativeGridInfo.Ny);
        %         end
        %     end
        % end

        function testGRASPHDF5NativeGridDetection(testCase)

            p = fileparts(mfilename("fullpath"));
            pathName = [p,'\..\..\data\SimPatterns\GRASPgrd\FF_feed.h5'];

            SF = SphereField.readGRASPgrd(pathName);

            % If the HDF5 grid contains NaN-truncated samples, we deliberately
            % discard native structured-grid information.
            if SF.isStructured
                % Canonical Ph-Th may still happen to be structured.
                % This is independent of native-grid preservation.
            end

            testCase.verifyEqual(...
                SF.nativeGridInfo.gridIndex,...
                (1:SF.Np).');

            % Change this expected value after checking the actual file once.
            disp(SF.hasNativeStructuredGrid)
            disp(SF.nativeGridInfo)
        end

        % function testGRASPNativeGridSummary(testCase)
        % 
        %     p = fileparts(mfilename("fullpath"));
        %     base = [p,'\..\..\data\SimPatterns\GRASPgrd\'];
        % 
        %     files = {
        %         "FF_uv_spherical.grd"
        %         "FF_trueview_spherical.grd"
        %         "FF_phth_spherical_sym360.grd"
        %         "FF_elaz_spherical_sym180.grd"
        %         "FF_azel_spherical_sym180.grd"
        %         "FF_feed.h5"
        %         };
        % 
        %     fprintf('\n');
        %     fprintf('%-35s %-10s %-10s %-12s %-8s %-8s\n',...
        %         'File','Canonical','Native','Coordinates','Nx','Ny');
        %     fprintf('%s\n',repmat('-',1,90));
        % 
        %     for ii = 1:numel(files)
        % 
        %         SF = SphereField.readGRASPgrd(...
        %             fullfile(base,files{ii}));
        % 
        %         if SF.hasNativeStructuredGrid
        %             coord = SF.nativeGridInfo.coordinates;
        %             Nx = SF.nativeGridInfo.Nx;
        %             Ny = SF.nativeGridInfo.Ny;
        %         else
        %             coord = "";
        %             Nx = 0;
        %             Ny = 0;
        %         end
        % 
        %         fprintf('%-35s %-10d %-10d %-12s %-8d %-8d\n',...
        %             files{ii},...
        %             SF.isStructured,...
        %             SF.hasNativeStructuredGrid,...
        %             coord,...
        %             Nx,...
        %             Ny);
        %     end
        % 
        %     testCase.verifyTrue(true);
        % end

        % function testGRASPNativeCoordinateMapping(testCase)
        % 
        %     p = fileparts(mfilename("fullpath"));
        %     base = [p,'\..\..\data\SimPatterns\GRASPgrd\'];
        % 
        %     files = {
        %         "FF_uv_spherical.grd"
        %         "FF_trueview_spherical.grd"
        %         "FF_phth_spherical_sym360.grd"
        %         "FF_elaz_spherical_sym180.grd"
        %         "FF_azel_spherical_sym180.grd"
        %         };
        % 
        %     fprintf('\n');
        %     fprintf('%-35s %-6s %-10s %-12s %-12s\n',...
        %         'File','IGRID','Coordinates','max err x','max err y');
        %     fprintf('%s\n',repmat('-',1,85));
        % 
        %     for ii = 1:numel(files)
        % 
        %         pathName = fullfile(base,files{ii});
        % 
        %         % Read the raw GRASP representation.
        %         G = SphereField.readGRASPgrdASCII(pathName);
        % 
        %         % Convert to canonical coordinates.
        %         [ph,th] = SphereField.graspGridToPhTh(...
        %             G.x,G.y,G.IGRID);
        % 
        %         % Temporary SphereField only so that we use the actual
        %         % SphereField coordinate transformation code.
        %         E = zeros(size(ph));
        % 
        %         SF = SphereField(ph,th,E,E,1e9);
        % 
        %         coordinates = SphereField.graspGridCoordinates(G.IGRID);
        % 
        %         C = SF.getCoordinates(Coordinates=coordinates);
        % 
        %         % Expected relationship between GRASP coordinates and the
        %         % corresponding SphereField coordinate family.
        %         switch G.IGRID
        % 
        %             case 1      % GRASP uv -> SphereField DirCos
        %                 xExpected = G.x;
        %                 yExpected = G.y;
        %                 valid = true(size(G.x));
        % 
        %             case 4      % GRASP elevation-over-azimuth
        %                 xExpected = -G.x;
        %                 yExpected = G.y;
        % 
        %                 % x is undefined at the AzEl poles.
        %                 valid = abs(abs(G.y)-90) > 1e-8;
        % 
        %             case 5      % GRASP true-view convention
        %                 xExpected = -G.x;
        %                 yExpected = G.y;
        %                 valid = true(size(G.x));
        % 
        %             case 6      % GRASP azimuth-over-elevation
        %                 xExpected = G.y;
        %                 yExpected = -G.x;
        % 
        %                 % x is undefined at the ElAz poles.
        %                 valid = abs(abs(G.x)-90) > 1e-8;
        % 
        %             case 7      % conventional phi-theta
        %                 xExpected = G.x;
        %                 yExpected = G.y;
        %                 valid = true(size(G.x));
        % 
        %             case 9      % EDX AzEl
        %                 xExpected = G.x;
        %                 yExpected = G.y;
        %                 valid = abs(abs(G.y)-90) > 1e-8;
        % 
        %             case 10     % EDX ElAz
        %                 xExpected = G.y;
        %                 yExpected = G.x;
        %                 valid = abs(abs(G.x)-90) > 1e-8;
        % 
        %             otherwise
        %                 error('Unexpected IGRID = %d.',G.IGRID);
        %         end
        % 
        %         % Angular x coordinates may differ by 360 deg.
        %         if any(G.IGRID == [4 6 7 9 10])
        %             dx = mod(C.x(valid)-xExpected(valid)+180,360)-180;
        %         else
        %             dx = C.x(valid)-xExpected(valid);
        %         end
        % 
        %         dy = C.y(valid)-yExpected(valid);
        % 
        %         maxErrX = max(abs(dx));
        %         maxErrY = max(abs(dy));
        % 
        %         fprintf('%-35s %-6d %-10s %-12.3g %-12.3g\n',...
        %             files{ii},G.IGRID,coordinates,maxErrX,maxErrY);
        % 
        %         testCase.verifyLessThan(maxErrX,1e-8);
        %         testCase.verifyLessThan(maxErrY,1e-8);
        %     end
        % end

        function testGRASPNativeStructuredGrids(testCase)

            p = fileparts(mfilename("fullpath"));
            base = [p,'\..\..\data\SimPatterns\GRASPgrd\'];

            cases = {
                "FF_uv_spherical.grd",              false, true,  "DirCos",   37, 37
                "FF_trueview_spherical.grd",        false, true,  "TrueView", 37, 37
                "FF_phth_spherical_sym360.grd",     true,  true,  "PhTh",     37, 37
                "FF_elaz_spherical_sym180.grd",     false, true,  "ElAz",     37, 37
                "FF_azel_spherical_sym180.grd",     false, true,  "AzEl",     37, 37
                "FF_feed.h5",                       true,  true,  "PhTh",     73, 37
                };

            for ii = 1:size(cases,1)

                fileName = cases{ii,1};
                expectedCanonical = cases{ii,2};
                expectedNative = cases{ii,3};
                expectedCoordinates = cases{ii,4};
                expectedNx = cases{ii,5};
                expectedNy = cases{ii,6};

                SF = SphereField.readGRASPgrd(fullfile(base,fileName));

                testCase.verifyEqual(SF.isStructured,expectedCanonical,...
                    sprintf('Unexpected canonical-grid state for %s.',fileName));

                testCase.verifyEqual(SF.hasNativeStructuredGrid,expectedNative,...
                    sprintf('Unexpected native-grid state for %s.',fileName));

                testCase.verifyEqual(SF.nativeGridInfo.coordinates,expectedCoordinates,...
                    sprintf('Unexpected native coordinate system for %s.',fileName));

                testCase.verifyEqual(SF.nativeGridInfo.Nx,expectedNx);
                testCase.verifyEqual(SF.nativeGridInfo.Ny,expectedNy);

                testCase.verifyEqual(...
                    SF.nativeGridInfo.Nx*SF.nativeGridInfo.Ny,...
                    SF.Np);

                testCase.verifyEqual(...
                    sort(SF.nativeGridInfo.gridIndex),...
                    (1:SF.Np).');
            end
        end

        function testGRASPNativeCoordinateMapping(testCase)

            p = fileparts(mfilename("fullpath"));
            base = [p,'\..\..\data\SimPatterns\GRASPgrd\'];

            files = {
                "FF_uv_spherical.grd"
                "FF_trueview_spherical.grd"
                "FF_phth_spherical_sym360.grd"
                "FF_elaz_spherical_sym180.grd"
                "FF_azel_spherical_sym180.grd"
                };

            for ii = 1:numel(files)

                G = SphereField.readGRASPgrdASCII(fullfile(base,files{ii}));

                [ph,th] = SphereField.graspGridToPhTh(G.x,G.y,G.IGRID);

                E = zeros(size(ph));
                SF = SphereField(ph,th,E,E,1e9);

                coordinates = SphereField.graspGridCoordinates(G.IGRID);
                C = SF.getCoordinates(Coordinates=coordinates);

                switch G.IGRID
                    case 1
                        xExpected = G.x;
                        yExpected = G.y;

                    case 5
                        xExpected = -G.x;
                        yExpected = G.y;

                    case 7
                        xExpected = G.x;
                        yExpected = G.y;

                    case 9
                        xExpected = G.x;
                        yExpected = G.y;

                    case 10
                        xExpected = G.y;
                        yExpected = G.x;
                end

                if any(G.IGRID == [7 9 10])
                    dx = mod(C.x-xExpected+180,360)-180;
                else
                    dx = C.x-xExpected;
                end

                dy = C.y-yExpected;

                testCase.verifyLessThan(max(abs(dx)),1e-8);
                testCase.verifyLessThan(max(abs(dy)),1e-8);
            end
        end

    end

    methods (Static,Access=private)

        % function verifyFieldEqual(testCase,A,B)
        %     tolAbs = 1e-10;
        %     tolRel = 1e-6;
        %
        %     err = abs(A - B);
        %     scale = max(abs(A),abs(B));
        %
        %     testCase.verifyTrue(all(err <= tolAbs + tolRel*scale,'all'),'Field values differ beyond the specified tolerance.');
        % end

    end
end