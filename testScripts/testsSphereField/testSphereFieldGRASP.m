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

        function testGRASPAgainstFarField(testCase,graspFile)

            pathName = fullfile(testCase.dataPath,graspFile);

            SF = SphereField.readGRASPgrd(pathName);
            FF = FarField.readGRASPgrd(pathName);

            % Convert legacy FarField to spherical coordinates and linear
            % polarization.
            FF = FF.coor2spherical(true);
            FF = FF.pol2linear;

            %--------------------------------------------------------------
            % Frequency
            %--------------------------------------------------------------
            testCase.verifyEqual(...
                SF.freqHz,...
                FF.freqHz,...
                'RelTol',1e-12);

            %--------------------------------------------------------------
            % Coordinates
            %--------------------------------------------------------------
            thSF = SF.th;
            phSF = SF.ph;

            thFF = rad2deg(FF.y);
            phFF = rad2deg(FF.x);

            %--------------------------------------------------------------
            % Remove pole samples.
            %
            % Phi and the spherical transverse basis are undefined at the
            % poles. The legacy and new implementations are therefore not
            % required to associate the same phi/basis with these samples.
            %--------------------------------------------------------------
            poleTol = 1e-8;

            regularSF = abs(sind(thSF)) > poleTol;
            regularFF = abs(sind(thFF)) > poleTol;

            thSF = thSF(regularSF);
            phSF = phSF(regularSF);
            EthSF = SF.Eth(regularSF,:);
            EphSF = SF.Eph(regularSF,:);

            thFF = thFF(regularFF);
            phFF = phFF(regularFF);
            EthFF = FF.E1(regularFF,:);
            EphFF = FF.E2(regularFF,:);

            testCase.verifyEqual(numel(thSF),numel(thFF));

            %--------------------------------------------------------------
            % Generate robust matching keys.
            %
            % FarField has undergone several radian/degree and coordinate
            % transformations, so round the angles before using them as
            % sorting keys.
            %--------------------------------------------------------------
            nDigits = 10;

            thKeySF = round(thSF,nDigits,'significant');
            thKeyFF = round(thFF,nDigits,'significant');

            phKeySF = round(phSF,nDigits,'significant');
            phKeyFF = round(phFF,nDigits,'significant');

            % Put phi into the same periodic interval.
            phKeySF = mod(phKeySF + 180,360) - 180;
            phKeyFF = mod(phKeyFF + 180,360) - 180;

            % Explicitly collapse values numerically at the +/-180 seam.
            seamTol = 1e-8;

            phKeySF(abs(abs(phKeySF) - 180) < seamTol) = -180;
            phKeyFF(abs(abs(phKeyFF) - 180) < seamTol) = -180;

            keySF = [thKeySF,phKeySF];
            keyFF = [thKeyFF,phKeyFF];

            [keySF,ia] = sortrows(keySF,[1 2]);
            [keyFF,ib] = sortrows(keyFF,[1 2]);

            %--------------------------------------------------------------
            % Verify that the same nominal angular samples were found.
            %--------------------------------------------------------------
            testCase.verifyEqual(...
                keySF,...
                keyFF,...
                'AbsTol',1e-8);

            %--------------------------------------------------------------
            % Compare the original coordinates after matching.
            %--------------------------------------------------------------
            thSF = thSF(ia);
            phSF = phSF(ia);

            thFF = thFF(ib);
            phFF = phFF(ib);

            testCase.verifyEqual(...
                thSF,...
                thFF,...
                'AbsTol',1e-6);

            % Compare phi using a periodic angular difference rather than
            % comparing its particular numerical representation.
            dph = mod(phSF - phFF + 180,360) - 180;

            testCase.verifyEqual(...
                dph,...
                zeros(size(dph)),...
                'AbsTol',1e-6);

            %--------------------------------------------------------------
            % Fields
            %--------------------------------------------------------------
            EthSF = EthSF(ia,:);
            EphSF = EphSF(ia,:);

            EthFF = EthFF(ib,:);
            EphFF = EphFF(ib,:);

            tolAbs = 1e-10;
            tolRel = 1e-9;

            % Treat tiny numerical residue as zero.
            EthSF(abs(EthSF) < tolAbs) = 0;
            EphSF(abs(EphSF) < tolAbs) = 0;
            EthFF(abs(EthFF) < tolAbs) = 0;
            EphFF(abs(EphFF) < tolAbs) = 0;

            testCase.verifyEqual(...
                EthSF,...
                EthFF,...
                'RelTol',tolRel,...
                'AbsTol',tolAbs);

            testCase.verifyEqual(...
                EphSF,...
                EphFF,...
                'RelTol',tolRel,...
                'AbsTol',tolAbs);
        end


        function testGRASPHDF5AgainstASCII(testCase)
            % The same GRASP calculation stored in ASCII and HDF5 should
            % produce equivalent SphereField objects.

            grdPath = fullfile(testCase.dataPath,'FF_feed.grd');
            h5Path = fullfile(testCase.dataPath,'FF_feed.h5');

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