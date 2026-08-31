function tests = testSphereField
%TESTSPHEREFIELD Unit tests for SphereField.

    tests = functiontests(localfunctions);

end


%% Basic construction

function testBasicConstructor(testCase)

    ph = [0; 90; 180];
    th = [45; 45; 45];

    Eph = [1; 2; 3];
    Eth = [4; 5; 6];

    SF = SphereField( ...
        ph, th, Eph, Eth, 10e9);

    verifyEqual(testCase, SF.ph, ph);
    verifyEqual(testCase, SF.th, th);

    verifyEqual(testCase, SF.Eph, Eph);
    verifyEqual(testCase, SF.Eth, Eth);

    verifyEqual(testCase, SF.freqHz, 10e9);

    verifyEqual(testCase, SF.Np, 3);
    verifyEqual(testCase, SF.Nf, 1);

end

function testDefaultConstructor(testCase)

    SF = SphereField;

    verifyEqual(testCase, SF.Nf, 1);
    verifyEqual(testCase, SF.freqHz, 1e9);

    verifyTrue(testCase, SF.isStructured);

    verifyEqual(testCase, SF.Nph, 36);
    verifyEqual(testCase, SF.Nth, 37);

    verifyEqual(testCase, SF.Np, 36*37);

    verifySize(testCase, SF.Eth, [SF.Np 1]);
    verifySize(testCase, SF.Eph, [SF.Np 1]);

    verifyFalse(testCase, SF.hasEr);

    verifyTrue(testCase, all(isfinite(SF.Eth)));
    verifyTrue(testCase, all(isfinite(SF.Eph)));

end

function testIncompleteConstructor(testCase)

    verifyError(testCase, ...
        @() SphereField([0;10], [0;10], [1;1]), ...
        'SphereField:NotEnoughInputs');

end


function testCoordinateShapeNormalisation(testCase)

    ph = [0 90 180];
    th = [45 45 45];

    Eph = [1 2 3];
    Eth = [4 5 6];

    SF = SphereField( ...
        ph, th, Eph, Eth, 10e9);

    verifySize(testCase, SF.ph, [3 1]);
    verifySize(testCase, SF.th, [3 1]);

    verifySize(testCase, SF.Eph, [3 1]);
    verifySize(testCase, SF.Eth, [3 1]);

end


%% Multiple frequencies

function testMultipleFrequencies(testCase)

    ph = [0; 90; 180];
    th = [45; 45; 45];

    freq = [1e9 2e9];

    Eph = [ ...
        1 10
        2 20
        3 30];

    Eth = [ ...
        4 40
        5 50
        6 60];

    SF = SphereField( ...
        ph, th, Eph, Eth, freq);

    verifyEqual(testCase, SF.Np, 3);
    verifyEqual(testCase, SF.Nf, 2);

    verifySize(testCase, SF.Eph, [3 2]);
    verifySize(testCase, SF.Eth, [3 2]);

    verifyEqual(testCase, SF.Eph, Eph);
    verifyEqual(testCase, SF.Eth, Eth);

end


%% Optional radial field

function testNoRadialField(testCase)

    SF = makeSimpleField();

    verifyFalse(testCase, SF.hasEr);
    verifyEmpty(testCase, SF.Er);

end


function testRadialField(testCase)

    ph = [0; 90; 180];
    th = [45; 45; 45];

    Eph = ones(3,1);
    Eth = 2*ones(3,1);
    Er  = 3*ones(3,1);

    SF = SphereField( ...
        ph, th, Eph, Eth, 10e9, ...
        Er = Er);

    verifyTrue(testCase, SF.hasEr);
    verifyEqual(testCase, SF.Er, Er);

end


%% Radius

function testRadius(testCase)

    SF = SphereField( ...
        0, 0, 1, 2, 10e9, ...
        r = 12.5);

    verifyEqual(testCase, SF.r, 12.5);

end


%% Radiated power

function testPrad(testCase)

    freq = [1e9 2e9 3e9];

    Eph = ones(2,3);
    Eth = ones(2,3);

    SF = SphereField( ...
        [0;90], ...
        [45;45], ...
        Eph, Eth, freq, ...
        Prad = [10 20 30]);

    verifyEqual(testCase, ...
        SF.Prad, ...
        [10 20 30]);

end


function testDefaultPradEmpty(testCase)

    SF = makeSimpleField();

    verifyEmpty(testCase, SF.Prad);

end


%% Structured grids

function testStructuredGrid(testCase)

    phVec = [0 90 180 270];
    thVec = [0 30 60];

    [PH,TH] = meshgrid(phVec,thVec);

    Np = numel(PH);

    SF = SphereField( ...
        PH, TH, ...
        ones(Np,1), ...
        2*ones(Np,1), ...
        10e9);

    verifyTrue(testCase, SF.isStructured);

    verifyEqual(testCase, SF.Nph, 4);
    verifyEqual(testCase, SF.Nth, 3);

    verifyEqual(testCase, ...
        SF.gridInfo.phVec, ...
        phVec);

    verifyEqual(testCase, ...
        SF.gridInfo.thVec, ...
        thVec);

end


function testStructuredGridOrdering(testCase)

    phVec = [0 90 180];
    thVec = [10 20];

    [PH,TH] = meshgrid(phVec,thVec);

    A = [ ...
        1 2 3
        4 5 6];

    SF = SphereField( ...
        PH, TH, ...
        A(:), A(:), ...
        10e9);

    Arec = zeros(SF.Nth,SF.Nph);

    Arec(SF.gridInfo.gridIndex) = SF.Eth(:,1);

    verifyEqual(testCase, Arec, A);

end


function testStructuredGridRandomOrdering(testCase)

    phVec = [0 90 180];
    thVec = [10 20];

    [PH,TH] = meshgrid(phVec,thVec);

    A = [ ...
        1 2 3
        4 5 6];

    ind = [4 1 6 3 2 5];

    SF = SphereField( ...
        PH(ind), ...
        TH(ind), ...
        A(ind), ...
        A(ind), ...
        10e9);

    verifyTrue(testCase, SF.isStructured);

    Arec = zeros(SF.Nth,SF.Nph);

    Arec(SF.gridInfo.gridIndex) = SF.Eth(:,1);

    verifyEqual(testCase, Arec, A);

end


function testUnstructuredGrid(testCase)

    ph = [0; 30; 60; 80];
    th = [10; 20; 30; 50];

    SF = SphereField( ...
        ph, th, ...
        ones(4,1), ...
        ones(4,1), ...
        10e9);

    verifyFalse(testCase, SF.isStructured);

    verifyEmpty(testCase, SF.Nph);
    verifyEmpty(testCase, SF.Nth);

end

%% Sphere coverage

function testFullSphereWithRepeatedPhiEndpoint(testCase)

    SF = makeGridField(0:10:360, 0:5:180);

    verifyTrue(testCase, SF.isStructured);

    verifyTrue(testCase, SF.hasPhiPeriodicEndpoint);
    verifyTrue(testCase, SF.coversFullPhi);
    verifyTrue(testCase, SF.coversFullTheta);
    verifyTrue(testCase, SF.isFullSphere);

end


function testFullSphereWithoutRepeatedPhiEndpoint(testCase)

    SF = makeGridField(0:10:350, 0:5:180);

    verifyTrue(testCase, SF.isStructured);

    verifyFalse(testCase, SF.hasPhiPeriodicEndpoint);
    verifyTrue(testCase, SF.coversFullPhi);
    verifyTrue(testCase, SF.coversFullTheta);
    verifyTrue(testCase, SF.isFullSphere);

end


function testIncompletePhiCoverage(testCase)

    SF = makeGridField(0:10:180, 0:5:180);

    verifyFalse(testCase, SF.coversFullPhi);
    verifyTrue(testCase, SF.coversFullTheta);

    verifyFalse(testCase, SF.isFullSphere);

end


function testIncompleteThetaCoverage(testCase)

    SF = makeGridField(0:10:360, 0:5:90);

    verifyTrue(testCase, SF.coversFullPhi);
    verifyFalse(testCase, SF.coversFullTheta);

    verifyFalse(testCase, SF.isFullSphere);

end


function testShiftedFullPhiGrid(testCase)

    SF = makeGridField(-180:10:180, 0:5:180);

    verifyTrue(testCase, SF.hasPhiPeriodicEndpoint);
    verifyTrue(testCase, SF.coversFullPhi);
    verifyTrue(testCase, SF.isFullSphere);

end


function testShiftedFullPhiWithoutEndpoint(testCase)

    SF = makeGridField(-180:10:170, 0:5:180);

    verifyFalse(testCase, SF.hasPhiPeriodicEndpoint);
    verifyTrue(testCase, SF.coversFullPhi);
    verifyTrue(testCase, SF.isFullSphere);

end


function testPhiSector(testCase)

    SF = makeGridField(-45:5:45, 0:5:180);

    verifyFalse(testCase, SF.coversFullPhi);
    verifyFalse(testCase, SF.isFullSphere);

end

function testRedundantPhiCoverage(testCase)

    SF = makeGridField(0:10:720, 0:5:180);

    verifyTrue(testCase, SF.isStructured);

    verifyTrue(testCase, SF.coversFullPhi);
    verifyTrue(testCase, SF.coversFullTheta);
    verifyTrue(testCase, SF.isFullSphere);

    verifyTrue(testCase, SF.hasRedundantPhi);

end

function testNoRedundantPhiForSingleRevolution(testCase)

    SF = makeGridField(0:10:360, 0:5:180);

    verifyTrue(testCase, SF.isFullSphere);
    verifyFalse(testCase, SF.hasRedundantPhi);

end

function testBothPolesPresent(testCase)

    SF = makeGridField(0:10:360, 0:5:180);

    verifyTrue(testCase, SF.hasNorthPole);
    verifyTrue(testCase, SF.hasSouthPole);

end


function testNoPolesPresent(testCase)

    SF = makeGridField(0:10:360, 5:5:175);

    verifyFalse(testCase, SF.hasNorthPole);
    verifyFalse(testCase, SF.hasSouthPole);

end


function testNorthPoleOnly(testCase)

    SF = makeGridField(0:10:360, 0:5:175);

    verifyTrue(testCase, SF.hasNorthPole);
    verifyFalse(testCase, SF.hasSouthPole);

end


function testSouthPoleOnly(testCase)

    SF = makeGridField(0:10:360, 5:5:180);

    verifyFalse(testCase, SF.hasNorthPole);
    verifyTrue(testCase, SF.hasSouthPole);

end


function testNorthPoleIndices(testCase)

    phVec = 0:10:360;
    thVec = 0:5:180;

    SF = makeGridField(phVec, thVec);

    idx = SF.gridInfo.northPoleIndex;

    verifyEqual(testCase, numel(idx), numel(phVec));

    verifyEqual(testCase, ...
        SF.th(idx), ...
        zeros(numel(phVec),1));

end


function testSouthPoleIndices(testCase)

    phVec = 0:10:360;
    thVec = 0:5:180;

    SF = makeGridField(phVec, thVec);

    idx = SF.gridInfo.southPoleIndex;

    verifyEqual(testCase, numel(idx), numel(phVec));

    verifyEqual(testCase, ...
        SF.th(idx), ...
        180*ones(numel(phVec),1));

end


function testPoleIndicesMatchPhiSamples(testCase)

    phVec = -180:20:180;
    thVec = 0:10:180;

    SF = makeGridField(phVec, thVec);

    idxN = SF.gridInfo.northPoleIndex;
    idxS = SF.gridInfo.southPoleIndex;

    verifyEqual(testCase, ...
        sort(SF.ph(idxN)), ...
        phVec(:));

    verifyEqual(testCase, ...
        sort(SF.ph(idxS)), ...
        phVec(:));

end


function testNoPoleIndicesWhenAbsent(testCase)

    SF = makeGridField(0:10:350, 5:5:175);

    verifyEmpty(testCase, SF.gridInfo.northPoleIndex);
    verifyEmpty(testCase, SF.gridInfo.southPoleIndex);

end

function SF = makeGridField(phVec, thVec)

    [PH,TH] = meshgrid(phVec,thVec);

    Eph = zeros(numel(PH),1);
    Eth = sind(TH(:));

    SF = SphereField( ...
        PH, TH, ...
        Eph, Eth, ...
        1e9);

end



%% Input validation

function testCoordinateMismatch(testCase)

    verifyError(testCase, ...
        @() SphereField( ...
        [0 10 20], ...
        [0 10], ...
        ones(3,1), ...
        ones(3,1), ...
        10e9), ...
        'SphereField:CoordinateSizeMismatch');

end


function testFieldSizeMismatch(testCase)

    verifyError(testCase, ...
        @() SphereField( ...
        [0;10;20], ...
        [0;10;20], ...
        ones(2,1), ...
        ones(3,1), ...
        10e9), ...
        'SphereField:FieldSizeMismatch');

end


function testMultiFrequencyFieldSizeMismatch(testCase)

    verifyError(testCase, ...
        @() SphereField( ...
        [0;10;20], ...
        [0;10;20], ...
        ones(3,2), ...
        ones(3,2), ...
        [1e9 2e9 3e9]), ...
        'SphereField:FieldSizeMismatch');

end


function testPradSizeMismatch(testCase)

    verifyError(testCase, ...
        @() SphereField( ...
        [0;10], ...
        [0;10], ...
        ones(2,2), ...
        ones(2,2), ...
        [1e9 2e9], ...
        Prad = [1 2 3]), ...
        'SphereField:PradSizeMismatch');

end


%% Metadata / superclass interaction

function testMetadata(testCase)

    metadata.name = "Synthetic pattern";

    provenance.solver = "Analytical";

    SF = SphereField( ...
        0, 0, 1, 2, 10e9, ...
        Metadata = metadata, ...
        Provenance = provenance);

    verifyEqual(testCase, ...
        SF.metadata, metadata);

    verifyEqual(testCase, ...
        SF.provenance, provenance);

end


%% Helper

function SF = makeSimpleField()

    SF = SphereField( ...
        [0; 90; 180], ...
        [45; 45; 45], ...
        [1; 2; 3], ...
        [4; 5; 6], ...
        10e9);

end

%% Cartesian representation

function testCartesianDirections(testCase)

    ph = [0; 90; 180; 270];
    th = 90*ones(4,1);

    SF = SphereField( ...
        ph, th, ...
        zeros(4,1), ...
        ones(4,1), ...
        1e9);

    [rHat,~,~] = SF.getCartesianDirections();

    expected = [ ...
         1  0 -1  0
         0  1  0 -1
         0  0  0  0];

    verifyEqual(testCase, ...
        rHat, expected, ...
        AbsTol = 1e-14);

end


function testCartesianUnitVectors(testCase)

    SF = makeGridField(0:30:330, 0:15:180);

    [rHat,thHat,phHat] = SF.getCartesianDirections();

    verifyEqual(testCase, ...
        sqrt(sum(rHat.^2,1)), ...
        ones(1,SF.Np), ...
        AbsTol = 1e-14);

    verifyEqual(testCase, ...
        sqrt(sum(thHat.^2,1)), ...
        ones(1,SF.Np), ...
        AbsTol = 1e-14);

    verifyEqual(testCase, ...
        sqrt(sum(phHat.^2,1)), ...
        ones(1,SF.Np), ...
        AbsTol = 1e-14);

end


function testCartesianBasisOrthogonal(testCase)

    SF = makeGridField(0:30:330, 0:15:180);

    [rHat,thHat,phHat] = SF.getCartesianDirections();

    verifyEqual(testCase, ...
        sum(rHat.*thHat,1), ...
        zeros(1,SF.Np), ...
        AbsTol = 1e-14);

    verifyEqual(testCase, ...
        sum(rHat.*phHat,1), ...
        zeros(1,SF.Np), ...
        AbsTol = 1e-14);

    verifyEqual(testCase, ...
        sum(thHat.*phHat,1), ...
        zeros(1,SF.Np), ...
        AbsTol = 1e-14);

end


function testCartesianThetaField(testCase)

    SF = SphereField( ...
        0, 90, ...
        0, ...
        2, ...
        1e9);

    E = SF.getCartesianE();

    FFfact = exp(-1i*SF.wavenumber*SF.r)/SF.r;

    verifyEqual(testCase, ...
        E(:,1,1), ...
        [0;0;-2*FFfact], ...
        AbsTol = 1e-14);

end


function testCartesianPhiField(testCase)

    SF = SphereField( ...
        0, 90, ...
        3, ...
        0, ...
        1e9);

    E = SF.getCartesianE();

    FFfact = exp(-1i*SF.wavenumber*SF.r)/SF.r;

    verifyEqual(testCase, ...
        E(:,1,1), ...
        [0;3*FFfact;0], ...
        AbsTol = 1e-14);

end

function testCartesianRadialField(testCase)

    SF = SphereField( ...
        0, 90, ...
        0, 0, ...
        1e9, ...
        Er = 4);

    E = SF.getCartesianE();

    FFfact = exp(-1i*SF.wavenumber*SF.r)/SF.r;

    verifyEqual(testCase, ...
        E(:,1,1), ...
        [4*FFfact;0;0], ...
        AbsTol = 1e-14);

end

function testGetEfieldFarFieldFactor(testCase)

    SF = SphereField( ...
        0, 90, ...
        3, 2, ...
        2e9, ...
        Er = 4, ...
        r = 5);

    [Eph,Eth,Er] = SF.getEfield();

    FFfact = exp(-1i*SF.wavenumber*SF.r)/SF.r;

    verifyEqual(testCase, Eph, 3*FFfact, AbsTol = 1e-14);
    verifyEqual(testCase, Eth, 2*FFfact, AbsTol = 1e-14);
    verifyEqual(testCase, Er,  4*FFfact, AbsTol = 1e-14);

end


function testCartesianMultipleFrequencies(testCase)

    SF = SphereField( ...
        [0;90], ...
        [90;90], ...
        zeros(2,2), ...
        [1 10; 2 20], ...
        [1e9 2e9]);

    E = SF.getCartesianE();

    verifySize(testCase, E, [3 2 2]);

end

%% Pole consistency

function testPoleConsistencyZeroAtPoles(testCase)

    SF = makeGridField(0:10:360, 0:5:180);

    [tf,details] = SF.checkPoleConsistency();

    verifyTrue(testCase, tf);
    verifyTrue(testCase, details.northConsistent);
    verifyTrue(testCase, details.southConsistent);

end


function testNorthPoleCartesianConsistency(testCase)

    ph = (0:30:360).';
    th = zeros(size(ph));

    % Same physical Cartesian field:
    %
    %   E = x_hat
    %
    % At theta = 0:
    %
    %   Eth = cos(phi)
    %   Eph = -sin(phi)

    Eth = cosd(ph);
    Eph = -sind(ph);

    SF = SphereField( ...
        ph, th, ...
        Eph, Eth, ...
        1e9);

    [tf,details] = SF.checkPoleConsistency();

    verifyTrue(testCase, tf);
    verifyTrue(testCase, details.northConsistent);

end


function testSouthPoleCartesianConsistency(testCase)

    ph = (0:30:360).';
    th = 180*ones(size(ph));

    % Same physical Cartesian field:
    %
    %   E = x_hat
    %
    % At theta = 180:
    %
    %   thHat = [-cos(phi); -sin(phi); 0]
    %   phHat = [-sin(phi);  cos(phi); 0]
    %
    % Therefore:
    %
    %   Eth = -cos(phi)
    %   Eph = -sin(phi)

    Eth = -cosd(ph);
    Eph = -sind(ph);

    SF = SphereField( ...
        ph, th, ...
        Eph, Eth, ...
        1e9);

    [tf,details] = SF.checkPoleConsistency();

    verifyTrue(testCase, tf);
    verifyTrue(testCase, details.southConsistent);

end


function testNorthPoleInconsistency(testCase)

    ph = (0:30:360).';
    th = zeros(size(ph));

    % Constant nonzero Eth is NOT physically constant at the pole,
    % because thHat rotates with phi.

    Eth = ones(size(ph));
    Eph = zeros(size(ph));

    SF = SphereField( ...
        ph, th, ...
        Eph, Eth, ...
        1e9);

    [tf,details] = SF.checkPoleConsistency();

    verifyFalse(testCase, tf);
    verifyFalse(testCase, details.northConsistent);

end


function testPoleConsistencyMultipleFrequencies(testCase)

    ph = (0:30:360).';
    th = zeros(size(ph));

    Eth1 = cosd(ph);
    Eph1 = -sind(ph);

    Eth = [Eth1, 2*Eth1];
    Eph = [Eph1, 2*Eph1];

    SF = SphereField( ...
        ph, th, ...
        Eph, Eth, ...
        [1e9 2e9]);

    [tf,details] = SF.checkPoleConsistency();

    verifyTrue(testCase, tf);

    verifySize(testCase, ...
        details.northMaxError, ...
        [1 2]);

end


function testDefaultFieldPoleConsistency(testCase)

    SF = SphereField;

    verifyTrue(testCase, SF.isFullSphere);

    tf = SF.checkPoleConsistency();

    verifyTrue(testCase, tf);

end

%% getGridView

function testStoredGridView(testCase)

    phVec = 0:10:360;
    thVec = 0:5:180;

    SF = makeGridField(phVec,thVec);

    G = SF.getGridView("stored");

    verifyTrue(testCase, G.isStructured);

    verifyEqual(testCase, G.phVec, phVec);
    verifyEqual(testCase, G.thVec, thVec);

    verifySize(testCase, ...
        G.Eth, ...
        [numel(thVec), numel(phVec), 1]);

end

function testSinglePeriodAlreadyClosed(testCase)

    phVec = 0:10:360;

    SF = makeGridField(phVec,0:10:180);

    G = SF.getGridView("singlePeriod");

    verifyEqual(testCase, ...
        G.phVec, ...
        phVec);

    verifyEqual(testCase, ...
        size(G.Eth,2), ...
        numel(phVec));

end

function testSinglePeriodAppendsEndpoint(testCase)

    SF = makeGridField(0:10:350,0:10:180);

    G = SF.getGridView("singlePeriod");

    verifyEqual(testCase, ...
        G.phVec, ...
        0:10:360);

    verifyEqual(testCase, ...
        G.Eth(:,end,:), ...
        G.Eth(:,1,:), ...
        AbsTol = 1e-14);

    verifyEqual(testCase, ...
        G.Eph(:,end,:), ...
        G.Eph(:,1,:), ...
        AbsTol = 1e-14);

end

function testSinglePeriodRemovesExtraRevolution(testCase)

    SF = makeGridField(0:10:720,0:10:180);

    G = SF.getGridView("singlePeriod");

    verifyEqual(testCase, ...
        G.phVec, ...
        0:10:360);

    verifyEqual(testCase, ...
        size(G.Eth,2), ...
        37);

end

function testSinglePeriodShiftedPhi(testCase)

    SF = makeGridField(-180:10:180,0:10:180);

    G = SF.getGridView("singlePeriod");

    verifyEqual(testCase, ...
        G.phVec, ...
        -180:10:180);

end

function testSinglePeriodShiftedWithoutEndpoint(testCase)

    SF = makeGridField(-180:10:170,0:10:180);

    G = SF.getGridView("singlePeriod");

    verifyEqual(testCase, ...
        G.phVec, ...
        -180:10:180);

end

function testSinglePeriodRejectsPartialPhi(testCase)

    SF = makeGridField(0:10:180,0:10:180);

    verifyError(testCase, ...
        @() SF.getGridView("singlePeriod"), ...
        'SphereField:SinglePeriodRequiresFullPhi');

end

function testStoredGridViewRandomSampleOrdering(testCase)

    phVec = [0 90 180 270 360];
    thVec = [0 45 90];

    [PH,TH] = meshgrid(phVec,thVec);

    A = reshape(1:numel(PH),size(PH));

    ind = randperm(numel(PH));

    SF = SphereField( ...
        PH(ind), ...
        TH(ind), ...
        zeros(numel(PH),1), ...
        A(ind), ...
        1e9);

    G = SF.getGridView("stored");

    verifyEqual(testCase, ...
        G.Eth(:,:,1), ...
        A);

end

%% Power integration

function testPowerIntegrationSinTheta(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    E0 = 3;
    r = 2;

    Eth = E0*sind(TH);
    Eph = zeros(size(Eth));

    SF = SphereField( ...
        PH, TH, ...
        Eph(:), Eth(:), ...
        1e9, ...
        r = r);

    P = SF.integratePower();

    Pexact = 4*pi*E0^2 / (3*SF.eta0);

    verifyEqual(testCase, ...
        P, Pexact, ...
        RelTol = 1e-4);

end

function testPowerIntegrationWithoutPhiEndpoint(testCase)

    phVec = 0:2:358;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    E0 = 3;
    r = 2;

    Eth = E0*sind(TH);
    Eph = zeros(size(Eth));

    SF = SphereField( ...
        PH, TH, ...
        Eph(:), Eth(:), ...
        1e9, ...
        r = r);

    P = SF.integratePower();

    Pexact = 4*pi*E0^2 / (3*SF.eta0);

    verifyEqual(testCase, ...
        P, Pexact, ...
        RelTol = 1e-4);

end

function testPowerIntegrationRedundantPhi(testCase)

    phVec = 0:2:720;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    E0 = 3;
    r = 2;

    Eth = E0*sind(TH);
    Eph = zeros(size(Eth));

    SF = SphereField( ...
        PH, TH, ...
        Eph(:), Eth(:), ...
        1e9, ...
        r = r);

    P = SF.integratePower();

    Pexact = 4*pi*E0^2 / (3*SF.eta0);

    verifyEqual(testCase, ...
        P, Pexact, ...
        RelTol = 1e-4);

end

function testPowerIntegrationMultipleFrequencies(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    r = 1.5;

    E1 = sind(TH);
    E2 = 2*sind(TH);
    E3 = 3*sind(TH);

    Eth = [E1(:), E2(:), E3(:)];
    Eph = zeros(size(Eth));

    SF = SphereField( ...
        PH, TH, ...
        Eph, Eth, ...
        [1e9 2e9 3e9], ...
        r = r);

    P = SF.integratePower();

    P1 = 4*pi/(3*SF.eta0);

    verifyEqual(testCase, ...
        P, ...
        P1*[1 4 9], ...
        RelTol = 1e-4);

end

function testPowerIntegrationRejectsPartialSphere(testCase)

    SF = makeGridField(0:10:180, 0:5:180);

    verifyError(testCase, ...
        @() SF.integratePower(), ...
        'SphereField:PowerIntegrationRequiresFullSphere');

end

function testPowerIntegrationIndependentOfRadius(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    Eth = sind(TH);
    Eph = zeros(size(Eth));

    SF1 = SphereField( ...
        PH,TH, ...
        Eph(:),Eth(:), ...
        1e9, ...
        r = 1);

    SF2 = SphereField( ...
        PH,TH, ...
        Eph(:),Eth(:), ...
        1e9, ...
        r = 10);

    P1 = SF1.integratePower();
    P2 = SF2.integratePower();

    verifyEqual(testCase, ...
        P1, P2, ...
        RelTol = 1e-14);

end

function testRadiationIntensityIndependentOfRadius(testCase)

    ph = [0;90;180];
    th = [45;45;45];

    Eph = [1;2;3];
    Eth = [4;5;6];

    SF1 = SphereField(ph,th,Eph,Eth,1e9,r=1);
    SF2 = SphereField(ph,th,Eph,Eth,1e9,r=10);

    U1 = SF1.getU();
    U2 = SF2.getU();

    verifyEqual(testCase,U1,U2,RelTol=1e-14);

end

function testDirectivitySinTheta(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    Eth = sind(TH);
    Eph = zeros(size(Eth));

    SF = SphereField( ...
        PH,TH, ...
        Eph(:),Eth(:), ...
        1e9);

    D = SF.getDirectivity( ...
        PowerSource="integrated");

    verifyEqual(testCase, ...
        max(D), ...
        1.5, ...
        RelTol=1e-4);

end