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

%% Radiation efficiency and gain

function testDefaultEtaRad(testCase)

    SF = makeSimpleField();

    verifyEqual(testCase, SF.etaRad, 1);

end


function testEtaRadConstructor(testCase)

    SF = SphereField( ...
        0, 90, ...
        0, 1, ...
        1e9, ...
        etaRad = 0.75);

    verifyEqual(testCase, SF.etaRad, 0.75);

end


function testEtaRadMultipleFrequencies(testCase)

    ph = [0; 90];
    th = [45; 45];

    freq = [1e9 2e9 3e9];

    Eph = zeros(2,3);
    Eth = ones(2,3);

    etaRad = [0.5 0.7 0.9];

    SF = SphereField( ...
        ph, th, ...
        Eph, Eth, ...
        freq, ...
        etaRad = etaRad);

    verifyEqual(testCase, SF.etaRad, etaRad);

end


function testEtaRadSizeMismatch(testCase)

    ph = [0; 90];
    th = [45; 45];

    Eph = ones(2,2);
    Eth = ones(2,2);

    verifyError(testCase, ...
        @() SphereField( ...
            ph, th, ...
            Eph, Eth, ...
            [1e9 2e9], ...
            etaRad = [0.5 0.6 0.7]), ...
        'SphereField:EtaRadSizeMismatch');

end


function testEtaRadGreaterThanOneRejected(testCase)

    verifyWarning(testCase, ...
        @() SphereField( ...
            0, 90, ...
            0, 1, ...
            1e9, ...
            etaRad = 1.1), ...
        'SphereField:InvalidEtaRad');

end


function testNegativeEtaRadRejected(testCase)

    verifyError(testCase, ...
        @() SphereField( ...
            0, 90, ...
            0, 1, ...
            1e9, ...
            etaRad = -0.1), ...
        'SphereField:InvalidEtaRad');

end


function testGainUnityEfficiency(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    Eth = sind(TH);
    Eph = zeros(size(Eth));

    SF = SphereField( ...
        PH, TH, ...
        Eph(:), Eth(:), ...
        1e9, ...
        etaRad = 1);

    D = SF.getDirectivity( ...
        PowerSource = "integrated");

    G = SF.getGain( ...
        PowerSource = "integrated");

    verifyEqual(testCase, ...
        G, D, ...
        RelTol = 1e-14);

end


function testGainRadiationEfficiency(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    Eth = sind(TH);
    Eph = zeros(size(Eth));

    etaRad = 0.6;

    SF = SphereField( ...
        PH, TH, ...
        Eph(:), Eth(:), ...
        1e9, ...
        etaRad = etaRad);

    D = SF.getDirectivity( ...
        PowerSource = "integrated");

    G = SF.getGain( ...
        PowerSource = "integrated");

    verifyEqual(testCase, ...
        G, etaRad*D, ...
        RelTol = 1e-14);

end


function testGainMultipleFrequencies(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    E1 = sind(TH);

    Eth = [ ...
        E1(:), ...
        2*E1(:), ...
        3*E1(:)];

    Eph = zeros(size(Eth));

    etaRad = [0.4 0.6 0.8];

    SF = SphereField( ...
        PH, TH, ...
        Eph, Eth, ...
        [1e9 2e9 3e9], ...
        etaRad = etaRad);

    D = SF.getDirectivity( ...
        PowerSource = "integrated");

    G = SF.getGain( ...
        PowerSource = "integrated");

    expected = D .* reshape(etaRad,1,[]);

    verifyEqual(testCase, ...
        G, expected, ...
        RelTol = 1e-14);

end


function testMaxGainSinTheta(testCase)

    phVec = 0:2:360;
    thVec = 0:2:180;

    [PH,TH] = meshgrid(phVec,thVec);

    Eth = sind(TH);
    Eph = zeros(size(Eth));

    etaRad = 0.8;

    SF = SphereField( ...
        PH, TH, ...
        Eph(:), Eth(:), ...
        1e9, ...
        etaRad = etaRad);

    G = SF.getGain( ...
        PowerSource = "integrated");

    % sin(theta) pattern has Dmax = 1.5.
    GmaxExact = etaRad * 1.5;

    verifyEqual(testCase, ...
        max(G), ...
        GmaxExact, ...
        RelTol = 1e-4);

end

%% Linear field basis

function testSphericalBasis(testCase)

    SF = SphereField( ...
        [0;45], ...
        [30;60], ...
        [1;2], ...      % Eph
        [3;4], ...      % Eth
        1e9);

    [E1,E2,E3,info] = ...
        SF.getLinearBasis(Basis = "spherical");

    verifyEqual(testCase,E1,SF.Eth);
    verifyEqual(testCase,E2,SF.Eph);
    verifyEmpty(testCase,E3);

    verifyEqual(testCase, ...
        info.componentNames, ...
        ["Eth","Eph"]);

end

function testLudwig1BasisAtZAxis(testCase)

    SF = SphereField( ...
        0, 0, ...
        3, ...      % Eph
        2, ...      % Eth
        1e9);

    [Ex,Ey,Ez,info] = ...
        SF.getLinearBasis(Basis = "Ludwig1");

    verifyEqual(testCase,Ex,2,AbsTol=1e-14);
    verifyEqual(testCase,Ey,3,AbsTol=1e-14);
    verifyEqual(testCase,Ez,0,AbsTol=1e-14);

    verifyEqual(testCase, ...
        info.componentNames, ...
        ["Ex","Ey","Ez"]);

end

function testLudwig1BasisAtXAxis(testCase)

    SF = SphereField( ...
        0, 90, ...
        3, ...      % Eph
        2, ...      % Eth
        1e9);

    [Ex,Ey,Ez] = ...
        SF.getLinearBasis(Basis = "Ludwig1");

    verifyEqual(testCase,Ex,0,AbsTol=1e-14);
    verifyEqual(testCase,Ey,3,AbsTol=1e-14);
    verifyEqual(testCase,Ez,-2,AbsTol=1e-14);

end

function testLudwig1MatchesCartesianPattern(testCase)

    ph = [0;30;90;170];
    th = [10;40;90;130];

    Eph = [1;2;3;4] + 1i*[2;1;4;3];
    Eth = [5;6;7;8] - 1i*[1;3;2;4];

    SF = SphereField( ...
        ph,th,Eph,Eth,1e9);

    [Ex,Ey,Ez] = ...
        SF.getLinearBasis(Basis = "Ludwig1");

    Ecart = SF.getCartesianE();

    FFfact = ...
        exp(-1i*SF.wavenumber*SF.r)/SF.r;

    Epattern = Ecart(:,:,1) ./ FFfact;

    verifyEqual(testCase, ...
        [Ex.';Ey.';Ez.'], ...
        Epattern, ...
        AbsTol=1e-13);

end

function testLudwig3PhiZero(testCase)

    SF = SphereField( ...
        0, 45, ...
        3, ...       % Eph
        2, ...       % Eth
        1e9);

    [Eh,Ev] = ...
        SF.getLinearBasis(Basis = "Ludwig3");

    verifyEqual(testCase,Eh,2,AbsTol=1e-14);
    verifyEqual(testCase,Ev,3,AbsTol=1e-14);

end

function testLudwig3Phi90(testCase)

    SF = SphereField( ...
        90, 45, ...
        3, ...
        2, ...
        1e9);

    [Eh,Ev] = ...
        SF.getLinearBasis(Basis = "Ludwig3");

    verifyEqual(testCase,Eh,-3,AbsTol=1e-14);
    verifyEqual(testCase,Ev, 2,AbsTol=1e-14);

end

function testLudwig3PreservesMagnitude(testCase)

    SF = makeComplexField();

    [Eh,Ev] = ...
        SF.getLinearBasis(Basis = "Ludwig3");

    magSph = ...
        abs(SF.Eth).^2 + abs(SF.Eph).^2;

    magL3 = ...
        abs(Eh).^2 + abs(Ev).^2;

    verifyEqual(testCase, ...
        magL3,magSph, ...
        AbsTol=1e-13);

end

function testLudwig2AEPhiZero(testCase)

    SF = SphereField( ...
        0, 40, ...
        3, ...
        2, ...
        1e9);

    [Eaz,Eel] = ...
        SF.getLinearBasis(Basis = "Ludwig2AE");

    verifyEqual(testCase,Eaz,2,AbsTol=1e-14);
    verifyEqual(testCase,Eel,3,AbsTol=1e-14);

end

function testLudwig2AESingularity(testCase)

    SF = SphereField( ...
        90, 90, ...
        3, ...       % Eph
        2, ...       % Eth
        1e9);

    [Eaz,Eel] = ...
        SF.getLinearBasis(Basis = "Ludwig2AE");

    verifyEqual(testCase,Eaz, 3,AbsTol=1e-14);
    verifyEqual(testCase,Eel,-2,AbsTol=1e-14);

    verifyTrue(testCase,isfinite(Eaz));
    verifyTrue(testCase,isfinite(Eel));

end

function testLudwig2EAPhi90(testCase)

    SF = SphereField( ...
        90, 40, ...
        3, ...
        2, ...
        1e9);

    [Eal,Eep] = ...
        SF.getLinearBasis(Basis = "Ludwig2EA");

    verifyEqual(testCase,Eal,-3,AbsTol=1e-14);
    verifyEqual(testCase,Eep, 2,AbsTol=1e-14);

end

function testLudwig2EASingularity(testCase)

    SF = SphereField( ...
        0, 90, ...
        3, ...
        2, ...
        1e9);

    [Eal,Eep] = ...
        SF.getLinearBasis(Basis = "Ludwig2EA");

    verifyEqual(testCase,Eal,-2,AbsTol=1e-14);
    verifyEqual(testCase,Eep,-3,AbsTol=1e-14);

    verifyTrue(testCase,isfinite(Eal));
    verifyTrue(testCase,isfinite(Eep));

end

function testLudwig2PreservesMagnitude(testCase)

    SF = makeComplexField();

    mag0 = ...
        abs(SF.Eth).^2 + abs(SF.Eph).^2;

    [E1,E2] = ...
        SF.getLinearBasis(Basis = "Ludwig2AE");

    verifyEqual(testCase, ...
        abs(E1).^2 + abs(E2).^2, ...
        mag0, ...
        AbsTol=1e-12);

    [E1,E2] = ...
        SF.getLinearBasis(Basis = "Ludwig2EA");

    verifyEqual(testCase, ...
        abs(E1).^2 + abs(E2).^2, ...
        mag0, ...
        AbsTol=1e-12);

end

function SF = makeComplexField()

    ph = [10;35;70;120;210;300];
    th = [15;35;55;75;110;145];

    Eph = ...
        [1;2;3;4;5;6] ...
        + 1i*[3;1;4;2;6;5];

    Eth = ...
        [6;5;4;3;2;1] ...
        - 1i*[1;3;2;5;4;6];

    SF = SphereField( ...
        ph,th,Eph,Eth,1e9);

end

%% Polarization projection

function testGetFieldLinear(testCase)
    SF = makeComplexField();
    [E1,E2,E3,info] = SF.getField(Basis="spherical",Polarization="linear");
    verifyEqual(testCase,E1,SF.Eth);
    verifyEqual(testCase,E2,SF.Eph);
    verifyEmpty(testCase,E3);
    verifyEqual(testCase,info.polarization,"linear");
end

function testCircularPureRHCP(testCase)
    % IEEE RHCP: E2 lags E1 by 90 deg for exp(+jwt).
    SF = SphereField(0,45,-1i,1,1e9);
    [ER,EL] = SF.getField(Basis="spherical",Polarization="circular");
    verifyEqual(testCase,ER,sqrt(2),AbsTol=1e-14);
    verifyEqual(testCase,EL,0,AbsTol=1e-14);
end

function testCircularPureLHCP(testCase)
    SF = SphereField(0,45,1i,1,1e9);
    [ER,EL] = SF.getField(Basis="spherical",Polarization="circular");
    verifyEqual(testCase,ER,0,AbsTol=1e-14);
    verifyEqual(testCase,EL,sqrt(2),AbsTol=1e-14);
end

function testCircularComponentNames(testCase)
    SF = makeSimpleField();
    [~,~,~,info] = SF.getField(Polarization="circular");
    verifyEqual(testCase,info.componentNames,["RHCP","LHCP"]);
    verifyEqual(testCase,info.convention,"IEEE");
end

function testCircularPreservesPower(testCase)
    SF = makeComplexField();
    [ER,EL] = SF.getField(Basis="Ludwig3",Polarization="circular");
    [E1,E2] = SF.getLinearBasis(Basis="Ludwig3");
    verifyEqual(testCase,abs(ER).^2+abs(EL).^2,abs(E1).^2+abs(E2).^2,AbsTol=1e-12);
end

function testSlantZero(testCase)
    SF = makeComplexField();
    [Es,Eo] = SF.getField(Basis="spherical",Polarization="slant",SlantAngle=0);
    verifyEqual(testCase,Es,SF.Eth);
    verifyEqual(testCase,Eo,SF.Eph);
end

function testSlant90(testCase)
    SF = makeComplexField();
    [Es,Eo] = SF.getField(Basis="spherical",Polarization="slant",SlantAngle=90);
    verifyEqual(testCase,Es,SF.Eph,AbsTol=1e-14);
    verifyEqual(testCase,Eo,-SF.Eth,AbsTol=1e-14);
end

function testSlant45(testCase)
    SF = SphereField(0,45,0,1,1e9);
    [Es,Eo] = SF.getField(Polarization="slant",SlantAngle=45);
    verifyEqual(testCase,Es,1/sqrt(2),AbsTol=1e-14);
    verifyEqual(testCase,Eo,-1/sqrt(2),AbsTol=1e-14);
end

function testSlantPreservesPower(testCase)
    SF = makeComplexField();
    [Es,Eo] = SF.getField(Basis="Ludwig3",Polarization="slant",SlantAngle=37);
    [E1,E2] = SF.getLinearBasis(Basis="Ludwig3");
    verifyEqual(testCase,abs(Es).^2+abs(Eo).^2,abs(E1).^2+abs(E2).^2,AbsTol=1e-12);
end

function testLudwig1CircularRejected(testCase)
    SF = makeSimpleField();
    verifyError(testCase,@() SF.getField(Basis="Ludwig1",Polarization="circular"),'SphereField:PolarizationRequiresTransverseBasis');
end

function testLudwig1SlantRejected(testCase)
    SF = makeSimpleField();
    verifyError(testCase,@() SF.getField(Basis="Ludwig1",Polarization="slant"),'SphereField:PolarizationRequiresTransverseBasis');
end

%% Coordinate views

function testCoordinatesPhTh(testCase)
    SF = SphereField([0;45;180],[10;30;90],[1;1;1],[1;1;1],1e9);
    C = SF.getCoordinates(Coordinates="PhTh");

    verifyEqual(testCase,C.x,SF.ph);
    verifyEqual(testCase,C.y,SF.th);
    verifyEqual(testCase,C.xName,"phi");
    verifyEqual(testCase,C.yName,"theta");
    verifyEqual(testCase,C.xUnit,"deg");
    verifyEqual(testCase,C.yUnit,"deg");
end

function testCoordinatesDirCosAxes(testCase)
    SF = SphereField([0;90;0],[90;90;0],zeros(3,1),ones(3,1),1e9);
    C = SF.getCoordinates(Coordinates="DirCos");

    verifyEqual(testCase,C.x,[1;0;0],AbsTol=1e-14);
    verifyEqual(testCase,C.y,[0;1;0],AbsTol=1e-14);
    verifyEqual(testCase,C.z,[0;0;1],AbsTol=1e-14);
end

function testCoordinatesDirCosMatchesDirections(testCase)
    SF = makeComplexField();
    C = SF.getCoordinates(Coordinates="DirCos");
    rHat = SF.getCartesianDirections();

    verifyEqual(testCase,[C.x.';C.y.';C.z.'],rHat,AbsTol=1e-14);
end

function testCoordinatesAzElKnownPoints(testCase)
    SF = SphereField([0;90;180],[90;90;90],zeros(3,1),ones(3,1),1e9);
    C = SF.getCoordinates(Coordinates="AzEl");

    verifyEqual(testCase,C.x,[90;0;-90],AbsTol=1e-12);
    verifyEqual(testCase,C.y,[0;90;0],AbsTol=1e-12);
end

function testCoordinatesAzElPreservesDirection(testCase)
    SF = makeComplexField();
    C = SF.getCoordinates(Coordinates="AzEl");

    u = cosd(C.y).*sind(C.x);
    v = sind(C.y);
    w = cosd(C.y).*cosd(C.x);

    rHat = SF.getCartesianDirections();
    verifyEqual(testCase,[u.';v.';w.'],rHat,AbsTol=1e-13);
end

function testCoordinatesElAzPreservesDirection(testCase)
    SF = makeComplexField();
    C = SF.getCoordinates(Coordinates="ElAz");

    u = sind(C.y);
    v = cosd(C.y).*sind(C.x);
    w = cosd(C.y).*cosd(C.x);

    rHat = SF.getCartesianDirections();
    verifyEqual(testCase,[u.';v.';w.'],rHat,AbsTol=1e-13);
end

function testCoordinatesTrueView(testCase)
    SF = SphereField([0;90;180;270],30*ones(4,1),zeros(4,1),ones(4,1),1e9);
    C = SF.getCoordinates(Coordinates="TrueView");

    verifyEqual(testCase,C.x,[30;0;-30;0],AbsTol=1e-13);
    verifyEqual(testCase,C.y,[0;30;0;-30],AbsTol=1e-13);
end

function testCoordinatesTrueViewRadius(testCase)
    SF = makeComplexField();
    C = SF.getCoordinates(Coordinates="TrueView");

    thRecovered = hypot(C.x,C.y);
    verifyEqual(testCase,thRecovered,SF.th,AbsTol=1e-13);
end

function testCoordinatesArcSin(testCase)
    SF = makeComplexField();
    C = SF.getCoordinates(Coordinates="ArcSin");
    D = SF.getCoordinates(Coordinates="DirCos");

    verifyEqual(testCase,sind(C.x),D.x,AbsTol=1e-14);
    verifyEqual(testCase,sind(C.y),D.y,AbsTol=1e-14);
end

function testCoordinatesPreserveSampleCount(testCase)
    SF = makeComplexField();
    coordinateTypes = ["PhTh","AzEl","ElAz","DirCos","TrueView","ArcSin"];

    for cc = coordinateTypes
        C = SF.getCoordinates(Coordinates=cc);
        verifySize(testCase,C.x,[SF.Np 1]);
        verifySize(testCase,C.y,[SF.Np 1]);
    end
end

function testCoordinateStructuredTopology(testCase)
    SF = makeGridField(0:10:360,0:10:180);

    coordinateTypes = ["PhTh","AzEl","ElAz","DirCos","TrueView","ArcSin"];
    for cc = coordinateTypes
        C = SF.getCoordinates(Coordinates=cc);
        verifyTrue(testCase,C.hasStructuredTopology);
    end
end

function testCoordinateUnstructuredTopology(testCase)
    SF = makeComplexField();
    C = SF.getCoordinates(Coordinates="AzEl");
    verifyFalse(testCase,C.hasStructuredTopology);
end

function testGetDataPatternComplex(testCase)
    SF = makeComplexField();
    D = SF.getData(FieldType="pattern",Quantity="complex");
    verifyEqual(testCase,D.values(:,:,1),SF.Eth);
    verifyEqual(testCase,D.values(:,:,2),SF.Eph);
    verifyEqual(testCase,D.valueUnit,"V");
    verifyEqual(testCase,D.fieldType,"pattern");
end

function testGetDataPhysicalField(testCase)
    SF = makeComplexField();
    D = SF.getData(FieldType="field",Quantity="complex");
    [Eph,Eth] = SF.getEfield();
    verifyEqual(testCase,D.values(:,:,1),Eth);
    verifyEqual(testCase,D.values(:,:,2),Eph);
    verifyEqual(testCase,D.valueUnit,"V/m");
    verifyEqual(testCase,D.fieldType,"field");
end

function testGetDataPatternReal(testCase)
    SF = makeComplexField();
    D = SF.getData(Quantity="real");
    verifyEqual(testCase,D.values(:,:,1),real(SF.Eth));
    verifyEqual(testCase,D.values(:,:,2),real(SF.Eph));
end

function testGetDataPatternImag(testCase)
    SF = makeComplexField();
    D = SF.getData(Quantity="imag");
    verifyEqual(testCase,D.values(:,:,1),imag(SF.Eth));
    verifyEqual(testCase,D.values(:,:,2),imag(SF.Eph));
end

function testGetDataPhysicalReal(testCase)
    SF = makeComplexField();
    D = SF.getData(FieldType="field",Quantity="real");
    [Eph,Eth] = SF.getEfield();
    verifyEqual(testCase,D.values(:,:,1),real(Eth));
    verifyEqual(testCase,D.values(:,:,2),real(Eph));
end

function testGetDataPhysicalImag(testCase)
    SF = makeComplexField();
    D = SF.getData(FieldType="field",Quantity="imag");
    [Eph,Eth] = SF.getEfield();
    verifyEqual(testCase,D.values(:,:,1),imag(Eth));
    verifyEqual(testCase,D.values(:,:,2),imag(Eph));
end

function testGetDataPatternMagnitude(testCase)
    SF = makeComplexField();
    D = SF.getData(FieldType="pattern",Quantity="magnitude");
    verifyEqual(testCase,D.values(:,:,1),abs(SF.Eth));
    verifyEqual(testCase,D.values(:,:,2),abs(SF.Eph));
end

function testGetDataPhysicalMagnitude(testCase)
    SF = makeComplexField();
    D = SF.getData(FieldType="field",Quantity="magnitude");
    [Eph,Eth] = SF.getEfield();
    verifyEqual(testCase,D.values(:,:,1),abs(Eth));
    verifyEqual(testCase,D.values(:,:,2),abs(Eph));
end

function testGetDataPatternIndependentOfRadius(testCase)
    SF1 = SphereField(0,45,2+1i,3-2i,1e9,r=1);
    SF2 = SphereField(0,45,2+1i,3-2i,1e9,r=10);

    D1 = SF1.getData(FieldType="pattern");
    D2 = SF2.getData(FieldType="pattern");
    verifyEqual(testCase,D1.values,D2.values);
end

function testGetDataFieldScalesWithRadius(testCase)
    SF1 = SphereField(0,45,2+1i,3-2i,1e9,r=1);
    SF2 = SphereField(0,45,2+1i,3-2i,1e9,r=10);

    D1 = SF1.getData(FieldType="field",Quantity="magnitude");
    D2 = SF2.getData(FieldType="field",Quantity="magnitude");
    verifyEqual(testCase,D2.values,D1.values/10,RelTol=1e-14);
end

function testGetDataComponentSelection(testCase)
    SF = makeComplexField();
    D = SF.getData(Basis="Ludwig3",ComponentIndex=2);
    [~,Ev] = SF.getLinearBasis(Basis="Ludwig3");
    verifyEqual(testCase,D.values,Ev);
    verifyEqual(testCase,D.componentNames,"Ev");
end

function testGetDataCircularComponentSelection(testCase)
    SF = makeComplexField();
    D = SF.getData(Polarization="circular",ComponentIndex=1);
    [ER,~] = SF.getField(Polarization="circular");
    verifyEqual(testCase,D.values,ER);
    verifyEqual(testCase,D.componentNames,"RHCP");
end

function testGetDataUnknownComponent(testCase)
    SF = makeComplexField();
    verifyError(testCase,@() SF.getData(ComponentIndex=3),'SphereField:InvalidComponentIndex');
end

function testGetDataMagnitudeDB20(testCase)
    SF = makeComplexField();
    Dlin = SF.getData(Quantity="magnitude");
    DdB = SF.getData(Quantity="magnitude",Scale="dB20");
    verifyEqual(testCase,DdB.values,20*log10(Dlin.values),AbsTol=1e-13);
    verifyEqual(testCase,DdB.valueUnit,"dB");
end

function testGetDataDirectivityDB10(testCase)
    SF = makeGridField(0:10:360,0:5:180);
    Dlin = SF.getData(Quantity="directivity");
    DdB = SF.getData(Quantity="directivity",Scale="dB10");
    verifyEqual(testCase,DdB.values,10*log10(Dlin.values),AbsTol=1e-13);
end

function testGetDataInvalidDB20(testCase)
    SF = makeComplexField();
    verifyError(testCase,@() SF.getData(Quantity="real",Scale="dB20"),'SphereField:InvalidScale');
end

function testGetDataInvalidDB10(testCase)
    SF = makeComplexField();
    verifyError(testCase,@() SF.getData(Quantity="magnitude",Scale="dB10"),'SphereField:InvalidScale');
end

function testFullSpherePositive180(testCase)
    SF = makeGridField(0:10:360,0:5:180);
    G = SF.getGridView("fullSphere",PhiRange="positive",ThetaRange="180");
    verifyEqual(testCase,G.phVec,(0:10:360));
    verifyEqual(testCase,G.thVec,(0:5:180));
end

function testFullSphereSymmetric180(testCase)
    SF = makeGridField(0:10:360,0:5:180);
    G = SF.getGridView("fullSphere",PhiRange="symmetric",ThetaRange="180");
    verifyEqual(testCase,G.phVec,(-180:10:180));
    verifyEqual(testCase,G.thVec,(0:5:180));
end

function testFullSpherePositive360(testCase)
    SF = makeGridField(0:10:360,0:5:180);
    G = SF.getGridView("fullSphere",PhiRange="positive",ThetaRange="360");
    verifyEqual(testCase,G.phVec,(0:10:180));
    verifyEqual(testCase,G.thVec,(-180:5:180));
end

function testFullSphereSymmetric360(testCase)
    SF = makeGridField(0:10:360,0:5:180);
    G = SF.getGridView("fullSphere",PhiRange="symmetric",ThetaRange="360");
    verifyEqual(testCase,G.phVec,(-90:10:90));
    verifyEqual(testCase,G.thVec,(-180:5:180));
end

function testInterpolationGridPadding(testCase)
    SF = makeGridField(0:10:360,0:5:180);
    G = SF.getGridView("interpolation",PaddingSamples=2);
    verifyEqual(testCase,G.phVec,(-20:10:380));
    verifyEqual(testCase,G.thVec,(-10:5:190));
    verifyEqual(testCase,G.paddingSamples,2);
end

function testInterpolationFullSphere(testCase)
    SF=makeGridField(0:10:360,0:5:180);
    G=SF.getGridView("interpolation",PaddingSamples=2);

    verifyEqual(testCase,G.phVec,-20:10:380);
    verifyEqual(testCase,G.thVec,-10:5:190);
    verifyTrue(testCase,G.extendedPhi);
    verifyTrue(testCase,G.extendedNorthPole);
    verifyTrue(testCase,G.extendedSouthPole);
end

function testInterpolationNorthernCap(testCase)
    ph=0:10:360; th=0:5:40;
    [PH,TH]=meshgrid(ph,th);
    SF=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    G=SF.getGridView("interpolation",PaddingSamples=2);

    verifyEqual(testCase,G.phVec,-20:10:380);
    verifyEqual(testCase,G.thVec,-10:5:40);
    verifyTrue(testCase,G.extendedPhi);
    verifyTrue(testCase,G.extendedNorthPole);
    verifyFalse(testCase,G.extendedSouthPole);
end

function testInterpolationSouthernCap(testCase)
    ph=0:10:360; th=140:5:180;
    [PH,TH]=meshgrid(ph,th);
    SF=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    G=SF.getGridView("interpolation",PaddingSamples=2);

    verifyEqual(testCase,G.phVec,-20:10:380);
    verifyEqual(testCase,G.thVec,140:5:190);
    verifyTrue(testCase,G.extendedPhi);
    verifyFalse(testCase,G.extendedNorthPole);
    verifyTrue(testCase,G.extendedSouthPole);
end

function testInterpolationFullPhiBand(testCase)
    ph=0:10:360; th=40:5:100;
    [PH,TH]=meshgrid(ph,th);
    SF=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    G=SF.getGridView("interpolation",PaddingSamples=2);

    verifyEqual(testCase,G.phVec,-20:10:380);
    verifyEqual(testCase,G.thVec,40:5:100);
    verifyTrue(testCase,G.extendedPhi);
    verifyFalse(testCase,G.extendedNorthPole);
    verifyFalse(testCase,G.extendedSouthPole);
end

function testInterpolationPartialPatch(testCase)
    ph=-60:10:60; th=20:5:70;
    [PH,TH]=meshgrid(ph,th);
    SF=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    G=SF.getGridView("interpolation",PaddingSamples=2);

    verifyEqual(testCase,G.phVec,ph);
    verifyEqual(testCase,G.thVec,th);
    verifyFalse(testCase,G.extendedPhi);
    verifyFalse(testCase,G.extendedNorthPole);
    verifyFalse(testCase,G.extendedSouthPole);
end

function testFullThetaFieldTransformation(testCase)
    ph=0:10:360; th=0:5:180;
    [PH,TH]=meshgrid(ph,th);

    Eth=cosd(TH(:)).*cosd(PH(:));
    Eph=-sind(PH(:));
    SF=SphereField(PH(:),TH(:),Eph,Eth,1e9);

    G=SF.getGridView("fullSphere",PhiRange="symmetric",ThetaRange="360");

    phNew=0; thNew=-10;
    phOld=180; thOld=10;

    iPhNew=find(abs(G.phVec-phNew)<1e-12,1);
    iThNew=find(abs(G.thVec-thNew)<1e-12,1);

    Enew=sphPatternToCartesian(phNew,thNew,G.Eph(iThNew,iPhNew),G.Eth(iThNew,iPhNew));

    % Get the equivalent original sample directly.
    iOld=find(abs(SF.ph-phOld)<1e-12 & abs(SF.th-thOld)<1e-12,1);
    Eold=sphPatternToCartesian(phOld,thOld,SF.Eph(iOld),SF.Eth(iOld));

    verifyEqual(testCase,Enew,Eold,AbsTol=1e-12);
end

function Ecart = sphPatternToCartesian(ph,th,Eph,Eth,Er)
    if nargin<5, Er=0; end
    sth=sind(th); cth=cosd(th); sph=sind(ph); cph=cosd(ph);
    rHat=[sth*cph;sth*sph;cth];
    thHat=[cth*cph;cth*sph;-sth];
    phHat=[-sph;cph;0];
    Ecart=Eth*thHat+Eph*phHat+Er*rHat;
end

%% Resample tests
function testResampleOriginalGrid(testCase)
    SF=makeComplexGridField();
    G=SF.getGridView("stored");

    S2=SF.resample(G.phVec,G.thVec,Method="linear");

    verifyEqual(testCase,S2.Eth,SF.Eth,AbsTol=1e-12);
    verifyEqual(testCase,S2.Eph,SF.Eph,AbsTol=1e-12);
end

function SF = makeComplexGridField()
    ph=0:10:360; th=0:5:180;
    [PH,TH]=meshgrid(ph,th);
    Eth=sind(TH).*(1+0.3i*cosd(PH));
    Eph=0.2*sind(TH).*exp(1i*deg2rad(PH));
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);
end

function testResampleOutputGrid(testCase)
    SF=makeComplexGridField();
    phNew=0:5:360; thNew=0:2.5:180;
    S2=SF.resample(phNew,thNew);

    verifyEqual(testCase,S2.Np,numel(phNew)*numel(thNew));
    G=S2.getGridView("stored");
    verifyEqual(testCase,G.phVec,phNew);
    verifyEqual(testCase,G.thVec,thNew);
end

function testResampleSymmetricPhi(testCase)
    SF=makeComplexGridField();

    S1=SF.resample(0:5:355,0:5:180);
    S2=SF.resample(-180:5:175,0:5:180);

    G1=S1.getGridView("stored");
    G2=S2.getGridView("stored");

    % Shift G2 into the same physical phi ordering as G1.
    i2=zeros(size(G1.phVec));
    for ii=1:numel(G1.phVec)
        d=mod(G2.phVec-G1.phVec(ii)+180,360)-180;
        i2(ii)=find(abs(d)<1e-10,1);
    end

    verifyEqual(testCase,G2.Eth(:,i2,:),G1.Eth,AbsTol=1e-11);
    verifyEqual(testCase,G2.Eph(:,i2,:),G1.Eph,AbsTol=1e-11);
end

function testResamplePartialField(testCase)
    ph=-60:10:60; th=20:5:70;
    [PH,TH]=meshgrid(ph,th);
    Eth=cosd(PH).*sind(TH);
    Eph=sind(PH).*sind(TH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    S2=SF.resample(-50:5:50,25:2.5:65,Method="linear");

    verifyEqual(testCase,S2.Np,21*17);
end

function testResamplePartialFieldRejectsExtrapolation(testCase)
    ph=-60:10:60; th=20:5:70;
    [PH,TH]=meshgrid(ph,th);
    SF=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    verifyError(testCase,@() SF.resample(-70:10:70,20:5:70),'SphereField:ResampleOutsidePhiSupport');
    verifyError(testCase,@() SF.resample(-60:10:60,10:5:80),'SphereField:ResampleOutsideThetaSupport');
end

function testResampleNorthernCap(testCase)
    ph=0:10:360; th=0:5:40;
    [PH,TH]=meshgrid(ph,th);

    % Globally x-directed transverse field.
    Eth=cosd(TH).*cosd(PH);
    Eph=-sind(PH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    S2=SF.resample(0:5:360,0:1:40,Method="spline");

    verifyFalse(testCase,any(isnan(S2.Eth),'all'));
    verifyFalse(testCase,any(isnan(S2.Eph),'all'));
end

function testResampleEquivalentPhi(testCase)
    SF=makeComplexGridField();

    A=SF.resample(0,10:2:60,Method="makima");
    B=SF.resample(360,10:2:60,Method="makima");

    verifyEqual(testCase,A.Eth,B.Eth,AbsTol=1e-12);
    verifyEqual(testCase,A.Eph,B.Eph,AbsTol=1e-12);
end

function testResampleMethods(testCase)
    SF=makeComplexGridField();
    methods=["nearest","linear","cubic","spline","makima"];

    for method=methods
        S2=SF.resample(2:8:358,2:4:178,Method=method);
        verifyFalse(testCase,any(isnan(S2.Eth),'all'));
        verifyFalse(testCase,any(isnan(S2.Eph),'all'));
    end
end

function testResampleProperties(testCase)
    SF=makeComplexGridField();
    S2=SF.resample(0:5:360,0:5:180,Method="makima");

    verifyEqual(testCase,S2.freqHz,SF.freqHz);
    verifyEqual(testCase,S2.r,SF.r);
    verifyEqual(testCase,S2.Prad,SF.Prad);
    verifyEqual(testCase,S2.etaRad,SF.etaRad);
    verifyEqual(testCase,S2.provenance.operation,"resample");
    verifyEqual(testCase,S2.provenance.interpolationMethod,"makima");
end

% Test resample interpolation accuracy
function [Eph,Eth] = analyticSmoothField(ph,th)
    Eth=cosd(th).*cosd(ph)+0.25*cosd(th).*sind(ph);
    Eph=-sind(ph)+0.25*cosd(ph);
end

function testResampleLinearAccuracy(testCase)
    ph=0:10:360; th=0:10:180;
    [PH,TH]=meshgrid(ph,th);
    [Eph,Eth]=analyticSmoothField(PH,TH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    phNew=5:10:355; thNew=5:10:175;
    S2=SF.resample(phNew,thNew,Method="linear");

    [PH2,TH2]=meshgrid(phNew,thNew);
    [EphExact,EthExact]=analyticSmoothField(PH2,TH2);

    errEth=max(abs(S2.Eth-EthExact(:)));
    errEph=max(abs(S2.Eph-EphExact(:)));

    verifyLessThan(testCase,errEth,1e-2);
    verifyLessThan(testCase,errEph,1e-2);
end

function testResampleHigherOrderAccuracy(testCase)
    ph=0:15:360; th=0:15:180;
    [PH,TH]=meshgrid(ph,th);
    [Eph,Eth]=analyticSmoothField(PH,TH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    phNew=7.5:7.5:352.5; thNew=7.5:7.5:172.5;
    [PH2,TH2]=meshgrid(phNew,thNew);
    [EphExact,EthExact]=analyticSmoothField(PH2,TH2);

    methods=["linear","cubic","spline","makima"];
    err=zeros(size(methods));

    for ii=1:numel(methods)
        S2=SF.resample(phNew,thNew,Method=methods(ii));
        err(ii)=max([abs(S2.Eth-EthExact(:));abs(S2.Eph-EphExact(:))]);
    end

    verifyLessThan(testCase,err(2),err(1));   % cubic
    verifyLessThan(testCase,err(3),err(1));   % spline
    verifyLessThan(testCase,err(4),err(1));   % makima
end

function testResampleAccuracyAcrossPhiSeam(testCase)
    ph=0:10:360; th=0:10:180;
    [PH,TH]=meshgrid(ph,th);
    [Eph,Eth]=analyticSmoothField(PH,TH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    phNew=[355 357.5 360 362.5 365];
    thNew=20:10:160;
    S2=SF.resample(phNew,thNew,Method="spline");

    [PH2,TH2]=meshgrid(phNew,thNew);
    [EphExact,EthExact]=analyticSmoothField(PH2,TH2);

    err=max([abs(S2.Eth-EthExact(:));abs(S2.Eph-EphExact(:))]);
    verifyLessThan(testCase,err,1e-3);
end

function testResampleAccuracyNearNorthPole(testCase)
    ph=0:10:360; th=0:5:40;
    [PH,TH]=meshgrid(ph,th);
    [Eph,Eth]=analyticSmoothField(PH,TH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    phNew=0:5:355; thNew=0:1:10;
    S2=SF.resample(phNew,thNew,Method="spline");

    [PH2,TH2]=meshgrid(phNew,thNew);
    [EphExact,EthExact]=analyticSmoothField(PH2,TH2);

    err=max([abs(S2.Eth-EthExact(:));abs(S2.Eph-EphExact(:))]);
    verifyLessThan(testCase,err,1e-3);
end

function testResampleAccuracyNearSouthPole(testCase)
    ph=0:10:360; th=140:5:180;
    [PH,TH]=meshgrid(ph,th);
    [Eph,Eth]=analyticSmoothField(PH,TH);
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);

    phNew=0:5:355; thNew=170:1:180;
    S2=SF.resample(phNew,thNew,Method="spline");

    [PH2,TH2]=meshgrid(phNew,thNew);
    [EphExact,EthExact]=analyticSmoothField(PH2,TH2);

    err=max([abs(S2.Eth-EthExact(:));abs(S2.Eph-EphExact(:))]);
    verifyLessThan(testCase,err,1e-3);
end

function testResampleConvergence(testCase)
    phFine=2.5:5:357.5; thFine=2.5:5:177.5;
    [PHf,THf]=meshgrid(phFine,thFine);
    [EphExact,EthExact]=analyticSmoothField(PHf,THf);

    steps=[20 10];
    err=zeros(size(steps));

    for ii=1:numel(steps)
        ph=0:steps(ii):360; th=0:steps(ii):180;
        [PH,TH]=meshgrid(ph,th);
        [Eph,Eth]=analyticSmoothField(PH,TH);
        SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);
        S2=SF.resample(phFine,thFine,Method="linear");
        err(ii)=max([abs(S2.Eth-EthExact(:));abs(S2.Eph-EphExact(:))]);
    end

    verifyLessThan(testCase,err(2),err(1));
end

function testResampleLike(testCase)
    SF=makeComplexGridField();

    phTarget=-180:5:175;
    thTarget=0:2.5:180;
    [PH,TH]=meshgrid(phTarget,thTarget);
    target=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    S1=SF.resampleLike(target,Method="makima");
    S2=SF.resample(phTarget,thTarget,Method="makima");

    verifyEqual(testCase,S1.ph,S2.ph);
    verifyEqual(testCase,S1.th,S2.th);
    verifyEqual(testCase,S1.Eth,S2.Eth,AbsTol=1e-13);
    verifyEqual(testCase,S1.Eph,S2.Eph,AbsTol=1e-13);
end

function testResampleLikeRejectsUnstructuredTarget(testCase)
    SF=makeComplexGridField();
    target=SphereField([0;20;55],[10;30;70],zeros(3,1),ones(3,1),1e9);

    verifyError(testCase,@() SF.resampleLike(target),'SphereField:ResampleLikeRequiresStructuredTarget');
end

%% Plotting tests
function testPlotMagnitude(testCase)
    SF=makeComplexGridField();
    fig=figure('Visible','off');
    cleanup=onCleanup(@() close(fig));
    ax=axes(fig);

    h=SF.plot(Axes=ax,Quantity="magnitude",ComponentIndex=1);

    verifyClass(testCase,h,'matlab.graphics.chart.primitive.Surface');
    verifyEqual(testCase,size(h.ZData),[SF.Nth SF.Nph]);
end

function testPlotDataMatchesGetData(testCase)
    SF=makeComplexGridField();
    fig=figure('Visible','off');
    cleanup=onCleanup(@() close(fig));
    ax=axes(fig);

    h=SF.plot(Axes=ax,Quantity="real",ComponentIndex=1);
    D=SF.getData(Quantity="real",ComponentIndex=1);
    G=SF.getGridView("stored");

    verifyEqual(testCase,h.ZData,reshape(D.values,numel(G.thVec),numel(G.phVec)));
end

function testPlotAzElCoordinates(testCase)
    SF=makeComplexGridField();
    fig=figure('Visible','off');
    cleanup=onCleanup(@() close(fig));
    ax=axes(fig);

    h=SF.plot(Axes=ax,Coordinates="AzEl",Quantity="magnitude",ComponentIndex=1);
    D=SF.getData(Coordinates="AzEl",Quantity="magnitude",ComponentIndex=1);
    G=SF.getGridView("stored");

    verifyEqual(testCase,h.XData,reshape(D.x,numel(G.thVec),numel(G.phVec)));
    verifyEqual(testCase,h.YData,reshape(D.y,numel(G.thVec),numel(G.phVec)));
end

function testPlotFullSphereSymmetric360(testCase)
    SF=makeComplexGridField();
    fig=figure('Visible','off');
    cleanup=onCleanup(@() close(fig));
    ax=axes(fig);

    h=SF.plot(Axes=ax,GridView="fullSphere",PhiRange="symmetric",ThetaRange="360",...
        Quantity="real",ComponentIndex=1);

    G=SF.getGridView("fullSphere",PhiRange="symmetric",ThetaRange="360");
    verifyEqual(testCase,size(h.ZData),[G.Nth G.Nph]);
    verifyEqual(testCase,h.XData,repmat(G.phVec,G.Nth,1));
    verifyEqual(testCase,h.YData,repmat(G.thVec.',1,G.Nph));
end

function testPlotFullSphereDataMatchesGetData(testCase)
    SF=makeComplexGridField();
    fig=figure('Visible','off');
    cleanup=onCleanup(@() close(fig));
    ax=axes(fig);

    h=SF.plot(Axes=ax,GridView="fullSphere",PhiRange="symmetric",ThetaRange="360",...
        Quantity="real",ComponentIndex=1);

    D=SF.getData(GridView="fullSphere",PhiRange="symmetric",ThetaRange="360",...
        Quantity="real",ComponentIndex=1);

    verifyEqual(testCase,h.ZData,reshape(D.values,D.Nth,D.Nph));
end

% Symmetry
function testSymmetryProperties(testCase)
    SF=makeHalfSphereField(Symmetry=struct("YZ","electric"));

    verifyEqual(testCase,SF.symmetry.YZ,"electric");
    verifyEqual(testCase,SF.symmetry.XZ,"none");
    verifyEqual(testCase,SF.symmetry.XY,"none");
    verifyTrue(testCase,SF.hasSymmetry);
    verifyEqual(testCase,SF.numberOfSymmetries,1);
    verifyEqual(testCase,SF.symmetryFactor,2);
    verifyEqual(testCase,SF.symmetrySide.YZ,1);
end

function testTwoSymmetryFactor(testCase)
    SF=makeQuarterSphereField(Symmetry=struct("YZ","electric","XZ","magnetic"));

    verifyEqual(testCase,SF.numberOfSymmetries,2);
    verifyEqual(testCase,SF.symmetryFactor,4);
end

function testInvalidSymmetryPlane(testCase)
    f=@() makeHalfSphereField(Symmetry=struct("AB","electric"));
    verifyError(testCase,f,'SphereField:InvalidSymmetryPlane');
end

function testInvalidSymmetryType(testCase)
    f=@() makeHalfSphereField(Symmetry=struct("YZ","periodic"));
    verifyError(testCase,f,'SphereField:InvalidSymmetryType');
end

function testInvalidSymmetryDomain(testCase)
    f=@() makeFullSphereField(Symmetry=struct("YZ","electric"));
    verifyError(testCase,f,'SphereField:InvalidSymmetryDomain');
end

function SF=makeHalfSphereField(options)
    arguments
        options.Symmetry struct = struct()
    end

    ph=-90:10:90;
    th=0:10:180;
    [PH,TH]=meshgrid(ph,th);

    % Cartesian pattern E=[1;0;0]
    Eth=cosd(TH).*cosd(PH);
    Eph=-sind(PH);

    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,Symmetry=options.Symmetry);
end

function SF=makeQuarterSphereField(options)
    arguments
        options.Symmetry struct = struct()
    end

    ph=0:10:90;
    th=0:10:180;
    [PH,TH]=meshgrid(ph,th);

    % Cartesian pattern E=[1;0;0]
    Eth=cosd(TH).*cosd(PH);
    Eph=-sind(PH);

    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,Symmetry=options.Symmetry);
end

function SF=makeFullSphereField(options)
    arguments
        options.Symmetry struct = struct()
    end

    ph=0:10:350;
    th=0:10:180;
    [PH,TH]=meshgrid(ph,th);

    % Transverse projection of Cartesian E=[1;0;0]
    Eth=cosd(TH).*cosd(PH);
    Eph=-sind(PH);

    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,Symmetry=options.Symmetry);
end

function SF=makeSymmetrySectorField(ph,th,symmetry)
    [PH,TH]=meshgrid(ph,th);
    Eph=zeros(numel(PH),1);
    Eth=zeros(numel(PH),1);
    SF=SphereField(PH(:),TH(:),Eph,Eth,1e9,Symmetry=symmetry);
end

function testYZHalfSphereCoverage(testCase)
    SF=makeSymmetrySectorField(-90:10:90,0:10:180,struct("YZ","electric"));

    verifyFalse(testCase,SF.isFullSphere);
    verifyTrue(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyTrue(testCase,SF.representsFullSphere);
    verifyEqual(testCase,SF.symmetrySide.YZ,1);
    verifyEqual(testCase,SF.symmetryFactor,2);
end

function testYZNegativeHalfSphereCoverage(testCase)
    SF=makeSymmetrySectorField(90:10:270,0:10:180,struct("YZ","magnetic"));

    verifyTrue(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyTrue(testCase,SF.representsFullSphere);
    verifyEqual(testCase,SF.symmetrySide.YZ,-1);
end

function testXZHalfSphereCoverage(testCase)
    SF=makeSymmetrySectorField(0:10:180,0:10:180,struct("XZ","electric"));

    verifyTrue(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyEqual(testCase,SF.symmetrySide.XZ,1);
    verifyEqual(testCase,SF.symmetryFactor,2);
end

function testXYHalfSphereCoverage(testCase)
    SF=makeSymmetrySectorField(0:10:350,0:10:90,struct("XY","magnetic"));

    verifyTrue(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyTrue(testCase,SF.representsFullSphere);
    verifyEqual(testCase,SF.symmetrySide.XY,1);
    verifyEqual(testCase,SF.symmetryFactor,2);
end

function testQuarterSphereCoverage(testCase)
    S=struct("YZ","electric","XZ","magnetic");
    SF=makeSymmetrySectorField(0:10:90,0:10:180,S);

    verifyTrue(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyTrue(testCase,SF.representsFullSphere);
    verifyEqual(testCase,SF.symmetrySide.YZ,1);
    verifyEqual(testCase,SF.symmetrySide.XZ,1);
    verifyEqual(testCase,SF.numberOfSymmetries,2);
    verifyEqual(testCase,SF.symmetryFactor,4);
end

function testOctantCoverage(testCase)
    S=struct("YZ","electric","XZ","magnetic","XY","electric");
    SF=makeSymmetrySectorField(0:10:90,0:10:90,S);

    verifyTrue(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyTrue(testCase,SF.representsFullSphere);
    verifyEqual(testCase,SF.numberOfSymmetries,3);
    verifyEqual(testCase,SF.symmetryFactor,8);
end

function testPartialHalfSphereIsNotFullCoverage(testCase)
    SF=makeSymmetrySectorField(-60:10:60,20:10:160,struct("YZ","electric"));

    verifyFalse(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyFalse(testCase,SF.representsFullSphere);
end

function testIncompleteThetaIsNotFullCoverage(testCase)
    SF=makeSymmetrySectorField(-90:10:90,0:10:150,struct("YZ","electric"));

    verifyFalse(testCase,SF.hasFullSphereSymmetryCoverage);
end

function testIncompletePhiIsNotFullCoverage(testCase)
    SF=makeSymmetrySectorField(-90:10:70,0:10:180,struct("YZ","electric"));

    verifyFalse(testCase,SF.hasFullSphereSymmetryCoverage);
end

function testExplicitFullSphereRepresentsFullSphere(testCase)
    SF=makeSymmetrySectorField(0:10:350,0:10:180,struct());

    verifyTrue(testCase,SF.isFullSphere);
    verifyFalse(testCase,SF.hasSymmetry);
    verifyFalse(testCase,SF.hasFullSphereSymmetryCoverage);
    verifyTrue(testCase,SF.representsFullSphere);
    verifyEqual(testCase,SF.symmetryFactor,1);
end

% Symmetry power integration
function testIntegratePowerFullSphere(testCase)
    ph=0:2:358;
    th=0:2:180;
    [PH,TH]=meshgrid(ph,th);

    Eth=sind(TH);
    Eph=zeros(size(Eth));

    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9);
    P=SF.integratePower();

    Pexact=4*pi/(3*SF.eta0);
    verifyEqual(testCase,P,Pexact,'RelTol',5e-4);
end

function testIntegratePowerYZSymmetry(testCase)
    ph=-90:2:90;
    th=0:2:180;
    [PH,TH]=meshgrid(ph,th);

    Eth=sind(TH);
    Eph=zeros(size(Eth));

    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,...
        Symmetry=struct("YZ","magnetic"));

    P=SF.integratePower();
    Pexact=4*pi/(3*SF.eta0);

    verifyEqual(testCase,SF.symmetryFactor,2);
    verifyEqual(testCase,P,Pexact,'RelTol',5e-4);
end

function testIntegratePowerQuarterSphere(testCase)
    ph=0:2:90;
    th=0:2:180;
    [PH,TH]=meshgrid(ph,th);

    Eth=sind(TH);
    Eph=zeros(size(Eth));

    S=struct("YZ","magnetic","XZ","magnetic");
    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,Symmetry=S);

    P=SF.integratePower();
    Pexact=4*pi/(3*SF.eta0);

    verifyEqual(testCase,SF.symmetryFactor,4);
    verifyEqual(testCase,P,Pexact,'RelTol',5e-4);
end

function testIntegratePowerXYSymmetry(testCase)
    ph=0:2:358;
    th=0:2:90;
    [PH,TH]=meshgrid(ph,th);

    Eth=sind(TH);
    Eph=zeros(size(Eth));

    SF=SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,...
        Symmetry=struct("XY","electric"));

    P=SF.integratePower();
    Pexact=4*pi/(3*SF.eta0);

    verifyEqual(testCase,SF.symmetryFactor,2);
    verifyEqual(testCase,P,Pexact,'RelTol',5e-4);
end

function testSymmetryIntegrationMatchesFullSphere(testCase)
    ph=0:2:358;
    th=0:2:180;
    [PH,TH]=meshgrid(ph,th);
    SFfull=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9);

    ph=0:2:90;
    th=0:2:180;
    [PH,TH]=meshgrid(ph,th);
    S=struct("YZ","magnetic","XZ","magnetic");
    SFquarter=SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);

    verifyEqual(testCase,SFquarter.integratePower(),SFfull.integratePower(),'RelTol',1e-12);
end

function testExpandSymmetryHalfSphere(testCase)
    ph = -90:10:90;
    th = 0:10:180;
    [PH,TH] = meshgrid(ph,th);

    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,...
        Symmetry=struct("YZ","magnetic"));

    SF2 = SF.expandSymmetry();

    verifyFalse(testCase,SF2.hasSymmetry);
    verifyTrue(testCase,SF2.isFullSphere);
    verifyTrue(testCase,SF2.representsFullSphere);
    verifyEqual(testCase,SF2.Nph,36);
    verifyEqual(testCase,SF2.Nth,19);
end

function testExpandSymmetryQuarterSphere(testCase)
    ph = 0:10:90;
    th = 0:10:180;
    [PH,TH] = meshgrid(ph,th);

    S = struct("YZ","magnetic","XZ","magnetic");
    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);

    SF2 = SF.expandSymmetry();

    verifyFalse(testCase,SF2.hasSymmetry);
    verifyTrue(testCase,SF2.isFullSphere);
    verifyEqual(testCase,SF2.Nph,36);
    verifyEqual(testCase,SF2.Nth,19);
end

function testExpandSymmetryOctant(testCase)
    ph = 0:10:90;
    th = 0:10:90;
    [PH,TH] = meshgrid(ph,th);

    S = struct("YZ","magnetic","XZ","magnetic","XY","electric");
    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);

    SF2 = SF.expandSymmetry();

    verifyFalse(testCase,SF2.hasSymmetry);
    verifyTrue(testCase,SF2.isFullSphere);
    verifyEqual(testCase,SF2.Nph,36);
    verifyEqual(testCase,SF2.Nth,19);
end

function testExpandSymmetryPartialPatch(testCase)
    ph = 0:10:40;
    th = 30:10:70;
    [PH,TH] = meshgrid(ph,th);

    Eth = ones(size(PH));
    Eph = zeros(size(PH));

    SF = SphereField(PH(:),TH(:),Eph(:),Eth(:),1e9,...
        Symmetry=struct("YZ","magnetic"));

    SF2 = SF.expandSymmetry();

    verifyFalse(testCase,SF2.hasSymmetry);
    verifyFalse(testCase,SF2.isFullSphere);
    verifyFalse(testCase,SF2.representsFullSphere);
    verifyGreaterThan(testCase,SF2.Np,SF.Np);
end

function testExpandSymmetryNoSymmetry(testCase)
    SF = makeFullSphereField();
    SF2 = SF.expandSymmetry();

    verifyEqual(testCase,SF2.ph,SF.ph);
    verifyEqual(testCase,SF2.th,SF.th);
    verifyEqual(testCase,SF2.Eph,SF.Eph);
    verifyEqual(testCase,SF2.Eth,SF.Eth);
end

function testExpandSymmetryPreservesPower(testCase)
    ph = 0:5:90;
    th = 0:5:90;
    [PH,TH] = meshgrid(ph,th);

    S = struct("YZ","magnetic","XZ","magnetic","XY","electric");
    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);

    P1 = SF.integratePower();
    SF2 = SF.expandSymmetry();
    P2 = SF2.integratePower();

    verifyEqual(testCase,P2,P1,'RelTol',1e-12);
end

function testPlotExpandedSymmetry(testCase)
    ph = -90:10:90;
    th = 0:10:180;
    [PH,TH] = meshgrid(ph,th);

    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,...
        Symmetry=struct("YZ","magnetic"));

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig);

    h = SF.plot(Axes=ax,Symmetry="expand",Quantity="magnitude",ComponentIndex=1);

    verifyEqual(testCase,size(h.ZData),[19 36]);
end

function testPlotExpandedSymmetryMatchesExplicitExpansion(testCase)
    ph = 0:10:90;
    th = 0:10:180;
    [PH,TH] = meshgrid(ph,th);

    S = struct("YZ","magnetic","XZ","magnetic");
    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);
    SF2 = SF.expandSymmetry();

    fig1 = figure('Visible','off');
    fig2 = figure('Visible','off');
    cleanup = onCleanup(@() close([fig1 fig2]));

    h1 = SF.plot(Axes=axes(fig1),Symmetry="expand",Quantity="real",ComponentIndex=1);
    h2 = SF2.plot(Axes=axes(fig2),Quantity="real",ComponentIndex=1);

    verifyEqual(testCase,h1.XData,h2.XData);
    verifyEqual(testCase,h1.YData,h2.YData);
    verifyEqual(testCase,h1.ZData,h2.ZData,'AbsTol',1e-13);
end

function testPlotExpandedSymmetryFullSphereView(testCase)
    ph = 0:10:90;
    th = 0:10:90;
    [PH,TH] = meshgrid(ph,th);

    S = struct("YZ","magnetic","XZ","magnetic","XY","electric");
    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));
    ax = axes(fig);

    h = SF.plot(Axes=ax,Symmetry="expand",GridView="fullSphere",...
        PhiRange="symmetric",ThetaRange="360",Quantity="magnitude",ComponentIndex=1);

    SF2 = SF.expandSymmetry();
    D = SF2.getData(GridView="fullSphere",PhiRange="symmetric",ThetaRange="360",...
        Quantity="magnitude",ComponentIndex=1);

    verifyEqual(testCase,h.XData,reshape(D.x,D.Nth,D.Nph));
    verifyEqual(testCase,h.YData,reshape(D.y,D.Nth,D.Nph));
    verifyEqual(testCase,h.ZData,reshape(D.values,D.Nth,D.Nph),'AbsTol',1e-13);
end

function testPlotSymmetryDoesNotModifyObject(testCase)
    ph = -90:10:90;
    th = 0:10:180;
    [PH,TH] = meshgrid(ph,th);

    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,...
        Symmetry=struct("YZ","magnetic"));

    Np = SF.Np;
    symmetry = SF.symmetry;

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    SF.plot(Axes=axes(fig),Symmetry="expand",Quantity="magnitude",ComponentIndex=1);

    verifyEqual(testCase,SF.Np,Np);
    verifyEqual(testCase,SF.symmetry,symmetry);
    verifyTrue(testCase,SF.hasSymmetry);
end

function testPlotCutAngle(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=0,Quantity="real",ComponentIndex=1);

    G = SF.getGridView("stored");
    expected = [-fliplr(G.thVec(2:end)) G.thVec].';

    verifyEqual(testCase,C.angle,expected);
end

function testPlotCutMatchesData(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=0,Quantity="real",ComponentIndex=1);

    G = SF.getGridView("stored");
    SFcut = SF.resample([0 180],G.thVec);
    D = SFcut.getData(Quantity="real",ComponentIndex=1);

    Z = reshape(D.values,D.Nth,D.Nph);
    expected = [flipud(Z(2:end,2)); Z(:,1)];

    verifyEqual(testCase,C.values,expected,'AbsTol',1e-13);
end

function testPlotCutInterpolatedPhi(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=45,Quantity="magnitude",ComponentIndex=1,Method="spline");

    verifyEqual(testCase,C.phi,45);
    verifyFalse(testCase,any(isnan(C.values)));
end

function testPlotCutExpandedSymmetry(testCase)
    ph = 0:10:90;
    th = 0:10:180;
    [PH,TH] = meshgrid(ph,th);

    S = struct("YZ","magnetic","XZ","magnetic");
    SF = SphereField(PH(:),TH(:),zeros(numel(PH),1),sind(TH(:)),1e9,Symmetry=S);

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=45,Symmetry="expand",Quantity="magnitude");

    verifyEqual(testCase,C.angle(1),-180);
    verifyEqual(testCase,C.angle(end),180);
    verifyFalse(testCase,any(isnan(C.values)));
end

function testGetDataNormalizationMax(testCase)
    SF = makeComplexGridField();

    D = SF.getData(Quantity="magnitude",ComponentIndex=0,Normalize="max");

    verifyEqual(testCase,max(D.values,[],'all'),1,'AbsTol',1e-13);
    verifyEqual(testCase,D.normalization,"max");
end

function testGetDataNormalizationSpecified(testCase)
    SF = makeComplexGridField();

    D0 = SF.getData(Quantity="magnitude",ComponentIndex=1);
    D = SF.getData(Quantity="magnitude",ComponentIndex=1,...
        Normalize="specified",ReferenceValue=2);

    verifyEqual(testCase,D.values,D0.values/2,'AbsTol',1e-13);
    verifyEqual(testCase,D.referenceValue,2);
end

function testGetDataComplexNormalization(testCase)
    SF = makeComplexGridField();

    ref = 2*exp(1i*pi/3);

    D0 = SF.getData(Quantity="complex",ComponentIndex=1);
    D = SF.getData(Quantity="complex",ComponentIndex=1,...
        Normalize="specified",ReferenceValue=ref);

    verifyEqual(testCase,D.values,D0.values/ref,'AbsTol',1e-13);
    verifyEqual(testCase,D.referenceValue,ref);
end

function testGetDataMagnitudeNormalizationUsesReferenceMagnitude(testCase)
    SF = makeComplexGridField();

    ref = 2*exp(1i*pi/3);

    D0 = SF.getData(Quantity="magnitude",ComponentIndex=1);
    D = SF.getData(Quantity="magnitude",ComponentIndex=1,...
        Normalize="specified",ReferenceValue=ref);

    verifyEqual(testCase,D.values,D0.values/abs(ref),'AbsTol',1e-13);
end

function testGetDataCommonComponentNormalization(testCase)
    SF = makeComplexGridField();

    D = SF.getData(Quantity="magnitude",ComponentIndex=0,Normalize="max");

    verifyEqual(testCase,max(D.values,[],'all'),1,'AbsTol',1e-13);

    m1 = max(D.values(:,:,1),[],'all');
    m2 = max(D.values(:,:,2),[],'all');

    verifyTrue(testCase,m1 == 1 || m2 == 1);
end

function testPlotNormalization(testCase)
    SF = makeComplexGridField();
    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    h = SF.plot(Axes=axes(fig),Quantity="magnitude",ComponentIndex=1,...
        Normalize="max");

    verifyEqual(testCase,max(h.ZData,[],'all'),1,'AbsTol',1e-13);
end

function testPlotCutNormalization(testCase)
    SF = makeComplexGridField();
    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=0,Quantity="magnitude",...
        ComponentIndex=1,Normalize="max");

    verifyEqual(testCase,max(C.values),1,'AbsTol',1e-13);
end

function testNormalizationBeforeComponentSelection(testCase)
    SF = makeComplexGridField();

    Dall = SF.getData(Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=0,Normalize="max");

    D2 = SF.getData(Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=2,Normalize="max");

    verifyEqual(testCase,D2.values,Dall.values(:,:,2),'AbsTol',1e-13);
    verifyEqual(testCase,D2.referenceValue,Dall.referenceValue,'AbsTol',1e-13);
end

function testPlotCutAllComponents(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=0,Basis="Ludwig3",...
        Quantity="magnitude",ComponentIndex=0);

    verifyEqual(testCase,numel(C.componentNames),2);
    verifyEqual(testCase,C.componentNames,["Eh","Ev"]);
    verifyEqual(testCase,size(C.values,3),2);
end

function testPlotCutUsesFullFieldNormalization(testCase)
    SF = makeComplexGridField();

    D = SF.getData(Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=0,Normalize="max");

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    [~,C] = SF.plotCut(Axes=axes(fig),Phi=0,Basis="Ludwig3",...
        Quantity="magnitude",ComponentIndex=0,Normalize="max");

    verifyEqual(testCase,C.referenceValue,D.referenceValue,'AbsTol',1e-13);
end

function testPlotCutSelectedComponentUsesCommonReference(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    ax1 = axes(fig);
    [~,Call] = SF.plotCut(Axes=ax1,Phi=0,Basis="Ludwig3",...
        Quantity="magnitude",ComponentIndex=0,Normalize="max");

    clf(fig);
    ax2 = axes(fig);
    [~,C2] = SF.plotCut(Axes=ax2,Phi=0,Basis="Ludwig3",...
        Quantity="magnitude",ComponentIndex=2,Normalize="max");

    verifyEqual(testCase,C2.values,Call.values(:,:,2),'AbsTol',1e-13);
    verifyEqual(testCase,C2.referenceValue,Call.referenceValue,'AbsTol',1e-13);
end

function testPlotPrincipalCutsScalar(testCase)
    SF = makeFullSphereField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    h = SF.plotPrincipalCuts(Axes=axes(fig),...
        Quantity="directivity",Scale="dB10",Normalize="none");

    verifySize(testCase,h,[3 1]);
    verifyTrue(testCase,all(isgraphics(h)));
end

function testPlotPrincipalCutsLudwig3Components(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    h = SF.plotPrincipalCuts(Axes=axes(fig),...
        Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=0,Scale="linear",Normalize="none");

    verifySize(testCase,h,[3 2]);
    verifyTrue(testCase,all(isgraphics(h),"all"));

    % Same cut -> same colour, different component -> different line style.
    for ii = 1:3
        verifyEqual(testCase,h(ii,1).Color,h(ii,2).Color,'AbsTol',1e-13);
        verifyNotEqual(testCase,string(h(ii,1).LineStyle),string(h(ii,2).LineStyle));
    end
end

function testPlotPrincipalCutsUsesFullFieldNormalization(testCase)
    SF = makeComplexGridField();

    D = SF.getData(Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=0,Normalize="max",Scale="linear");

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    h = SF.plotPrincipalCuts(Axes=axes(fig),...
        Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=0,Normalize="max",Scale="linear");

    % Reconstruct each plotted cut independently using the already-known
    % full-field reference and compare directly with plotted YData.
    phi = [0 45 90];

    for ii = 1:numel(phi)
        ax2 = axes('Parent',fig,'Visible','off');

        [~,C] = SF.plotCut(Axes=ax2,Phi=phi(ii),...
            Basis="Ludwig3",Quantity="magnitude",...
            ComponentIndex=0,Normalize="specified",...
            ReferenceValue=D.referenceValue,Scale="linear");

        Y = reshape(C.values,numel(C.angle),[]);

        for jj = 1:size(Y,2)
            verifyEqual(testCase,h(ii,jj).YData(:),Y(:,jj),'AbsTol',1e-13);
        end

        delete(ax2);
    end
end

function testPlotPrincipalCutsExpandedSymmetry(testCase)
    S=struct("YZ","electric","XZ","magnetic","XY","electric");
    SF=makeSymmetrySectorField(0:10:90,0:10:90,S);

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    h = SF.plotPrincipalCuts(Axes=axes(fig),...
        Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=0,Normalize="none",Scale="linear",...
        Symmetry="expand");

    verifySize(testCase,h,[3 2]);
    verifyTrue(testCase,all(isgraphics(h),"all"));

    for ii = 1:numel(h)
        verifyTrue(testCase,all(isfinite(h(ii).XData)));
        verifyTrue(testCase,all(isfinite(h(ii).YData)));
    end
end

function testPlotPrincipalCutsSelectedComponent(testCase)
    SF = makeComplexGridField();

    fig = figure('Visible','off');
    cleanup = onCleanup(@() close(fig));

    h = SF.plotPrincipalCuts(Axes=axes(fig),...
        Basis="Ludwig3",Quantity="magnitude",...
        ComponentIndex=2,Normalize="max",Scale="dB20");

    verifySize(testCase,h,[3 1]);
    verifyTrue(testCase,all(isgraphics(h)));
end

function testPowerPatternStoredU(testCase)
    ph = [0;90;180];
    th = [30;30;30];
    U = [1;2;4];

    SF = SphereField(ph,th,[],[],1e9,PowerPattern=U);

    D = SF.getData(Quantity="U");

    testCase.verifyEqual(D.values,U,'AbsTol',1e-12);
    testCase.verifyTrue(SF.isPowerPattern);
    testCase.verifyEqual(SF.powerPattern,U,'AbsTol',1e-12);
end

function testPowerPatternBlocksGetField(testCase)
    SF = SphereField([0;90],[30;30],[],[],1e9,...
        PowerPattern=[1;2]);

    testCase.verifyError(@() SF.getField(),...
        'SphereField:PowerPatternVectorOperation');
end

function testPowerPatternBlocksCartesianField(testCase)
    SF = SphereField([0;90],[30;30],[],[],1e9,...
        PowerPattern=[1;2]);

    testCase.verifyError(@() SF.getCartesianPattern(),...
        'SphereField:PowerPatternVectorOperation');
end

function testPowerPatternIntegratedIsotropic(testCase)
    ph = 0:5:360;
    th = 0:5:180;
    [PH,TH] = meshgrid(ph,th);

    U0 = 2;
    U = U0*ones(size(PH));

    SF = SphereField(PH(:),TH(:),[],[],1e9,...
        PowerPattern=U(:));

    P = SF.integratePower();

    testCase.verifyEqual(P,4*pi*U0,'RelTol',2e-3);
end

function testPowerPatternDirectivityIsotropic(testCase)
    ph = 0:5:360;
    th = 0:5:180;
    [PH,TH] = meshgrid(ph,th);

    U = 3*ones(size(PH));

    SF = SphereField(PH(:),TH(:),[],[],1e9,...
        PowerPattern=U(:));

    D = SF.getData(Quantity="directivity",...
        PowerSource="integrated");

    testCase.verifyEqual(D.values,...
        ones(size(D.values)),...
        'RelTol',2e-3);
end

function testReadGRASPPhThCircular(testCase)
    p = fileparts(mfilename("fullpath"));
    pathName = [p,'\..\..\data\SimPatterns\GRASPgrd\FF_phth_circular_pos180.grd'];

    SF = SphereField.readGRASPgrd(pathName);

    % Basic dimensions
    testCase.verifyEqual(SF.Nf,11);
    testCase.verifyEqual(SF.freqHz,...
        (1.00:0.05:1.50)*1e9,...
        'AbsTol',1);

    % Structured PhTh grid
    testCase.verifyTrue(SF.isStructured);

    G = SF.getGridView("stored");

    testCase.verifyEqual(G.phVec,0:10:360,'AbsTol',1e-12);
    testCase.verifyEqual(G.thVec,0:5:180,'AbsTol',1e-12);

    testCase.verifySize(G.Eth,[37 37 11]);
    testCase.verifySize(G.Eph,[37 37 11]);
end

function testReadGRASPPeriodicPhi(testCase)
    p = fileparts(mfilename("fullpath"));
    pathName = [p,'\..\..\data\SimPatterns\GRASPgrd\FF_phth_circular_pos180.grd'];

    SF = SphereField.readGRASPgrd(pathName);
    G = SF.getGridView("stored");

    testCase.verifyEqual(G.Eth(:,1,:),G.Eth(:,end,:),...
        'RelTol',1e-10,'AbsTol',1e-12);

    testCase.verifyEqual(G.Eph(:,1,:),G.Eph(:,end,:),...
        'RelTol',1e-10,'AbsTol',1e-12);
end

function testReadGRASPAgainstFarField(testCase)
    p = fileparts(mfilename("fullpath"));
    pathName = [p,'\..\..\data\SimPatterns\GRASPgrd\FF_phth_circular_pos180.grd'];

    SF = SphereField.readGRASPgrd(pathName);
    FF = FarField.readGRASPgrd(char(pathName));

    % Frequencies
    testCase.verifyEqual(SF.freqHz,FF.freqHz,...
        'RelTol',1e-12);

    % Compare canonical spherical fields.
    %
    % Adjust this part depending on exactly how FarField exposes
    % spherical linear components.
    FFs = FF.coor2spherical(true);
    FFs = FFs.pol2linear;

    testCase.verifyEqual(SF.th,rad2deg(FFs.y),...
        'AbsTol',1e-10);

    testCase.verifyEqual(SF.ph,rad2deg(FFs.x),...
        'AbsTol',1e-10);

    testCase.verifyEqual(SF.Eth,FFs.E1,...
        'RelTol',1e-10,'AbsTol',1e-12);

    testCase.verifyEqual(SF.Eph,FFs.E2,...
        'RelTol',1e-10,'AbsTol',1e-12);
end

function testReadGRASPPowerInvariant(testCase)
    p = fileparts(mfilename("fullpath"));
    pathName = [p,'\..\..\data\SimPatterns\GRASPgrd\FF_phth_circular_pos180.grd'];

    SF = SphereField.readGRASPgrd(pathName);
    FF = FarField.readGRASPgrd(char(pathName));

    Usf = SF.getU;

    FFs = FF.coor2spherical(true);
    FFs = FFs.pol2linear;

    Uff = (abs(FFs.E1).^2 + abs(FFs.E2).^2)/(2*FF.eta0);

    testCase.verifyEqual(Usf,Uff,...
        'RelTol',1e-9,'AbsTol',1e-10);
end

function testNativeDirCosGridDetected(testCase)

    uVec = -0.2:0.1:0.2;
    vVec = -0.15:0.05:0.15;

    [U,V] = meshgrid(uVec,vVec);
    W = sqrt(1 - U.^2 - V.^2);

    ph = atan2d(V,U);
    th = atan2d(hypot(U,V),W);

    ph = ph(:);
    th = th(:);

    Eth = ones(size(ph));
    Eph = zeros(size(ph));

    SF = SphereField(ph,th,Eph,Eth,1e9,...
        NativeGrid="DirCos");

    testCase.verifyTrue(SF.hasNativeStructuredGrid);
end

function testNativeAzElGridDetected(testCase)

    azVec = -40:10:40;
    elVec = -30:10:30;

    [AZ,EL] = meshgrid(azVec,elVec);

    u = sind(AZ).*cosd(EL);
    v = sind(EL);
    w = cosd(AZ).*cosd(EL);

    ph = atan2d(v,u);
    th = atan2d(hypot(u,v),w);

    ph = ph(:);
    th = th(:);

    Eth = ones(size(ph));
    Eph = zeros(size(ph));

    SF = SphereField(ph,th,Eph,Eth,1e9,...
        NativeGrid="AzEl");

    testCase.verifyTrue(SF.hasNativeStructuredGrid);
end

function testNativeAzElPoleWarning(testCase)

    azVec = -180:30:180;
    elVec = -90:30:90;

    [AZ,EL] = meshgrid(azVec,elVec);

    u = sind(AZ).*cosd(EL);
    v = sind(EL);
    w = cosd(AZ).*cosd(EL);

    ph = atan2d(v,u);
    th = atan2d(hypot(u,v),w);

    ph = ph(:);
    th = th(:);

    Eth = ones(size(ph));
    Eph = zeros(size(ph));

    lastwarn("");
    SF = SphereField(ph,th,Eph,Eth,1e9,...
        NativeGrid="AzEl");

    [~,warnID] = lastwarn;

    testCase.verifyEqual(warnID,...
        'SphereField:NativeGridPoleAmbiguity');

    testCase.verifyFalse(SF.hasNativeStructuredGrid);
end

function testNativeGridUnstructuredFallback(testCase)

    ph = [0;17;43;91;137];
    th = [10;24;39;61;83];

    Eth = ones(size(ph));
    Eph = zeros(size(ph));

    lastwarn("");
    SF = SphereField(ph,th,Eph,Eth,1e9,...
        NativeGrid="DirCos");

    [~,warnID] = lastwarn;

    testCase.verifyEqual(warnID,...
        'SphereField:NativeGridNotStructured');

    testCase.verifyFalse(SF.hasNativeStructuredGrid);
end

function testNoNativeGridByDefault(testCase)

    phVec = 0:30:330;
    thVec = 0:15:180;

    [PH,TH] = meshgrid(phVec,thVec);

    Eth = sind(TH(:));
    Eph = zeros(size(Eth));

    SF = SphereField(PH(:),TH(:),Eph,Eth,1e9);

    testCase.verifyTrue(SF.isStructured);
    testCase.verifyFalse(SF.hasNativeStructuredGrid);
end

function testNativeAzElGridWithPolesAndAxes(testCase)

    azVec = -180:30:180;
    elVec = -90:30:90;

    [AZ,EL] = meshgrid(azVec,elVec);

    u = sind(AZ).*cosd(EL);
    v = sind(EL);
    w = cosd(AZ).*cosd(EL);

    ph = atan2d(v,u);
    th = atan2d(hypot(u,v),w);

    Eth = ones(numel(ph),1);
    Eph = zeros(numel(ph),1);

    SF = SphereField(ph(:),th(:),Eph,Eth,1e9,...
        NativeGrid="AzEl",...
        NativeXVec=azVec,...
        NativeYVec=elVec);

    testCase.verifyTrue(SF.hasNativeStructuredGrid);
    testCase.verifyEqual(SF.nativeGridInfo.coordinates,"AzEl");
    testCase.verifyEqual(SF.nativeGridInfo.xVec,azVec);
    testCase.verifyEqual(SF.nativeGridInfo.yVec,elVec);
end

function testNativeGridAxesMismatch(testCase)

    azVec = -40:10:40;
    elVec = -30:10:30;

    [AZ,EL] = meshgrid(azVec,elVec);

    u = sind(AZ).*cosd(EL);
    v = sind(EL);
    w = cosd(AZ).*cosd(EL);

    ph = atan2d(v,u);
    th = atan2d(hypot(u,v),w);

    Eth = ones(numel(ph),1);
    Eph = zeros(numel(ph),1);

    badAz = azVec;
    badAz(3) = badAz(3) + 1;

    f = @() SphereField(ph(:),th(:),Eph,Eth,1e9,...
        NativeGrid="AzEl",...
        NativeXVec=badAz,...
        NativeYVec=elVec);

    testCase.verifyError(f,...
        'SphereField:NativeGridDirectionMismatch');
end