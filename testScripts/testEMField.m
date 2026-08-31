function tests = testEMField
%TESTEMFIELD Unit tests for EMField.

    tests = functiontests(localfunctions);

end


%% Constructor

function testFrequencyStorage(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    verifyEqual(testCase, ...
        f.freqHz, ...
        [1e9 2e9 3e9]);

    verifyEqual(testCase, f.Nf, 3);

end


function testColumnFrequencyConvertedToRow(testCase)

    f = TestEMFieldImpl([1e9; 2e9; 3e9]);

    verifyEqual(testCase, ...
        f.freqHz, ...
        [1e9 2e9 3e9]);

    verifyTrue(testCase, isrow(f.freqHz));

end


function testMetadata(testCase)

    metadata.description = "Test field";
    metadata.ID = 42;

    f = TestEMFieldImpl(1e9, ...
        Metadata = metadata);

    verifyEqual(testCase, ...
        f.metadata, metadata);

end


function testProvenance(testCase)

    provenance.solver = "FEKO";
    provenance.filename = "test.ffe";

    f = TestEMFieldImpl(1e9, ...
        Provenance = provenance);

    verifyEqual(testCase, ...
        f.provenance, provenance);

end


%% Physical quantities

function testWavelength(testCase)

    freq = [1e9 2e9 10e9];

    f = TestEMFieldImpl(freq);

    expected = EMField.c0 ./ freq;

    verifyEqual(testCase, ...
        f.wavelength, expected, ...
        AbsTol = 1e-15);

end


function testWavenumber(testCase)

    freq = [1e9 2e9 10e9];

    f = TestEMFieldImpl(freq);

    expected = 2*pi*freq / EMField.c0;

    verifyEqual(testCase, ...
        f.wavenumber, expected, ...
        RelTol = 1e-14);

end


%% Frequency indexing

function testFrequencyIndex(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    verifyEqual(testCase, ...
        f.frequencyIndex(2e9), 2);

end


function testFrequencyIndexVector(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    idx = f.frequencyIndex([3e9 1e9]);

    verifyEqual(testCase, idx, [3 1]);

end


function testFrequencyTolerance(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    idx = f.frequencyIndex( ...
        2e9 + 1, ...
        AbsoluteTolerance = 2);

    verifyEqual(testCase, idx, 2);

end


function testFrequencyNotFound(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    verifyError(testCase, ...
        @() f.frequencyIndex(2.5e9), ...
        'EMField:FrequencyNotFound');

end


function testNearestFrequency(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    verifyEqual(testCase, ...
        f.nearestFrequencyIndex(2.4e9), 2);

    verifyEqual(testCase, ...
        f.nearestFrequencyIndex(2.6e9), 3);

end


function testGetFrequency(testCase)

    f = TestEMFieldImpl([1e9 2e9 3e9]);

    verifyEqual(testCase, ...
        f.getFrequency([1 3]), ...
        [1e9 3e9]);

end


%% Input validation

function testNegativeFrequencyRejected(testCase)

    verifyError(testCase, ...
        @() TestEMFieldImpl(-1e9), ...
        ?MException);   % We can tighten this later

end


function testZeroFrequencyRejected(testCase)

    verifyError(testCase, ...
        @() TestEMFieldImpl(0), ...
        ?MException);

end


function testInvalidFrequencyIndex(testCase)

    f = TestEMFieldImpl([1e9 2e9]);

    verifyError(testCase, ...
        @() f.getFrequency(3), ...
        'EMField:FrequencyIndexOutOfRange');

end

