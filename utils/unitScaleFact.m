function scaleFact = unitScaleFact(unitText)
% Very basic unitconversion - will not always work!

if any(strcmp(unitText,{'Hz','F','H'})) || strcmp(unitText,'m')
    scaleFact = 1;
else
    switch unitText(1)
        case 'p'
            scaleFact = 1e-12;
        case 'n'
            scaleFact = 1e-9;
        case 'u'
            scaleFact = 1e-6;
        case 'm'
            scaleFact = 1e-3;
        case 'k'
            scaleFact = 1e3;
        case 'M'
            scaleFact = 1e6;
        case 'G'
            scaleFact = 1e9;
        case 'T'
            scaleFact = 1e12;
    end
end