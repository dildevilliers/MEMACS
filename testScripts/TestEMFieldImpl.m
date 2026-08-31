%% Concrete implementation
classdef TestEMFieldImpl < EMField

    methods

        function obj = TestEMFieldImpl(freqHz, options)

            arguments
                freqHz double

                options.Metadata struct = struct()
                options.Provenance struct = struct()
            end

            obj@EMField(freqHz, ...
                Metadata = options.Metadata, ...
                Provenance = options.Provenance);

        end

    end

end
