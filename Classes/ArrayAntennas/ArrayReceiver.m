classdef ArrayReceiver
    %ARRAYRECEIVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        noisePower(1,1) double {mustBeReal, mustBeFinite} = -110 % at the input and in dBm
        LNAGain(1,1) double {mustBeReal, mustBeFinite} = 20    % LNA gain in dB
        IFGain(1,1) double {mustBeReal, mustBeFinite} = 19    % IF gain in dB
        freqLO(1,1) double {mustBeReal, mustBeFinite} = 1        % LO frequency in Hz
        couplingMatrix(:,:) double {mustBeFinite} = 1 % Mutual coupling matrix between channels
    end
    
    properties (SetAccess = private)
        P     % Noise power in Watt
    end
    
    methods
        function obj = ArrayReceiver(noisePower,LNAGain,IFGain,freqLO,couplingMatrix)
            if nargin >= 1
                obj.noisePower = noisePower;
            end
            if nargin >= 2
                obj.LNAGain = LNAGain;
            end
            if nargin >= 3
                obj.IFGain = IFGain;
            end
            if nargin >= 4
                obj.freqLO = freqLO;
            end
            if nargin >= 5
                obj.couplingMatrix = couplingMatrix;
            end
            obj.P = lin10(obj.noisePower-30);      % Noise power in W
        end
        
        function Nt = getNoise(obj,Nsamp,Nchan)
            % Return complex noise matrix with [NsampxNchan] elements at
            % input of LNA
            if nargin < 3
                Nchan = 1;
            end
            if nargin < 2
                Nsamp = 1001;
            end
            Nt = randn(Nsamp,Nchan) + 1i.*randn(Nsamp,Nchan); % Complex noise
            NtP = rms(Nt);
            Nt = bsxfun(@rdivide,Nt.*sqrt(obj.P),NtP);
        end
        
        function [SN,S,N] = sigLNA(obj,St)
            % Calculates the signal and noise after the LNA - adds noise and applies
            % gain.  St is a signal matrix, with Nrows the number of
            % samples and Ncol the number of channels
            % SN = S + N;
            [Nsamp,Nchan] = size(St);
            Nt = getNoise(obj,Nsamp,Nchan);
            N = lin20(obj.LNAGain).*Nt;     % Voltage gain - lin20
            S = lin20(obj.LNAGain).*St;
            SN = S + N;
        end
        
        function SNif = sigIF(obj,SNrf,t)
            % Mixes down the signal after the LNA, and applies the IF gain
            % Test size of t
            assert(min(size(t)) == 1,'Error, t must be a vector, not a matrix');
            assert(isreal(t),'Error, t must be real');
            
            t = t(:).';
            SNif = lin20(obj.IFGain).*bsxfun(@times,SNrf,exp(-1i*2*pi*obj.freqLO.*t)); % Voltage gain - lin20
        end
        
        function [sn,si,sq] = sigRec(obj,portSigMat,t)
            % Full effect of receiver on signal - outputs the IF signal
            % plus noise from the port signal matrix
            snCoupled = obj.couplingMatrix*portSigMat;
            snLNA = obj.sigLNA(snCoupled);
            sn = obj.sigIF(snLNA,t);
            si = real(sn);
            sq = imag(sn);
        end
    end
end

