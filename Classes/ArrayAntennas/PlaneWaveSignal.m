classdef PlaneWaveSignal
    %PLANEWAVESIGNAL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ph(1,1) double {mustBeReal, mustBeFinite} = 0       % azimuth angle in rad
        th(1,1) double {mustBeReal, mustBeFinite} = 0       % polar angle in rad
        freq(1,1) double {mustBeReal, mustBeFinite} = 1     % in Hz
        sigPower(1,1) double {mustBeReal, mustBeFinite} = 0 % in dBm (Isotropic source at origin would recieve this power from this signal)
        sigPhase(1,1) double {mustBeReal, mustBeFinite} = 0 % in rad
        sigType(1,:) char {mustBeMember(sigType,{'compExp','noise','GPS_L1', 'narrow_noise'})} = 'compExp'
        b(1,1) double {mustBeReal, mustBeFinite} = 0       % bandwidth of signal
        n(1,1) double {mustBeReal, mustBeFinite} = 1       % number of frequency elements
    end
    
    properties (SetAccess = private)
       k    % Propagation constant of the signal [3x1]
       lambda % Wavelength of the signal in M
       P     % Signal power in Watt
    end
    
    methods
        function obj = PlaneWaveSignal(sigType,freq,th,ph,sigPower,sigPhase,b,n)
            if nargin >= 1
                obj.sigType = sigType;
            end
            if nargin >= 2
                obj.freq = freq;
            end
            if nargin >= 3
                obj.th = th;
            end
            if nargin >= 4
                obj.ph = ph;
            end
            if nargin >= 5
                obj.sigPower = sigPower;
            end
            if nargin >= 6
                obj.sigPhase = sigPhase;
            end
            if nargin >= 7
                obj.b = b;
                if obj.b == 0
                    n = 1;
                end
            end
            if nargin >= 8
                obj.n = n;
            end
            fr = linspace(freq-obj.b,freq+obj.b,obj.n);
            obj.lambda = physconst('lightspeed')./fr;
            [u,v,w] = PhTh2DirCos(obj.ph,obj.th);
            obj.k = 2*pi./obj.lambda.*[u;v;w];
%             obj.k = obj.k(:);
%             obj.k = 2*pi./obj.lambda.*[cos(obj.ph).*sin(obj.th);sin(obj.ph).*sin(obj.th);cos(obj.th)];
            obj.P = lin10(obj.sigPower-30);      % Signal power in W
        end
        
        function St = getSignal(obj,t)
            % Test size of t
            assert(min(size(t)) == 1,'Error, t must be a vector, not a matrix');
            assert(isreal(t),'Error, t must be real');
            
            t = t(:);
            switch obj.sigType
                case 'compExp'
                    St = sqrt(obj.P).*exp(1i.*(2*pi*obj.freq.*t + obj.sigPhase));
                case 'noise'
                    St = randn(size(t))+1j*randn(size(t));
                    StP = rms(St);
                    St = St.*sqrt(obj.P)./StP;
                case 'narrow_noise'              
                    if obj.b == 0
                        range = 0;
                        L = 1;
                    else
                        L = obj.n;
                        range = linspace(-obj.b,obj.b,L);
                    end
                    phase_deg = -180 + 2*180.*randn(L,1);
                    phase_rad = deg2rad(phase_deg);
                    w = 2*pi*range;
                    omega_c = 2*pi*obj.freq;
                    St = sqrt(obj.P/L)*exp(1j*((omega_c + w).'*t.' + phase_rad));                   
                case 'GPS_L1'
                     error(['sigType: ', obj.sigType, ' not implemented yet'])
            end
        end
    end
end

