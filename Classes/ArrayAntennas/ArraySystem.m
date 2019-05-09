classdef ArraySystem

    properties
        % Timing
        freqSamp(1,1) double {mustBeReal, mustBeFinite} = 16.368e6   % Sample rate in Hz
        Nt(1,1) double {mustBeReal, mustBeFinite} = 1600;            % Number of time samples
    end
    
    properties (SetAccess = private)
        % Elements
        antPos(1,:) pnt3D = pnt3D([0,1],0,0)   % Antenna positions in 3D points - internal vector
        channelPhasors(1,:) double {mustBeFinite} = 1 % Vector of complex channel phasors for scaling/calibration
        couplingMatrix(:,:) double {mustBeFinite} = 1 % Mutual coupling matrix between channels
        % Receiver
        noisePower(1,1) double {mustBeReal, mustBeFinite} = -110 % at the input and in dBm
        LNAGain(1,1) double {mustBeReal, mustBeFinite} = 20    % LNA gain in dB
        IFGain(1,1) double {mustBeReal, mustBeFinite} = 130    % IF gain in dB
        freqLO(1,1) double {mustBeReal, mustBeFinite} = 1.571324e9        % LO frequency in Hz
        % ADC
        Nbits(1,1) double {mustBeInteger, mustBeFinite} = 1 
        maxV(1,1) double {mustBeReal, mustBeFinite} = 1
        minV(1,1) double {mustBeReal, mustBeFinite} = -1
        ADCtype(1,:) char {mustBeMember(ADCtype,{'binary','sign_mag','2scomp'})} = 'binary' 
        
        elements
        receiver
        adc
        delT
        t
        Nant
    end
    
    methods
         % Constructor
         function obj = ArraySystem(antPos,channelPhasors,couplingMatrix,noisePower,LNAGain,IFGain,freqLO,freqSamp,Nt,Nbits,maxV,minV,ADCtype)
             if nargin >= 1
                 obj = obj.setAntPos(antPos);
             end
             if nargin >= 2
                 obj = obj.setChannelPhasors(channelPhasors); 
             end
             if nargin >= 3
                 obj = obj.setCouplingMatrix(couplingMatrix);
             end
             if nargin >= 4
                 obj = obj.setNoisePower(noisePower);
             end
             if nargin >= 5
                 obj = obj.setLNAGain(LNAGain);
             end
             if nargin >= 6
                 obj = obj.setIFGain(IFGain);
             end
             if nargin >= 7
                 obj = obj.setFreqLO(freqLO);
             end
             if nargin >= 8, obj.freqSamp = freqSamp; end
             if nargin >= 9, obj.Nt = Nt; end
             if nargin >= 10
                 obj = obj.setNbits(Nbits);
             end
             if nargin >=11
                 obj = obj.setMaxV(maxV);
             end
             if nargin >=12
                 obj = obj.setMinV(minV);
             end
             if nargin >=13
                 obj = obj.setADCtype(ADCtype); end
             
             obj.elements = ArrayElements(obj.antPos,obj.channelPhasors,obj.couplingMatrix);
             obj.receiver = ArrayReceiver(obj.noisePower,obj.LNAGain,obj.IFGain,obj.freqLO,obj.couplingMatrix);
             obj.adc = ArrayADC(obj.Nbits,obj.maxV,obj.minV,obj.ADCtype);
             obj.delT = 1/obj.freqSamp;
             t0 = 0;   % Hardcode for now - no reason to change I think...
             obj.t = t0:obj.delT:(t0+obj.delT*(obj.Nt-1));
             obj.Nant = size(obj.antPos.pointMatrix,2);
         end
         
         % Parameter setters
         function obj = setAntPos(obj,antPos)
             [obj.antPos,obj.elements.antPos] = deal(antPos);
         end
         function obj = setChannelPhasors(obj,channelPhasors)
             [obj.channelPhasors,obj.elements.channelPhasors] = deal(channelPhasors);
         end
         function obj = setCouplingMatrix(obj,couplingMatrix)
             [obj.couplingMatrix,obj.elements.couplingMatrix,obj.receiver.couplingMatrix] = deal(couplingMatrix);
         end
         function obj = setNoisePower(obj,noisePower)
             [obj.noisePower,obj.receiver.noisePower] = deal(noisePower);
         end
         function obj = setLNAGain(obj,LNAGain)
             [obj.LNAGain,obj.receiver.LNAGain] = deal(LNAGain);
         end
         function obj = setIFGain(obj,IFGain)
             [obj.IFGain,obj.receiver.IFGain] = deal(IFGain);
         end
         function obj = setFreqLO(obj,freqLO)
             [obj.freqLO,obj.receiver.freqLO] = deal(freqLO);
         end
         function obj = setNbits(obj,Nbits)
             [obj.Nbits,obj.adc.Nbits] = deal(Nbits);
         end
         function obj = setMaxV(obj,maxV)
             [obj.maxV,obj.adc.maxV] = deal(maxV);
         end
         function obj = setMinV(obj,minV)
             [obj.minV,obj.adc.minV] = deal(minV);
         end
         function obj = setADCtype(obj,ADCtype)
             [obj.ADCtype,obj.adc.ADCtype] = deal(ADCtype);
         end
         
         function x = getPortSignal(obj,s,Qtype)
             % Returns the digitised signals at the antenna ports as a
             % matrix of size [Nant, Nt].  s is an array of PlaneWaveSignal
             % objects
             % Qtype describes how to handle the Q channel:
             %  0: use the actual Q
             %  [...-2,-1,1,2,...]: Integer shift by this number
             %  +-inf: Hilbert transform
             
             if nargin < 3
                 Qtype = 0;
             end
             assert(isscalar(Qtype),'Error: Qtype must be scalar');
             
             % Get signals after elements
             portSigMat = obj.elements.portSignals(s,obj.t);
             % Put them through receiver
             sn = obj.receiver.sigRec(portSigMat,obj.t);
             % Digitise
             xi = obj.adc.ADC(real(sn));
             if Qtype == 0
                 xq = obj.adc.ADC(imag(sn));
             elseif isfloat(Qtype)
                 xq = circshift(xi,Qtype,2);
             elseif isinf(Qtype)
                 xq = hilbert(xi);
             else
                 error(['Unknown Qtype'])
             end
             x = xi + 1i.*xq;
         end
         
         function plotPortSignal(obj,s,portNumber)
             if nargin < 3
                 portNumber = 1;
             end
             x = obj.getPortSignal(s);
             xP = real(x(portNumber,:));
             plot(obj.t,xP), grid on, hold on
             xlabel('t (s)')
         end
    end
    
end
