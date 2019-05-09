classdef ArrayADC
    %ARRAYA2D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nbits(1,1) double {mustBeInteger, mustBeFinite} = 1 
        maxV(1,1) double {mustBeReal, mustBeFinite} = 1
        minV(1,1) double {mustBeReal, mustBeFinite} = -1
        ADCtype(1,:) char {mustBeMember(ADCtype,{'binary','sign_mag','2scomp'})} = 'binary' % Default of unsigned binary
    end
    
    properties (SetAccess = private)
        delV
        Nthresholds
        thresholdLevels
        outLevels
        outBinary  % Binary output
    end
    
    methods
        function obj = ArrayADC(Nbits,maxV,minV,ADCtype)
            if nargin >= 1
                obj.Nbits = Nbits;
            end
            if nargin >=2 
                obj.maxV = maxV;
                obj.minV = -maxV;
            end
            if nargin >= 3
                obj.minV = minV;
            end
            if nargin >= 4
                obj.ADCtype = ADCtype;
            end
            assert(obj.maxV > obj.minV, 'Error: maxV must be larger than minV');
            obj.delV = obj.maxV - obj.minV;
            obj.Nthresholds = 2^obj.Nbits-1;
            delThNorm = 1/(obj.Nthresholds+1);
            obj.thresholdLevels = ((delThNorm:delThNorm:obj.Nthresholds*delThNorm)*obj.delV + obj.minV).';
            obj.outLevels = linspace(obj.minV,obj.maxV,obj.Nthresholds+1).';
            switch obj.ADCtype
                case 'binary'
                    obj.outBinary = dec2bin(0:(2^obj.Nbits - 1),Nbits);
                otherwise
                    error(['ADCtype: ', obj.ADCtype, ' not implemented yet']);
            end
            
        end
        
        function b = getBinarySignal(obj,s)
            % Returns the binary samples of the analog signals in s as b.
            % The size of b is [size(s)xNbits].  s has size [Nant x Nsamples]
            binNumber = zeros(size(s));
            binNumber(s < obj.thresholdLevels(1)) = 1;
            for tt = 2:obj.Nthresholds
                binNumber(s >= obj.thresholdLevels(tt-1) & s <= obj.thresholdLevels(tt)) = tt;
            end
            binNumber(s > obj.thresholdLevels(obj.Nthresholds)) = obj.Nthresholds+1;
            b = obj.outBinary(binNumber,:);
            b = reshape(b,size(s,1),size(s,2),obj.Nbits);
        end
        
        function x = getDigitalSignal(obj,b)
            % Returns the digital signal of the binary samples b in x.
            % x is of size [Nsig x Nsamp]. The size of b is [size(s) x Nbits].  
            [Nsig,Nsamp,Nbit] = size(b);
            assert(Nbit == obj.Nbits,'Error: Last dimension of b should be of size obj.Nbits');
            bV = reshape(b,Nsig*Nsamp,obj.Nbits);
            switch obj.ADCtype
                case 'binary'
                    binNumber = bin2dec(bV) + 1;
            end
            x = obj.outLevels(binNumber);
            x = reshape(x,Nsig,Nsamp);
        end
        
        function x = ADC(obj,s)
           % Full sampling of signal s into digital signal x
           b = getBinarySignal(obj,s);
           x = getDigitalSignal(obj,b);
        end
        
        function plot(obj,showBin)
            if nargin < 2
                showBin = 1;
            end
            x = obj.thresholdLevels;
            % Build the line segments for the plot
            [xp,yp] = deal(zeros(1,2^(obj.Nbits+1)));
            xp(1) = obj.minV;
            yp(1) = obj.outLevels(1);
            for ii = 2:(2^(obj.Nbits+1) - 1)
                if mod(ii,2)
                    xp(ii) = xp(ii-1);
                    yp(ii) = obj.outLevels((ii-1)/2+1);
                else
                    xp(ii) = x(ii/2);
                    yp(ii) = yp(ii-1);
                end
            end
            xp(end) = obj.maxV;
            yp(end) = obj.outLevels(end);
            plot(xp,yp,'k','LineWidth',1.5), grid on, hold on
            xlabel('Input Level (V)')
            ylabel('Output Level (V)')
            if showBin && obj.Nbits <= 3
                for ii = 1:2^obj.Nbits
                    text(xp(2*ii - 1) + obj.delV./obj.Nbits/8,yp(2*ii)+obj.delV/25,obj.outBin(ii,:))
                end
            end
        end
    end
end
