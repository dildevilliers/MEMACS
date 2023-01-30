classdef ArrayDBE
    % Class for the ArraySystem digital back-end
    
    properties
        arraySystem(1,1) ArraySystem
        calVect = 1
    end
    
    methods
        % Very simple constructor
        function obj = ArrayDBE(arraySystem)
            if nargin >= 1, obj.arraySystem = arraySystem; end
            
        end
        
        
        function [X,freq] = getPortFFT(obj,x,overSampFact)
            % x is a typical output from ArraySystem.getPortSignal
            % A signal matrix of size [Nant, Nsamp]
            
            if nargin < 3
                overSampFact = 1;
            end
            
            [Nant,Nt] = size(x);
            Lfft = overSampFact*Nt;
            S = zeros(Nant,Lfft);
            X = zeros(Nant,Lfft/2+1);
            for ii = 1:Nant
                S(ii,:) = fft(x(ii,:),Lfft);
                % Single sided
                Xa = abs(S(ii,:)./Lfft);
                X(ii,:) = Xa(1:Lfft/2+1);
                X(ii,2:end-1) = 2*X(ii,2:end-1);
            end
            freq = obj.arraySystem.freqSamp*(0:(Lfft/2))/Lfft;
        end
        
        function [P] = scanBeam(obj,freqRF,th,ph,x)
            % freqRF is the nominal RF frequency of the system
            % th and ph are the scanning directions in spherical
            % coordinates from the origin. one can be scalar.
            % x is a typical output from ArraySystem.getPortSignal
            % A signal matrix of size [Nant, Nsamp]
            % calVect is a vecotr of calibration phasors of size [1, Nant]
            
            % stretch ph and th if required
            th = (ph+eps(realmin))./(ph+eps(realmin)).*th;
            ph = (th+eps(realmin))./(th+eps(realmin)).*ph;
            
            [u,v,w] = PhTh2DirCos(ph,th);
            k_scan = 2*pi.*freqRF./physconst('lightspeed').*[u;v;w]; % [3 x Nscan]
            A_scan = exp(1i*obj.arraySystem.antPos.pointMatrix.'*k_scan); % [Nant x Nscan]
            AC_scan = bsxfun(@times,A_scan,obj.calVect(:));
            y = AC_scan.'*x;  % [Nscan x Nsamp]
            P = portPower(y); % [Nscan x 1]
        end
        
         function obj = calLinScan(obj,x)
            % Run this to calculate the calibration vector
            for aa = 1:size(x,1)
                obj.calVect(aa) = mean(x(1,:)./x(aa,:));
            end
            obj.arraySys = obj.arraySys.setChannelPhasors(obj.calVect);
         end
        
        %% Plotting
        function plotPortPSD(obj,x,portNumber,overSampFact)
            if nargin < 4
                overSampFact = 1;
            end
            if nargin < 3
                portNumber = 1;
            end
            [X,freq] = getPortFFT(obj,x,overSampFact);
            plot(freq + obj.arraySystem.freqLO,abs(X(portNumber,:)),'k'), grid on, hold on
            xlabel('Frequency [Hz]')
        end
        
        function plotScanBeam(obj,freqRF,th,ph,x,calVect)
            % Much to do here - this is a skelaton
            if nargin < 6
                calVect = 1;
            end
            P = scanBeam(obj,freqRF,th,ph,x,calVect);
            plot(rad2deg(th),P), grid on
            xlabel('\theta^\circ')
        end
        
        function plotCalVect(obj)
            colorVect = 'krbm';
            markerVect = '*osd';
            for pp = 1:size(obj.calVect)
                plot(real(obj.calVect(pp)),imag(obj.calVect(pp)),[colorVect(pp),markerVect(pp)]), hold on, grid on
            end
            t = linspace(0,2*pi,1000);
            x = sin(t) + 1i*cos(t);
            plot(real(x),imag(x));
        end
    end
end



