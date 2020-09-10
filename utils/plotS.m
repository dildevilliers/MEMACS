function [] = plotS(S,freq)

% Function to plot the S-parameters in the standard S-parameter block 'S'
% Optional frequency vector can be included in 'freq'
% S should be a nxnxf matrix, with n the number of ports

% Author: Dirk de Villiers
% Date  : 2015/12/18

[n1,n2,nf] = size(S);

if n1 ~= n2, error('S should be square for along the 3rd dimension'); end

if nargin == 1, freq = 1:nf; end

for nn1 = 1:n1
    for nn2 = 1:n2
        Snn = reshape(S(nn1,nn2,:),1,nf);
        figure
        plot(freq,dB20(Snn),'k'), grid on
        xlabel('Frequency')
        ylabel(['S_{',num2str(nn1),num2str(nn2),'} (dB)'])
    end
end
