function [S,freq,Z0] = touchread(FileName,NrPorts,Fmin,Fmax)

% TOUCHREAD(...) loads the measured s-parameter data from a text file with Touchstone
% format.  FileName is a text variable containing the full path and filename of 
% the file to be loaded, NrPorts is the number of ports of the S-parameter matrix,
% Fmin is the minimum frequency of measurement points that must be considered and Fmax
% is the maximum frequency for which data should be loaded from file.  
% Any s-parameters measured at frequencies below Fmin and
% above Fmax will not be returned to the user.  Fmin and Fmax are optional parameters
% if they are neglected from the parameter list, then values of Fmin = 1 Hz and 
% Fmax = 500 GHz will be set as the default values.
%
% format: [S,freq,Z0]=touchread(FileName,NrPorts,Fmin,Fmax)
%
% The multidimensional matrix S contains the loaded s-parameters in complex format.  The full
% S-parameter matrix for each frequency point is returned as S(Srows,Scolomns,freq).  
% The variable freq contains the measurement frequencies in Hz.
%
% Author: Dirk de Villiers
% Date  : 2008/06/07


% Set the default values of the optional parameters

if (nargin < 4), Fmax=500E9; end
if (nargin < 3), Fmin=1; end

% Read the header of the Touchstone file by checking for comments
% and the unit identifier line.  If the unit ID line is read, 
% scan the contents to determine the format in which the two port
% parameters are stored, as well as the unit of the measured
% frequencies.

nl = 1;
textline_count = 0;
FUnit = 0;
fid = fopen(FileName);
while nl == 1,
    line = fgetl(fid);
    textline_count = textline_count + 1;
    nl = 0;
    if strcmp(line(1),'!') == 1, nl = 1; end;
    if strcmp(line(1),'#') == 1
        nl = 1;
        if (findstr(line,'HZ')  > 0), FUnit = 1E0; end
        if (findstr(line,'Hz')  > 0), FUnit = 1E0; end
        if (findstr(line,'GHZ') > 0), FUnit = 1E9; end
        if (findstr(line,'GHz') > 0), FUnit = 1E9; end
        if (findstr(line,'MHZ') > 0), FUnit = 1E6; end
        if (findstr(line,'MHz') > 0), FUnit = 1E6; end
        if (findstr(line,'KHZ') > 0), FUnit = 1E3; end
        if (findstr(line,'KHz') > 0), FUnit = 1E3; end

        if findstr(line,'MA') > 0, DatFor=1; end
        if findstr(line,'RI') > 0, DatFor=2; end
        infoLineNr = textline_count;
    end
end
status = fseek(fid, 0, 'bof');
form= '%s %s %s %s %s %f';
snp_info = textscan(fid, form, 1, 'HeaderLines', infoLineNr-1);
Z0 = snp_info{6};
% Close the file
fclose('all');

% Get the full matrix from the file
Sfmat = dlmread(FileName,'',textline_count-1,0);

% Sort the matrix
[M,N] = size(Sfmat);
if N > 9
    Sfmat = Sfmat(:,1:9);
end
[M,N] = size(Sfmat);
lines_port = ceil(NrPorts/4);
if NrPorts == 2
    lines_freq = 1;
else
    lines_freq = lines_port*NrPorts;
end
NrFreqs = M/lines_freq;

Sfmat_col1 = Sfmat(:,1);
freq = Sfmat_col1(1:lines_freq:M).*FUnit;

S = zeros(NrPorts,NrPorts,NrFreqs);

if NrPorts == 2
    NrFreqs = M;
    freq = Sfmat(:,1);
    for ii = 1:NrFreqs
        if DatFor == 1
            S(1,1,ii) = Sfmat(ii,2)*exp(j*Sfmat(ii,3)*pi/180);
            S(2,1,ii) = Sfmat(ii,4)*exp(j*Sfmat(ii,5)*pi/180);
            S(1,2,ii) = Sfmat(ii,6)*exp(j*Sfmat(ii,7)*pi/180);
            S(2,2,ii) = Sfmat(ii,8)*exp(j*Sfmat(ii,9)*pi/180);
        elseif DatFor == 2
            S(1,1,ii) = Sfmat(ii,2) + j*Sfmat(ii,3);
            S(2,1,ii) = Sfmat(ii,4) + j*Sfmat(ii,5);
            S(1,2,ii) = Sfmat(ii,6) + j*Sfmat(ii,7);
            S(2,2,ii) = Sfmat(ii,8) + j*Sfmat(ii,9);
        end
    end
else
    for ii = 1:NrFreqs
        for jj = 1:NrPorts
            first_row = ((ii-1)*lines_port*NrPorts) + (lines_port*jj - (lines_port-1));
            if jj == 1
                submat_temp = Sfmat(first_row:first_row+lines_port-1,:);
                submat = [submat_temp(1,2:N);submat_temp(2:lines_port,1:N-1)];
                if lines_port > 1
                    submat_row = reshape(submat',1,lines_port*(N-1));
                else
                    submat_row = submat;
                end
            else
                submat = Sfmat(first_row:first_row+lines_port-1,1:N-1);
                if lines_port > 1
                    submat_row = reshape(submat',1,lines_port*(N-1));
                else
                    submat_row = submat;
                end
            end
            for kk = 1:NrPorts
                if DatFor == 1
                    S_row(kk) = submat_row(2*kk-1)*exp(j*submat_row(2*kk)*pi/180);
                elseif DatFor == 2
                    S_row(kk) = submat_row(2*kk-1) + j*submat_row(2*kk);
                end
            end
            S(jj,:,ii) = S_row;
        end
    end
end


freq_range = find(freq >= Fmin & freq <= Fmax);
freq = freq(freq_range);
S = S(:,:,freq_range);