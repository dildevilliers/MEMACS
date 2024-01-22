function [S, freq] = touchread(file, NrPorts)
%
% TOUCHREAD(...) loads the measured s-parameter data from a text file with Touchstone
% format.  FileName is a text variable containing the full path and filename of 
% the file to be loaded, NrPorts is the number of ports of the S-parameter matrix.
% NrPorts is optional and will be mostly read from the file extension of metadata (if available)
%
% format: [S, freq] = touchread(FileName, NrPorts)
%
% The multidimensional matrix S contains the loaded s-parameters in complex format.  The full
% S-parameter matrix for each frequency point is returned as S(Srows,Scolumns,freq).  
% The variable freq contains the measurement frequencies in Hz.
%
% Author: Robert Lehmensiek and Dirk de Villiers
% Date  : 2017/11
% Updated: 2023/11


if nargin < 1, error('touchread -> Not enough input parameters'); end
if nargin < 2, NrPorts = 0; end

[freq,S] = read(file,NrPorts);


function [f,Smat] = read(file,NrPorts)

ind = findstr(file,'.');
if isempty(ind), error('Add extension'); end

fid = fopen(file, 'r');
if fid==-1, error(['Load error -> Cannot open: ', file]); end

S = fscanf(fid, '%c');
fclose(fid);

% Get all end-of-line indexes
in2 = strfind(S,char([10]));     % for DOS: [13 10] 13 = carriage return; 10 = line feed for DOS (remember +2 in ln = )
% for UNIX: [10] (remember +1 in ln = )

[~,~,ext] = fileparts(file);
n = sscanf(lower(ext),'.s%dp');
if isempty(n)
    % Try to find from the header
    % This works for some TICRA files - add more when we find other versions
    in1 = strfind(S,'[Number of Ports]');

    ind = find(in2>in1); ind = ind(1);
    ln = S(in1+17:in2(ind));

    n = sscanf(ln,'%f');

    if isempty(n), error('Wrong extension, or provide number of ports'); end
end
if NrPorts ~= 0, assert(NrPorts == n,'Input number of ports different from file extension - please check'); end


% Read the # info
in1 = strfind(S,'#');

ind = find(in2>in1); ind = ind(1);
ln = S(in1+1:in2(ind));

if contains(upper(ln),'GHZ'), fu = 1e9;
elseif contains(upper(ln),'MHZ'), fu = 1e6;
elseif contains(upper(ln),'KHZ'), fu = 1e3;
elseif contains(upper(ln),'HZ'), fu = 1;
else, error('Wrong Unit')
end

if contains(upper(ln),'MA'), tp = 1;
elseif contains(upper(ln),'RI'), tp = 2;
elseif contains(upper(ln),'DB'), tp = 3;
else, error('Wrong Type')
end

% Find the end of the header - can be last ! or #
in1 = sort([strfind(S,'!'),strfind(S,'#')]);
in1_row = strfind(S,'! row');    % Find the end of row comment positions and remove them from the vector
in1 = setdiff(in1, in1_row);
in1 = find(in2>in1(end)); in1 = in2(in1(1));
% Now remove the row line comments from the actual data
if ~isempty(in1_row)
    diffInd = in2 - repmat(in1_row(:),1,length(in2));  % Find the number of indexes between the line ends and start of row comments
    diffInd(diffInd < 0) = nan;  % Throw away negative values
    nrChars = min(diffInd,[],2);  % Each row here is the number of characters between the start of the row index and end of line (should all be the same)
    clear diffInd
    assert(max(abs(nrChars - max(nrChars))) < eps,'Error - unknown file format where row comments are not consistantly spaced before line ends')  % Make sure...
    % Make matrix with start of comment until one before line end
    in3 = in1_row(:) + [0:max(nrChars-1)];
    % Sort and reshape
    in3 = in3.';
    in3 = in3(:).';
    S(in3) = ' ';
end


S = sscanf(S(in1+1:end),'%f');
ind_freq = 1:2*n^2+1:length(S);

f = S(ind_freq).*fu;

ind = setdiff(1:length(S),ind_freq);
S = reshape(S(ind),2*n^2,[]).';

switch tp
  case 1
    S = S(:,[1:2:end]) .* exp(1i.*deg2rad(S(:,[2:2:end]))); 
    
  case 2
    S = S(:,[1:2:end]) + 1i.*S(:,[2:2:end]); 
    
  case 3
    S = lin20(S(:,[1:2:end])) .* exp(1i.*deg2rad(S(:,[2:2:end]))); 
    
end

if n==2
  S = [S(:,1) S(:,3) S(:,2) S(:,4)];
end

% Build up the S-matrix
Smat = zeros(n,n,length(f));
for rr = 1:n
    for cc = 1:n
        Smat(rr,cc,:) = S(:,(rr - 1).*n+cc);
    end
end


