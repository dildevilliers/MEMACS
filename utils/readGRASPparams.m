function [ParamStruct] = readGRASPparams(GRASPfilePathName,outFilePathName)
% READGRASPPARAMS reads all real_variables from a GRASP.tor file
%
% [ParamStruct] = readGRASPparams(GRASPfilePathName,outFilePathName) reads
% the variables from a GRASP.tor file at GRASPfilePathName and puts all the
% parameters in the structure ParamStruct. The parameters are sorted so
% that the functional ones are evaluated, and constant parameters appear first.
%
% Inputs
% - GRASPfilePathName: full path and name of the .tor file
% - outFilePathName: full path and name of the output .m file (params.m)
%
% Outputs
% - ParamStruct: Vector structure containing the following fields
%   -- descr: Parameter description
%   -- expr:  Parameter expression
%   -- name:  Parameter name
%   -- value: Parameter value
%   -- type:  Parameter type
%
% Dependencies
% -
%
% Created: 2023-03-08, Dirk de Villiers
% Updated: 2023-03-08, Dirk de Villiers
%
% Tested : Matlab R2022b
%  Level : 1
%   File :


if nargin == 0
    [name,path] = uigetfile('*.tor');
    pathName = [path,name];
else
    pathName = GRASPfilePathName;
end
if ~strcmp(pathName(end-3:end),'.tor')
    pathName = [pathName,'.tor'];
end
writeOutput = true;
if nargin < 2 || isempty(outFilePathName)
    outName = 'params';
elseif islogical(outFilePathName) && ~outFilePathName
    outName = [];
    writeOutput = false;
else
    outName = outFilePathName;
end
if writeOutput&& ~strcmp(outName(end-1:end),'.m')
    outName = [outName,'.m'];
end

fid = fopen(pathName,'r');

parCount = 0;
keepReading = true;
while keepReading

    line = fgetl(fid);

    if line == -1
        keepReading = false; 
        break; 
    elseif contains(line,'real_variable')
        parCount = parCount + 1;
        % data = textscan(line, '%s %s');
        % paramName = data{1}{1};
        data = textscan(line, '%s');
        paramName = data{1}{end-1};
        skipLines = true;
        while skipLines
            line = fgetl(fid);  % Dummy past the bracket or comments
            if ismember(line(1),{'(','/'})
                skipLines = true;
            else
                skipLines = false;
            end
        end
        data = textscan(line, '%s %s %s');
        paramValue = str2double(data{3}{1});
        if ~isnan(paramValue)
            paramExpr = paramValue;
            paramType = 0;  % Number
        else
            paramType = 1; % Expression
            paramExpr = data{3}{1};
            paramExpr = paramExpr(2:end-1);  % Strip quotations
            remChar = false(size(paramExpr));
            cc = 1;
            while cc <= numel(paramExpr)
                % Remove the ref() strings
                if cc < (numel(paramExpr) - 3) && strcmp(paramExpr(cc:cc+3),'ref(')
                    remChar(cc:cc+3) = [true,true,true,true];
                    cc = cc+4;
                    nextBrackRem = true;
                elseif paramExpr(cc) == ')' && nextBrackRem
                    remChar(cc) = true;
                    cc = cc + 1;
                    nextBrackRem = false;
                else
                    cc = cc+1;
                end
            end
            paramExpr(remChar) = '';
        end
        ParStruct(parCount).descr = '';  % For compatibility and later upgrade maybe
        ParStruct(parCount).name = paramName;
        ParStruct(parCount).value = paramValue;
        ParStruct(parCount).expr = paramExpr;
        ParStruct(parCount).type = paramType;
    end
end
fclose(fid);

Npar = length(ParStruct);
[~,iSort] = sort([ParStruct.type]);
ParamStruct = ParStruct(iSort);


% Binary sort the rest
for pp = 1:Npar
    numFlag(pp) = isnumeric(ParamStruct(pp).expr);
end
iSplit = find(~numFlag,1);
swopped = true;
whileCounter = 0;
while swopped
    whileCounter = whileCounter + 1;

    swopped = false;
    for iStep = iSplit:Npar-1
        s = strrep(ParamStruct(iStep).expr,' ','');
        s = split(s,{'(',')','+','*','-','/',','});
        mustMove = false;
        for ss = 1:length(s)
            mustMove = mustMove || any(strcmp(s(ss),{ParamStruct(iStep:end).name}));
        end
        if ~mustMove
            swopped = swopped || false;
        else
            swopped = true;
            temp = ParamStruct(iStep+1);
            ParamStruct(iStep+1) = ParamStruct(iStep);
            ParamStruct(iStep) = temp;
        end
    end
end

% Save the m-file, and calculate the unknown values
if writeOutput, fid = fopen(outName,'wt'); end
for pp = 1:Npar
    if isnumeric(ParamStruct(pp).expr)
        val = num2str(ParamStruct(pp).expr);
    else
        val = ParamStruct(pp).expr;
    end
    printLine = [ParamStruct(pp).name,' = ',val,'; %',ParamStruct(pp).descr];
    if writeOutput, fprintf(fid,'%s\n',printLine); end
    
    % Evaluate all variables to get them in the workspace (in order)
    eval(printLine);
    if isnan(ParamStruct(pp).value)
        % Just evcaluate the variable to put it in the value field
        ParamStruct(pp).value = eval(ParamStruct(pp).name);
    end
end
if writeOutput, fclose(fid); end

end