function [ParamStruct] = readCSTparams(CSTfilePathName,outFilePathName)
% READCSTPARAMS reads a CST parameters.json file
% 
% [ParamStruct] = readCSTparams(CSTfilePathName,outFilePathName) reads the
% .json file at CSTfilePathName and puts all the parameters in the
% structure ParamStruct. The parameters are sorted so that they can be
% functional ones are evaluated, and constant parameters appear first.
%
% Inputs
% - CSTfilePathName: full path and name of the .json file
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
% Created: 2019-10-01, Dirk de Villiers
% Updated: 2019-10-01, Dirk de Villiers
%
% Tested : Matlab R2018b
%  Level : 1
%   File : 
%
% Example
%   F = FarField;
%   F.plot

if nargin == 0
    [name,path] = uigetfile('*.json');
    pathName = [path,name];
else
    pathName = CSTfilePathName;
end
if ~strcmp(pathName(end-4:end),'.json')
    pathName = [pathName,'.json'];
end
writeOutput = true;
if nargin < 1 || isempty(outFilePathName)
    outName = 'params';
elseif islogical(outFilePathName) && ~outFilePathName
    outName = [];
    writeOutput = false;
else
    outName = outFilePathName;
end
if writeOutput && ~strcmp(outName(end-1:end),'.m')
    outName = [outName,'.m'];
end

%% Read the parameter file
fid = fopen(pathName);
keepReading = true;
% Skip first 2 lines
fgetl(fid);     
fgetl(fid);
parCount = 0;
while keepReading
    a1 = fgetl(fid);
    a1Clean = strtrim(a1);
    b1 = strsplit(a1,'"');    % Break into parts
    
    if strncmp(a1Clean,'{',1)
        parCount = parCount + 1;
        [ParStruct(parCount).descr,ParStruct(parCount).expr,ParStruct(parCount).name,ParStruct(parCount).value] = deal([]);
        readPar = true;
        while readPar
            a = fgetl(fid);
            aClean = strtrim(a);
            b = strsplit(a,'"');    % Break into parts
            if strncmp(aClean,'"descr"',7)
                ParStruct(parCount).descr = b{4};
            elseif strncmp(aClean,'"expr"',6)
                % Fix the CST function names
                expr = b{4};
                expr = strrep(expr,'atn','atan');
                expr = strrep(expr,'sqr','sqrt');
                ParStruct(parCount).expr = expr;
            elseif strncmp(aClean,'"name"',6)
                ParStruct(parCount).name = b{4};
            elseif strncmp(aClean,'"value"',7)
                ParStruct(parCount).value = b{4};
            elseif strncmp(aClean,'}',1)
                readPar = false;
            end
        end
    elseif strncmp(a1Clean,'],',2)
        keepReading = false;
    end
end
fclose(fid);

Npar = length(ParStruct);
for pp = 1:Npar
    ParStruct(pp).value = str2num(ParStruct(pp).value);
    if ~isempty(str2num(ParStruct(pp).expr))
        ParStruct(pp).type = 0;
        ParStruct(pp).expr = str2num(ParStruct(pp).expr);
    else
        ParStruct(pp).type = 1;
    end
end

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

if writeOutput
    fid = fopen(outName,'wt');
    for pp = 1:Npar
        if isnumeric(ParamStruct(pp).expr)
            val = num2str(ParamStruct(pp).expr);
        else
            val = ParamStruct(pp).expr;
        end
        fprintf(fid,'%s\n',[ParamStruct(pp).name,' = ',val,'; %',ParamStruct(pp).descr]);
    end
    fclose(fid);
end