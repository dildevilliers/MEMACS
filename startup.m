% STARTUP Sets paths for MEMACS package
% [] = startup() Sets the current path as the root for the package, and
% links all internal folders
%
% Inputs
% - 
%
% Outputs
%
% Dependencies
% -
%
% Created: 2019-05-10, Dirk de Villiers
% Updated: 2019-05-10, Dirk de Villiers
%
% Tested : Matlab R2018b, Dirk de Villiers
%  Level : 1
%   File : 
%

function startup

% get location of this file (toolbox root path)
p = mfilename('fullpath');
MEMACSroot = p(1:end-7);

% add matlab class paths
addpath(genpath([MEMACSroot,'Classes']))
addpath(genpath([MEMACSroot,'data']))
addpath(genpath([MEMACSroot,'developerScratch']))
addpath(genpath([MEMACSroot,'testScripts']))
addpath(genpath([MEMACSroot,'utils']))
addpath(genpath([MEMACSroot,'contrib']))
end
