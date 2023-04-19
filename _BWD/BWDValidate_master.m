%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% "BWDValidate_master"
%   Written by Wilfried Beslin
%   For TWD version 1.3
%   Revised 2023/04/19
%
%   Despcription:
%   Interface script for running WB's version of the interactive beaked 
%   whale click validation/species identification procedure. Edit the 
%   variables in this file as needed, then run to launch the species ID
%   app.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to analysis folder
dirPath_analysis = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\Testing_BWD_1month';

% name of results folder
dirName_detResults = 'results02';

% name of detection target
targetName = 'Beaked';

% event merging option
%%% only used once, i.e. if validation spreadsheet doesn't exist yet
eventMergeOpt = 'timebin'; %%% Options are 'none', 'timegap', 'calendar', or 'timebin'
eventMergeVal = duration(0,1,0); %%% MATLAB duration object (for timegap and timebin) or string (for calendar)

% max number of clicks to display
nClicksMax = 10000;

% default frequency display range
fRange = [0, 100];

% ID of monitor on which to show figures
%%% as shown in Windows "Display settings"
monID = 2;

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% determine path to code directory
dirPath_root = mfilename('fullpath');
[dirPath_root,~,~] = fileparts(dirPath_root);

% run process
BWD.validate(...
    dirPath_root,...
    dirPath_analysis,...
    dirName_detResults,...
    targetName,...
    eventMergeOpt,...
    eventMergeVal,...
    nClicksMax,...
    fRange,...
    monID);