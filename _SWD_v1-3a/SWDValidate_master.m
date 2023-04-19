%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% "SWDValidate_master"
%   Written by Wilfried Beslin
%   For SWD version 1.3a 
%   Revised 2022/10/13
%
%   Despcription:
%   Interface script for running interactive sperm whale click 
%   validation/species identification procedure. Edit the variables in this 
%   file as needed, then run to launch the species ID app.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to analysis folder
dirPath_analysis = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\Testing_SWD_1-3';

% name of results folder
dirName_detResults = 'results03';

% event merging option
%%% only used once, i.e. if validation spreadsheet doesn't exist yet
eventMergeOpt = 'none'; %%% Options are 'none', 'timegap', 'calendar', or 'timebin'
eventMergeVal = ''; %%% MATLAB duration object (for timegap and timebin) or string (for calendar)

% max number of clicks to display
nClicksMax = 10000;

% default frequency display range
fRange = [0, 50];

% ID of monitor on which to show figures
%%% as shown in Windows "Display settings"
monID = 2;

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% determine path to code directory
dirPath_root = mfilename('fullpath');
[dirPath_root,~,~] = fileparts(dirPath_root);

% run process
SWD.validate(...
    dirPath_root,...
    dirPath_analysis,...
    dirName_detResults,...
    eventMergeOpt,...
    eventMergeVal,...
    nClicksMax,...
    fRange,...
    monID);