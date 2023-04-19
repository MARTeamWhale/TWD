%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "replaceValidations"
%   Original written by Joy Stanistreet
%   Edited by Wilfried Beslin
%   Last Updated Oct. 18, 2019
%
%   Description:
%   Replace original validation codes with new ones for a subset of 
%   re-validated events.
%   Comments will also be replaced if they have been edited.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to results folder containing spreadsheets
dirPath = 'C:\Users\BeslinW\Documents\BeakedWhales_Analyses\BWD_TEST\results';

% Name of original validation spreadsheet
fileName_Master = 'TEST_Beaked_Validated_Full.xlsx';

% Name of revised subset spreadsheet
fileName_Revised = 'TEST_Beaked_Validated_Sub1.xlsx';

% Name to give to new combined spreadsheet
fileName_Out = 'TEST_Beaked_Validated_final.xlsx';

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set full paths
filePath_Master = fullfile(dirPath,fileName_Master);
filePath_Revised = fullfile(dirPath,fileName_Revised);
filePath_Out = fullfile(dirPath,fileName_Out);

% load master and revised spreadsheets
dataMaster = readtable(filePath_Master);
dataRevised = readtable(filePath_Revised);

% match datasets based on start time
% find which rows in dataMaster match rows in dataRevised
[~,iMaster,iRevised] = intersect(dataMaster.StartTime,dataRevised.StartTime);

% create new table
dataOut = dataMaster;
dataOut.Species(iMaster) = dataRevised.Species(iRevised);
dataOut.Comment(iMaster) = dataRevised.Comment(iRevised);

% save new master spreadsheet
writetable(dataOut,filePath_Out)
