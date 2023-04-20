%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "replaceValidations"
%   Original written by Joy Stanistreet
%   Edited by Wilfried Beslin
%   Last updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   Replace original validation codes with new ones for a subset of 
%   re-validated events.
%   Comments will also be replaced if they have been edited.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUT - CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% path to results folder containing spreadsheets
dirPath = 'D:\TWD_results\results01';

% Name of original validation spreadsheet
fileName_Master = 'DEP_Target_Validated.xlsx';

% Name of revised subset spreadsheet
fileName_Revised = 'DEP_Target_Validated_Sub1.xlsx';

% Name to give to new combined spreadsheet
fileName_Out = 'DEP_Target_Validated_Final.xlsx';

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


newValidations = do_replaceValidations(dirPath, fileName_Master, fileName_Revised, fileName_Out);


%% ------------------------------------------------------------------------
function dataOut = do_replaceValidations(dirPath, fileName_Master, fileName_Revised, fileName_Out)
    % set full paths
    filePath_Master = fullfile(dirPath, fileName_Master);
    filePath_Revised = fullfile(dirPath, fileName_Revised);
    filePath_Out = fullfile(dirPath, fileName_Out);

    % load master and revised spreadsheets
    dataMaster = readtable(filePath_Master);
    dataRevised = readtable(filePath_Revised);

    % match datasets based on start time
    %%% find which rows in dataMaster match rows in dataRevised
    [~,iMaster,iRevised] = intersect(dataMaster.StartTime, dataRevised.StartTime);

    % create new table
    dataOut = dataMaster;
    dataOut.Species(iMaster) = dataRevised.Species(iRevised);
    dataOut.Comment(iMaster) = dataRevised.Comment(iRevised);

    % save new master spreadsheet
    writetable(dataOut, filePath_Out)
    fprintf('Saved new validations table to:\n"%s"\n', filePath_Out)
end