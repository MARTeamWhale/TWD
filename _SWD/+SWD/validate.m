%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "validate" (SWD)
%   Written by WB, based on "beaked_automatic_overlay_combined_det_AMAR250" 
%   by SBP and JS.
%   Last updated Apr. 14, 2023, using MATLAB R2018b
%
%   Description:
%   Prepares and launches the interactive event validation/species ID
%   process for SWD. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function validate(...
    dirPath_root,...
    dirPath_analysis,...
    dirName_detResults,...
    eventMergeOpt,...
    eventMergeVal,...
    nClicksMax,...
    fRange,...
    monID)

    % set subfolder paths
    dirPath_mat = fullfile(dirPath_analysis,'mat');
    dirPath_out = fullfile(dirPath_analysis,dirName_detResults);
    
    % set target name
    targetName = 'Pm';
    
    % get path to 1st MAT file, and extract deployment name and Fs
    [~, matNames] = TWD_Common.Utilities.listFiles(dirPath_mat, 'mat');
    matName1 = matNames{1};
    %%% depName
    depExpr = ['.*(?=_',targetName,'_\d{8}_\d{6})']; % regular expression for extracting 1st part of MAT file name, i.e. deployment
    depName = regexp(matName1,depExpr,'match');
    assert(isscalar(depName),'MAT file names not recognized')
    depName = depName{1};
    %%% Fs
    data1 = load(fullfile(dirPath_mat,matName1),'fs','specClick');
    Fs = data1.fs;
    clear data1
    
    % determine detection protocol by reading text file
    detProtocol = TWD_Common.readProtocolTextFile(dirPath_out);
    dirPath_detCrit = fullfile(dirPath_root,'DetectionCriteria',detProtocol);

    % create click discriminator object
    clickDiscrimFilePath = fullfile(dirPath_detCrit,'ClickDiscrimParams.xlsx');
    clickDiscrimData = struct();
    clickDiscrimData.Obj = TWD_Common.ClickDiscriminator(clickDiscrimFilePath);
    clickDiscrimData.Names = {targetName};
    clickDiscrimData.Defaults = {true, false, false};
    
    % set other variables
    refSpecData = struct.empty();

    % run validation app
    TWD_Common.runValidation(...
        'SWD',...
        dirPath_root,...
        dirPath_analysis,...
        dirName_detResults,...
        depName,...
        targetName,...
        eventMergeOpt,...
        eventMergeVal,...
        nClicksMax,...
        fRange,...
        monID,...
        clickDiscrimData,...
        refSpecData,...
        Fs)
end