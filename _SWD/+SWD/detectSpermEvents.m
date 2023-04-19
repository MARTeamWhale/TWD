%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "detectSpermEvents"
%   Written by Wilfried Beslin
%   Last updated Apr. 14, 2023, using MATLAB R2018b
%
%   Description:
%   Runs event detection in SWD. analyses SWD MAT files and returns a
%   "RawEvents.xlsx" file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES: 
% - this is basically a simplified spermwhale-centric version of detectBeaked

function detectSpermEvents(dirPath_root,dirPath_analysis,detProtocol,dirName_detResults,depName)

    % get path to detection criteria directory
    dirPath_SpermDetCrit = fullfile(dirPath_root,'DetectionCriteria',detProtocol);

    % create EventDetector object
    spermEventCritFilePath = fullfile(dirPath_SpermDetCrit,'EventDetParams.xlsx');
    eDetObj = TWD_Common.EventDetector(spermEventCritFilePath);
    
    % Query MAT files
    dirPath_mat = fullfile(dirPath_analysis,'mat');
    [~, matFileNames] = TWD_Common.Utilities.listFiles(dirPath_mat, 'mat');
    
    % get unique recording files
    matExpr = '.*(?=_\d*\.mat)';
    recNames = regexp(matFileNames,matExpr,'match');
    recNames = vertcat(recNames{:});
    [recNames,~,iMatRec] = unique(recNames);
    nRecs = numel(recNames);
    
    % initialize container for output parts
    segRangesCell = cell(nRecs,1);
    nClicksCell = cell(nRecs,1);
    isEventCell = cell(nRecs,1);
    
    % loop through each recording
    for ii = 1:nRecs
        matsii = iMatRec == ii;
        matFileNamesii = matFileNames(matsii);
        
         % get event detection results for this recording
        fprintf('Recording %d/%d: Detecting sperm whales\n', ii,nRecs)
        [segRangesii,nClicksii,isEventii] = detectSpermEvent_Rec(matFileNamesii,dirPath_mat,eDetObj);
        
        % store results
        segRangesCell{ii} = segRangesii;
        nClicksCell{ii} = nClicksii;
        isEventCell{ii} = isEventii;
    end
    
    % concatenate results from all recordings
    segRanges = vertcat(segRangesCell{:});
    nClicks = vertcat(nClicksCell{:});
    isEvent = vertcat(isEventCell{:});
    
    % convert segment datetimes to strings
    dtFormat = 'yyyyMMdd_HHmmss';
    segRanges.Format = dtFormat;
    segRangeStr = cellstr(segRanges);
    
    % create results directory 
    % (use an alternative name if results exist already)
    dirName_detResultsFinal = TWD_Common.Utilities.uniqueName(dirPath_analysis,dirName_detResults);
    if ~strcmp(dirName_detResults,dirName_detResultsFinal)
        warning('Folder "%s" already exists;\ncreating alternative "%s"',dirName_detResults,dirName_detResultsFinal);
    end
    dirPath_out = fullfile(dirPath_analysis,dirName_detResultsFinal);
    mkdir(dirPath_out);
    
    % create and save spreadsheet
    outTable = table(...
        segRangeStr(isEvent,1),...
        segRangeStr(isEvent,2),...
        nClicks(isEvent),...
        'VariableNames',{'StartTime','EndTime','NClicks'});
    xlsName = sprintf('%s_Pm_RawEvents.xlsx',depName);
    xlsPath = fullfile(dirPath_out,xlsName);
    writetable(outTable,xlsPath);
    
    % add a plain text file marking the name of the detection protocol
    TWD_Common.writeProtocolTextFile(dirPath_out, detProtocol)
end

%% detectSpermEvent_Rec ---------------------------------------------------
function [recRange,nClicks,isEvent] = detectSpermEvent_Rec(matFileNames,dirPath_mat,eDetObj)
% Returns event detection results for one recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get number of MAT files and detection targets
    nMat = numel(matFileNames);
    
    % loop through each MAT file and compile sperm whale click times
    tClickCell = cell(nMat,1);
    for ii = 1:nMat
        % load file
        matFileNameii = matFileNames{ii};
        matFilePathii = fullfile(dirPath_mat,matFileNameii);
        dataii = load(matFilePathii);
        
        % get click start times
        tClickii = seconds(dataii.pos(:,1));
        
        % save loop output to containers
        tClickCell{ii} = tClickii;
    end
    
    % concatenate MAT file data
    tClick = vertcat(tClickCell{:});
    nClicks = numel(tClick);
    
    % get time range
    %%% remember for SWD, 1 segment = 1+ MAT file(s). So here Rec really 
    %%% refers to Segment.
    dtRecStart = datetime(dataii.rawStart);
    dtRecEnd = dtRecStart + seconds(dataii.rawDur);
    recRange = [dtRecStart,dtRecEnd];
    
    % do event detection
    isSperm = true(nClicks,1);
    isEvent = eDetObj.evalTargetEvent_logical(isSperm,tClick);
end