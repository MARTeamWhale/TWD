%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "detectBeaked"
%   Written by WB, based on "analyze_now_AMAR250" and 
%   "combine_bw_det_AMAR" by SBP and JS.
%   Last updated Jul. 24, 2019, using MATLAB R2016b
%
%   Description:
%   Runs automatic detection of beaked whale events based on discrimination
%   of clicks. Though this function is based on "analyze_now_AMAR250" and 
%   "combine_bw_det_AMAR", it works completely differently. Click
%   discrimination is no longer performed directly with hardcoded
%   thresholds. Instead, it is managed using an object-oriented structure, 
%   where thresholds are set via external Excel spreadsheets. 
%
%   Furthermore, it is now possible to use multiple sets of discrimination 
%   criteria (each implemented as a ClickDiscriminator object). This can be
%   done on two levels. At the first level (called the "target" level), 
%   each discrimination scheme will result in seperate event detection 
%   results (similar to rerunning the code to try out different schemes, 
%   except you don't have to do it manually). At the second level
%   (within-target level), each discriminator contributes to the decision
%   on a particular click. In other words, a click will "pass" if it meets
%   at least one set of criteria.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
%   - May be desirable to implement automatic resumption here too, because 
%   the routine can take some time, and there is no information saved until 
%   it finishes completely. Thus interruption results in total loss of 
%   progress.
% - To do this, will need to update spreadsheet(s) periodically


function detectBeaked(dirPath_root,depName,dirPath_analysis,detProtocol,dirName_detResults)

    % get path to detection criteria directory
    dirPath_detCrit = fullfile(dirPath_root,'DetectionCriteria',detProtocol);

    % create beaked whale detector objects
    disp('Preparing for beaked whale detection')
    detTargets = createDetectors(dirPath_detCrit);
    nTargets = numel(detTargets);
    
    % Query MAT files
    dirPath_mat = fullfile(dirPath_analysis,'mat');
    matInfo = dir([dirPath_mat,filesep,'*.mat']);
    matFileNames = {matInfo.name}';
    
    % get unique recording files
    matExpr = '.*(?=_\d*\.mat)';
    recNames = regexp(matFileNames,matExpr,'match');
    recNames = vertcat(recNames{:});
    [recNames,~,iMatRec] = unique(recNames);
    nRecs = numel(recNames);
    
    % initialize container for output parts
    segRangesCell = cell(nRecs,1);
    nClicksAllCell = cell(nRecs,1);
    nClicksBeakedCell = cell(nRecs,1);
    isEventCell = cell(nRecs,1);
    
    % loop through each recording
    for ii = 1:nRecs
        matsii = iMatRec == ii;
        matFileNamesii = matFileNames(matsii);
        
        % get event detection results for this recording
        fprintf('Recording %d/%d: Detecting beaked whales\n', ii,nRecs)
        [segRangesii,nClicksAllii,nClicksBeakedii,isEventii] = detectBeaked_Rec(matFileNamesii,dirPath_mat,detTargets);
        
        % store results
        segRangesCell{ii} = segRangesii;
        nClicksAllCell{ii} = nClicksAllii;
        nClicksBeakedCell{ii} = nClicksBeakedii;
        isEventCell{ii} = isEventii;
    end
    
    % concatenate results from all recordings
    segRanges = vertcat(segRangesCell{:});
    nClicksAll = vertcat(nClicksAllCell{:});
    nClicksBeaked = vertcat(nClicksBeakedCell{:});
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
    
    % for each target, create tables and save them as 
    % spreadsheets
    for ii = 1:nTargets
        isEventii = isEvent(:,ii);
        
        % create table
        outTableii = table(...
            segRangeStr(isEventii,1),...
            segRangeStr(isEventii,2),...
            nClicksAll(isEventii),...
            nClicksBeaked(isEventii,ii),...
            'VariableNames',{'StartTime','EndTime','NClicksAll','NClicksBeaked_EventDet'});
        
        % save spreadsheet
        targetNameii = detTargets(ii).name;
        fprintf('Saving all events for target "%s"\n',targetNameii)
        
        xlsNameii = sprintf('%s_%s_RawEvents.xlsx',depName,targetNameii);
        xlsPathii = fullfile(dirPath_out,xlsNameii);
        writetable(outTableii,xlsPathii);
    end
    
    % add a plain text file marking the name of the detection protocol
    txtFileName = 'DetectionProtocol.txt';
    txtFilePath = fullfile(dirPath_out,txtFileName);
    [txtFileID,txtFileError] = fopen(txtFilePath,'w');
    if txtFileID ~= -1
        fprintf(txtFileID,'%s',detProtocol);
        fclose(txtFileID);
    else
        warning('Could not write detection protocol file (%s):\n%s',txtFileName,txtFileError)
    end
end

%% createDetectors --------------------------------------------------------
function detTargets = createDetectors(dirPath_detCrit)
% Creates event detector and click discriminator objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define file naming conventions
    eventCritFileFlag = 'EventDetParams';
    clickCritFileFlag = 'ClickDiscrimParams_EventDet';
    
    % get all subdirectories in detection criteria folder
    targetNames = TWD_Common.Utilities.getSubDirs(dirPath_detCrit);
    nTargets = numel(targetNames);
    
    % initialize output (detection target struct)
    detTargets(nTargets,1) = struct(...
        'name',[],...
        'eventDet',[],...
        'clickDiscrim',[]);
    
    % loop through each target
    for ii = 1:nTargets
        targetNameii = targetNames{ii};
        targetPathii = fullfile(dirPath_detCrit,targetNameii);
        
        % get all spreadsheets in target directory
        fileNamesii = TWD_Common.Utilities.getFileNames(targetPathii,'xlsx');
        
        % create event detector
        iEventDetii = find(contains(fileNamesii,eventCritFileFlag));
        assert(isscalar(iEventDetii),'In "%s":\nExpected exactly one Excel file containing the name "%s"',targetPathii,eventCritFileFlag)
        eventDetFileii = fileNamesii{iEventDetii};
        eventDetFilePathii = fullfile(targetPathii,eventDetFileii);
        eventDetii = TWD_Common.EventDetector(eventDetFilePathii);
        
        % create click discriminator(s)
        iClickDiscrimii = find(contains(fileNamesii,clickCritFileFlag));
        assert(isvector(iClickDiscrimii),'In "%s":\nExpected at least one Excel file containing the name "%s"',targetPathii,clickCritFileFlag)
        nClickDiscrimii = numel(iClickDiscrimii);
        clickDiscrimCellii = cell(nClickDiscrimii,1);
        for jj = 1:nClickDiscrimii
            iClickDiscrimjj = iClickDiscrimii(jj);
            clickDiscrimFilejj = fileNamesii{iClickDiscrimjj};
            clickDiscrimFilePathjj = fullfile(targetPathii,clickDiscrimFilejj);
            clickDiscrimjj = BWD.ClickDiscriminator(clickDiscrimFilePathjj);
            clickDiscrimCellii{jj} = clickDiscrimjj;
        end
        clickDiscrimii = vertcat(clickDiscrimCellii{:});
        
        % save iteration output to struct
        detTargets(ii).name = targetNameii;
        detTargets(ii).eventDet = eventDetii;
        detTargets(ii).clickDiscrim = clickDiscrimii;
    end
end

%% detectBeaked_Rec -------------------------------------------------------
function [segRanges,nClicksAll,nClicksBeaked,isEvent] = detectBeaked_Rec(matFileNames,dirPath_mat,detTargets)
% Returns event detection results for one recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get number of MAT files and detection targets
    nMat = numel(matFileNames);
    nTargets = numel(detTargets);
    
    % loop through each MAT file and compile beaked whale clicks
    isBeakedCell = cell(nMat,1);
    tClickCell = cell(nMat,1);
    for ii = 1:nMat
        % load file
        matFileNameii = matFileNames{ii};
        matFilePathii = fullfile(dirPath_mat,matFileNameii);
        dataii = load(matFilePathii);
        
        % get click start times
        tClickii = seconds(dataii.pos(:,1));
        nClicksii = numel(tClickii);
        
        % loop through each discrimination scheme and detect beaked whale
        % clicks
        isBeakedii = zeros(nClicksii,nTargets);
        for jj = 1:nTargets
            clickDiscrimjj = detTargets(jj).clickDiscrim;
            nDiscrimjj = numel(clickDiscrimjj);
            
            % loop through each discriminator. Clicks are considered beaked
            % whale if they pass criteria for at least one discriminator.
            isBeakedMatjj = false(nClicksii,nDiscrimjj);
            for kk = 1:nDiscrimjj
                clickDiscrimkk = clickDiscrimjj(kk);
                isBeakedMatjj(:,kk) = clickDiscrimkk.evalBeakedClicks(dataii);
            end
            isBeakedjj = any(isBeakedMatjj,2);
            isBeakedii(:,jj) = isBeakedjj;
        end
        
        % save loop output to containers
        tClickCell{ii} = tClickii;
        isBeakedCell{ii} = isBeakedii;
    end
    
    % concatenate MAT file data
    tClick = vertcat(tClickCell{:});
    isBeaked = vertcat(isBeakedCell{:});
    
    % get time segments
    dtSegStart = datetime(dataii.rawStart);
    dtSegEnd = dtSegStart + seconds(dataii.rawDur);
    segRanges = [dtSegStart,dtSegEnd];
    nSeg = numel(dtSegStart);
    dtRecStart = dtSegStart(1);
    
    % loop through each segment and compute total number of clicks, number
    % of beaked whale clicks, and whether an event exists or not
    nClicksAll = zeros(nSeg,1);
    nClicksBeaked = zeros(nSeg,nTargets);
    isEvent = false(nSeg,nTargets);
    for ii = 1:nSeg
        tSegStartii = dtSegStart(ii) - dtRecStart;
        tSegEndii = dtSegEnd(ii) - dtRecStart;
        
        % isolate clicks in segment
        isInSegii = (tClick >= tSegStartii) & (tClick < tSegEndii);
        isBeakedii = isBeaked(isInSegii,:);
        tClickii = tClick(isInSegii);
        
        % get total number of clicks and number of beaked whale clicks
        nClicksAllii = sum(isInSegii);
        nClicksBeakedii = sum(isBeakedii,1);
        
        % loop through each target and determine whether or not an event 
        % was detected
        isEventii = false(1,nTargets);
        for jj = 1:nTargets
            eventDetjj = detTargets(jj).eventDet;
            isBeakedjj = isBeakedii(:,jj);
            isEventjj = eventDetjj.evalBeakedEvent_logical(isBeakedjj,tClickii);
            isEventii(jj) = isEventjj;
        end
        
        % save iteration output to container
        nClicksAll(ii) = nClicksAllii;
        nClicksBeaked(ii,:) = nClicksBeakedii;
        isEvent(ii,:) = isEventii;
    end
end