%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "runValidation"
%   Written by WB, based on "beaked_automatic_overlay_combined_det_AMAR250" 
%   by SBP and JS.
%   Last updated Apr. 14, 2023, using MATLAB R2018b
%
%   Description:
%   Runs the interactive event validation/species ID process. This
%   includes both construction and display of event summary, spectral 
%   overlay, and individual click plots.
%   While the original plots and basic workflow of the validation 
%   process have not changed much, this function has been very heavily 
%   modified since "beaked_automatic_overlay_combined_det_AMAR250",
%   offering much more flexibility and some new figures.
%
%   Upon first run on a particular dataset, this function will create a
%   "Validated" spreadsheet that will be updated directly (though saving to
%   disk only occurs on demand or upon closure). When run again on a
%   partially processed dataset, the routine will pick up where it left
%   off.
%
%   NOTE: Do not terminate the process by pressing Ctrl+C - use the "q" 
%   option instead. If it is abruptly terminated, there may remain hidden
%   figures that cannot be closed with a "close all" command (you will get 
%   an error). If this happens, run the "closeGhostWin" script included 
%   with this package, or do the following:
%       - enter "a = gcf"
%       - enter "delete(a)"
%       - repeat until "gcf" causes a new figure to open
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEV NOTES
% 2019/10/11
% - This function is to be used by both BWD and SWD but still retains some
% "Beaked" terminology. Be sure to change variable names as appropriate.
% - Creation and maintenance of discriminator switcher window even in SWD 
% is a waste since discriminator switching is useless in SWD

function runValidation(...
    caller,... % string specifying 'BWD' or 'SWD'
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

    % 1) INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Initializing...')
    
    % set subfolder paths
    dirPath_mat = fullfile(dirPath_analysis,'mat');
    dirPath_results = fullfile(dirPath_analysis,dirName_detResults);
    dirPath_plots = fullfile(dirPath_results,['plots_',targetName]);
    dirPath_plotsClicks = fullfile(dirPath_plots,'clicks');
    
    % get MAT file names and start times
    %%% Remember, multiple MATs may have the same start time.
    [~, matNames] = TWD_Common.Utilities.listFiles(dirPath_mat, 'mat');
    matStarts = TWD_Common.Utilities.readDateTime(matNames);
    
    % Get ID and/or event spreadsheets
    dtFormat = 'yyyyMMdd_HHmmss';
    [xlsOutData,xlsOutPath,iStart,eventsTable,validInput] = getSpreadsheets(dirPath_results,depName,targetName,dtFormat,eventMergeOpt,eventMergeVal,caller);
    %%% if spreadsheet(s) could not be loaded, end program
    if ~validInput
        disp('Exiting program')
        return
    end
    
    % load species ID codes and set up ID prompt string
    spCodeFilePath = fullfile(dirPath_root,'SpeciesCodes.xlsx');
    spCodes = readtable(spCodeFilePath);
    sidPrompt = strcat(num2str(spCodes.Code),{' ='},spCodes.Species);
    sidPrompt = strjoin(sidPrompt,', ');
    sidPrompt = ['SPECIES: ',sidPrompt];
    
    % define variables
    xlsOutVars = xlsOutData.Properties.VariableNames;
    nEvents = height(xlsOutData);
    eStartTimes = datetime(xlsOutData.StartTime,'Format',dtFormat);
    eEndTimes = datetime(xlsOutData.EndTime,'Format',dtFormat);
    targetTimeStr = xlsOutData.(xlsOutVars{contains(xlsOutVars,'TimeWith')});
    eSegStartTimes = datetime(eventsTable.StartTime,'Format',dtFormat);
    eSegEndTimes = datetime(eventsTable.EndTime,'Format',dtFormat);
    %xlsOutPath = xlsPath_IDs;
    nDiscrim = numel(clickDiscrimData.Obj);
    
    % close existing figures
    close all

    % initialize main figure handles
    hfMain = createMainWindows(caller,clickDiscrimData.Obj,monID);
    doOverlayPlot = isfield(hfMain,'Overlay');
    
    % define data for Click Params figure (stored within UserData)
    createClickParamsFigData(hfMain.ClickParams,clickDiscrimData.Obj); %%% consider calling within createMainWindows
    
    % initialize switcher figure and table handles
    [hfDiscrim,htDiscrim] = createDiscrimSwitchWindow(clickDiscrimData.Names,clickDiscrimData.Defaults);
    [hfEvent,htEvent] = createEventSwitchWindow();
    
    % initialize param pass display status
    showParamPass = false;
    
    % create plot directories if they don't exist
    if ~isfolder(dirPath_plots)
        mkdir(dirPath_plots);
    end
    if ~isfolder(dirPath_plotsClicks)
        mkdir(dirPath_plotsClicks);
    end
    
    
    % 2) EXECUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % loop through each event
    try
        ii = iStart;
        doQuit = false;
        while ii <= nEvents
            fprintf('Setting up event %d/%d...\n',ii,nEvents)

            % get event times
            eStartTimeii = eStartTimes(ii);
            eEndTimeii = eEndTimes(ii);
            
            % get duration with beaked/sperm whales
            targetTimeStrii = targetTimeStr{ii};
            targetTimeStrPartsii = strsplit(targetTimeStrii,'_');
            targetTimeNumii = str2double(targetTimeStrPartsii);
            targetTimeii = duration(targetTimeNumii);
            
            % query and run all discriminators on all clicks within event
            [clickStartsAllii,clickMatIdxAllii,clickPassedDiscrimAllii,clickParamValsAllii,clickParamPassAllii,clickQAllii] = inspectAllEventClicks(eStartTimeii,eEndTimeii,eSegStartTimes,eSegEndTimes,matStarts,matNames,dirPath_mat,clickDiscrimData.Obj);
            clickStartsAllii.Format = 'yyyy-MM-dd   HH:mm:ss.SSS'; % for displaying in params window (is this still used?)
            
            %%% CONSIDER SETTING FIGURE NAMES HERE
            
            % reset active discrimination criteria
            htDiscrim.Data = hfDiscrim.UserData.DefaultSettings;

            % process file
            %%% there are multiple possibilities here. Use a state machine 
            %%% (while loop with switch-case block) to determine which 
            %%% processes to execute.
            doNextEvent = false;
            state = 'reset';
            isNewEvent = true;
            while ~doNextEvent
                switch state
                    case 'reset'
                        % get active discriminators and isolate those
                        % clicks that meet the criteria
                        isDiscrimPositive = cell2mat(htDiscrim.Data(:,1));
                        isDiscrimNegative = cell2mat(htDiscrim.Data(:,2));
                        isDiscrimActive = isDiscrimPositive | isDiscrimNegative;
                        if any(isDiscrimPositive)
                            isTarget = any(clickPassedDiscrimAllii(:,isDiscrimPositive),2) & ~any(clickPassedDiscrimAllii(:,isDiscrimNegative),2);
                        else
                            % if no discriminators are positive, show all 
                            % clicks in event, unless they meet criteria
                            % for negative discriminators
                            isTarget = ~any(clickPassedDiscrimAllii(:,isDiscrimNegative),2);
                        end
                        
                        % get subset of clicks if there are too many
                        subset = limitNumClicks(clickStartsAllii,eSegStartTimes,eSegEndTimes,isTarget,clickQAllii,nClicksMax,caller);
                        clickStartTarget = clickStartsAllii(subset);
                        clickMatIdxTarget = clickMatIdxAllii(subset);
                        clickParamValsTarget = clickParamValsAllii(subset,:);
                        clickParamPassTarget = clickParamPassAllii;
                        for jj = 1:nDiscrim
                            clickParamPassTarget{jj} = clickParamPassTarget{jj}(subset,:);
                        end
                        clickPassedDiscrimTarget = clickPassedDiscrimAllii(subset,:);
                        
                        % get data for clicks of interest
                        data = getClicks(clickStartTarget,clickMatIdxTarget,matStarts,matNames,dirPath_mat);      
                        
                        % clear all plots 
                        %%% some apparent redundancy because this is 
                        %%% already performed in 'done' state, however
                        %%% 'reset' does not always follow 'done' (notably
                        %%% on discriminator change)
                        structfun(@clf,hfMain);
                        
                        % do summary plots
                        plotSummaryFigs(hfMain.Summary,data,depName,eStartTimeii,fRange,caller);
                        if doOverlayPlot
                            plotSpectrumOverlay(hfMain.Overlay,data,refSpecData,ii,nEvents,targetTimeii,sum(isTarget));
                        end
                        haConSpec = findobj(hfMain.Summary, 'Type','axes', 'Tag','ConcatenatedSpectrogramAxes');
                        
                        % Save plots and number of clicks if this is first
                        % reset for this event
                        if isNewEvent
                            figNameBaseii = sprintf('%s_%s_%s',depName,char(eStartTimeii),targetName);

                            % summary plots
                            summaryFigPathii = fullfile(dirPath_plots,[figNameBaseii,'_summary']);
                            saveas(hfMain.Summary,summaryFigPathii,'jpg');
                            disp('Summary figure saved')

                            % overlay plot
                            if doOverlayPlot
                                overlayFigPathii = fullfile(dirPath_plots,[figNameBaseii,'_overlay']);
                                saveas(hfMain.Overlay,overlayFigPathii,'jpg');
                                disp('Spectral overlay figure saved')
                            end

                            % add number of target clicks to output (BWD only)
                            if strcmp(caller,'BWD')
                                xlsOutData.NClicksBeaked_Validation(ii) = sum(isTarget); %%% This will take defaults into account (not the case in older versions)
                            end

                            % print message and disable saving flag
                            isNewEvent = false;
                        end
                        
                        % set click data for click viewing
                        FPeakTarget = data.peakFr;
                        ppSignalTarget = data.ppSignal;
                        %%% get indicies of clicks sorted by p-p magnitude and peak frequency
                        [~, ippSort]=sort(ppSignalTarget,'descend');
                        [~, iFPeakSort]=sort(FPeakTarget,'descend');
                        iFPeakSortR = flipud(iFPeakSort);
                        iChronSort = sort(ippSort);
                        %%% other stuff
                        nClicksTarget = sum(subset);
                        clickImg = 1;
                        clickImgStep = 1;
                        clickStepStr = 'Next';
                        clickSortOpt = 'peak-to-peak';
                        click_iSort = ippSort; 

                        % change state
                        state = 'viewclicks';
                        
                    case 'viewclicks'
                        % click viewing and input processing.
                        % State changes depending on input:
                        %   'done' if species ID has been set
                        %   'reset' if discriminator has been changed
                        %   'exit' if user wants to quit
                        if nClicksTarget > 0
                            plotClick(hfMain.Click,hfMain.ClickSeq,data,clickStartTarget,clickSortOpt,click_iSort,clickImg,Fs,fRange,figNameBaseii,caller);
                            showClickParams(hfMain.ClickParams,clickParamValsTarget,clickParamPassTarget,clickPassedDiscrimTarget,clickDiscrimData.Names,isDiscrimActive,clickStartTarget,click_iSort,clickImg,showParamPass,figNameBaseii)
                            highlightConSpecClick(haConSpec, find(iFPeakSortR==click_iSort(clickImg)), [1,0.5,0], 'SpecHighlight', nClicksTarget);
                        end

                        % get user input
                        eventPrompt = sprintf('EVENT %d/%d',ii,nEvents);
                        clickPrompt = sprintf('CLICK PLOTS: enter =%s, d =Change Direction, j =Jump to Click (Number Input), jc =Jump to Click (Select from Plot)\n             o =Change Order, r =Reset, v =Change Freq Range, s =Save',clickStepStr);
                        switch caller
                            case 'BWD'
                                actPrompt = 'OTHER ACTIONS: f =Toggle Param Status, p =Change Discriminator, c =Add Comment, u =Update Spreadsheet, U =New Spreadsheet, e =Change Event, q =Quit';
                            case 'SWD'
                                actPrompt = 'OTHER ACTIONS: f =Toggle Param Status, c =Add Comment, u =Update Spreadsheet, U =New Spreadsheet, e =Change Event, q =Quit';
                        end
                        prompt = sprintf('\n%s\n%s\n%s\n%s\n',eventPrompt,clickPrompt,actPrompt,sidPrompt);
                        opt = input(prompt,'s');
                        switch opt
                            case ''
                                % next/previous plot
                                clickImg = clickImg + clickImgStep;
                                if clickImg > nClicksTarget
                                    clickImg = 1;
                                elseif clickImg < 1
                                    clickImg = nClicksTarget;
                                end

                            case 'd'
                                % change direction, then go to next
                                clickImgStep = clickImgStep*-1;
                                if clickImgStep > 0
                                    clickStepStr = 'Next';
                                else
                                    clickStepStr = 'Previous';
                                end
                                clickImg = clickImg + clickImgStep;
                                if clickImg > nClicksTarget
                                    clickImg = 1;
                                elseif clickImg < 1
                                    clickImg = nClicksTarget;
                                end
                                
                            case 'j'
                                % jump to another click by typing the click
                                % number into the command line (relative to
                                % sorting order)
                                jumpOpt = input('Enter click number: ');
                                if ismember(jumpOpt, click_iSort)
                                    clickImg = jumpOpt;
                                else
                                    disp('Invalid input')
                                end
                                
                            case 'jc'
                                % jump to another click by clicking on the
                                % concatenated spectrogram plot
                                disp('Select a click on the concatenated spectrogram plot')
                                
                                figure(hfMain.Summary);
                                hfMain.Summary.UserData = 'HoverState';
                                hfMain.Summary.WindowButtonMotionFcn = @FigButtonMotionFcn_summary;
                                hfMain.Summary.WindowButtonDownFcn = @FigButtonDownFcn_summary;
                                
                                % wait for callback functions to update the
                                % 'UserData' property of the summary
                                % figure. This property will store the
                                % index of the selected click.
                                waitfor(hfMain.Summary, 'UserData');
                                
                                % determine the correct click to jump to
                                % (if any)
                                if ~isempty(hfMain.Summary.UserData)
                                    
                                    % clicks on the concatenated
                                    % spectrogram are sorted from lowest to
                                    % highest peak frequency. This is the
                                    % opposite of how they are sorted when
                                    % the FPeak sort method is used.
                                    clickImg = find(click_iSort == iFPeakSortR(hfMain.Summary.UserData)); 
                                
                                    % reset UserData property
                                    hfMain.Summary.UserData = [];
                                end
                                
                                % reset callback functions
                                hfMain.Summary.WindowButtonMotionFcn = '';
                                hfMain.Summary.WindowButtonDownFcn = '';

                            case 'o'
                                % change sorting order
                                if all(click_iSort==ippSort)
                                    % switch to pf
                                    click_iSort = iFPeakSort;
                                    clickSortOpt = 'F_{peak}';
                                    clickImg = find(click_iSort == ippSort(clickImg));
                                elseif all(click_iSort==iFPeakSort)
                                    % switch to chronological sorting
                                    click_iSort = iChronSort;
                                    clickSortOpt = 'chronological';
                                    clickImg = find(click_iSort == iFPeakSort(clickImg));
                                else
                                    % switch to pp
                                    click_iSort = ippSort;
                                    clickSortOpt = 'peak-to-peak';
                                    clickImg = find(click_iSort == iChronSort(clickImg));
                                end
                                % reset
                                %clickImg = 1;
                                %clickImgStep = 1;
                                %clickStepStr = 'Next';

                            case 'r'
                                % reset (go to 1st click)
                                clickImg = 1;
                                clickImgStep = 1;
                                clickStepStr = 'Next';

                            case 'v'
                                % Change frequency range
                                freqPrompt = {sprintf('Enter new frequency display range\n(NOTE: event will be reloaded)\n\nFMax:');'FMin:'};
                                fRangeCellstr = flip(strsplit(num2str(fRange))');
                                fRangeCellstrNew = inputdlg(freqPrompt,'Freq Range',1,fRangeCellstr);
                                if ~isempty(fRangeCellstrNew)
                                    fRangeNew = flip(str2double(fRangeCellstrNew))';
                                    if fRangeNew(1) < fRangeNew(2)
                                        fRange = fRangeNew;
                                        state = 'reset';
                                    else
                                        disp('Bad input')
                                    end
                                end

                            case 's'
                                % save image
                                if all(click_iSort==ippSort)
                                    sortStr = 'pp';
                                elseif all(click_iSort==iFPeakSort)
                                    sortStr = 'pf';
                                else
                                    sortStr = 'ch';
                                end
                                clickFigName = sprintf('%s_%s_click_%d_%s', depName, char(eStartTimeii) ,clickImg, sortStr);
                                clickFigPath = fullfile(dirPath_plotsClicks,clickFigName);
                                saveas(hfMain.Click, clickFigPath, 'jpg')
                                disp('Click figure saved')
                                % + sequence
                                clickSeqFigName = [clickFigName,'_sequence'];
                                clickSeqFigPath = fullfile(dirPath_plotsClicks,clickSeqFigName);
                                saveas(hfMain.ClickSeq, clickSeqFigPath, 'jpg')
                                disp('Click Sequence figure saved')
                                % + parameters
                                clickParamFigName = [clickFigName,'_params'];
                                clickParamFigPath = fullfile(dirPath_plotsClicks,clickParamFigName);
                                saveas(hfMain.ClickParams, clickParamFigPath, 'jpg')
                                disp('Click Params figure saved')

                            case 'f'
                                % toggle feature discrimination display
                                showParamPass = ~showParamPass;

                            case 'p'
                                % change discrimination scheme (BWD only)
                                if strcmp(caller,'BWD')
                                    % keep current discriminator settings in
                                    % memory
                                    hfDiscrim.UserData.CurrentSettings = htDiscrim.Data;

                                    % change species discriminator
                                    hfDiscrim.Visible = 'on';
                                    uiwait(hfDiscrim);

                                    % reset if OK was pressed, otherwise do nothing
                                    if hfDiscrim.UserData.Proceed
                                        disp('Changing discriminators...')
                                        state = 'reset';
                                    end
                                else
                                    disp('Bad input')
                                end

                            case 'c'
                                % add comment
                                comment = xlsOutData.Comment(ii);
                                newComment = inputdlg('Enter comment','Comment',[1 50],comment);
                                if ~isempty(newComment)
                                    comment = newComment;
                                    xlsOutData.Comment(ii) = comment;
                                end

                            case 'u'
                                % update spreadsheet (manually save to disk)
                                %%% Ask user to create spreadsheet if it 
                                %%% doesn't exist yet
                                if isempty(xlsOutPath)
                                    [xlsOutPath,~] = saveTableAs(xlsOutData,dirPath_results,depName,targetName);
                                else
                                    saveTable(xlsOutData,xlsOutPath,false)
                                end

                            case 'U'
                                % save validation table as new spreadsheet
                                [xlsOutPathNew,newSaveCanceled] = saveTableAs(xlsOutData,dirPath_results,depName,targetName);
                                if ~newSaveCanceled
                                    xlsOutPath = xlsOutPathNew;
                                end

                            case 'e'
                                % change event

                                % update event switcher table and activate
                                % window
                                updateEventSwitchTable(htEvent,xlsOutData,ii)
                                drawnow
                                hfEvent.Visible = 'on';
                                uiwait(hfEvent);

                                % check if OK or Last button was pressed
                                if hfEvent.UserData.Proceed
                                    % check if new selection differs from
                                    % last one. If so, change ii and set
                                    % state to "done"
                                    iSelectedEvent = find(cell2mat(htEvent.Data(:,3)));
                                    if iSelectedEvent ~= ii
                                        ii = iSelectedEvent - 1;
                                        state = 'done';
                                    end
                                end

                            case 'q'
                                % quit
                                state = 'exit';

                            otherwise
                                % expect species ID
                                evalID = str2double(opt);
                                if isnan(evalID)
                                    % bad input
                                    disp('Bad input')
                                else
                                    % species ID - add to spreadsheet
                                    xlsOutData.Species(ii) = evalID;
                                    state = 'done';
                                end
                        end
                        
                    case 'done'
                        % clear all plots and mark end of processing for 
                        % current event
                        structfun(@clf,hfMain)
                        drawnow
                        doNextEvent = true;

                    case 'exit'
                        % break from event loop
                        doQuit = true;
                        break

                    otherwise
                        error('Unrecognized state')
                end
            end
            
            % end of processing for current file.
            % Check if termination was requested, and break if so.
            if doQuit
                break
            end
            
            % increment ii
            ii = ii + 1;
        end

        % save spreadsheet
        if isempty(xlsOutPath)
            processSave = true;
            while processSave
                [xlsOutPath,newSaveCanceled] = saveTableAs(xlsOutData,dirPath_results,depName,targetName);
                if newSaveCanceled
                    saveOpt = questdlg('Are you sure you want to quit without saving?','Warning: unsaved data','Yes','No','Yes');
                    if strcmp(saveOpt,'Yes')
                        processSave = false;
                    end
                else
                    processSave = false;
                end
            end
        else
            saveTable(xlsOutData,xlsOutPath,false)
        end
        
    % emergency saving for unexpected errors
    catch ME
        warning('Encountered error:\n %s\n\nAttempting to save backup table.',ME.getReport)
        dt = datetime();
        dt.Format = 'yyyyMMddHHmmss';
        if isempty(xlsOutPath)
            xlsOutPath = fullfile(dirPath_results,sprintf('%s_%s_Validated.xlsx',depName,targetName));
        end
        [xlsOutDir,xlsOutName,xlsOutExt] = fileparts(xlsOutPath);
        xlsOutEName = [xlsOutName,'_',char(dt)];
        xlsOutEPath = fullfile(xlsOutDir,[xlsOutEName,xlsOutExt]);
        saveTable(xlsOutData,xlsOutEPath,false)   

    end
    
    % delete all windows
    structfun(@delete,hfMain);
    %%% switcher windows MUST be deleted, close won't work
    delete(hfDiscrim)
    delete(hfEvent)
end


%% getSpreadsheets -------------------------------------------------------
function [xlsOutData,xlsOutPath,iStart,eventsTable,validInput] = getSpreadsheets(dirPath_results,depName,targetName,dtFormat,eventMergeOpt,eventMergeVal,caller)
% Prompts user for an input "Validated" or "RawEvents" spreadsheet, and
% loads the files as needed.
% If "RawEvents" sheet is entered, a new set of validations will be made.
% If "Validated" sheet is entered, the existing validations will be loaded,
% but a corresponding RawEvents sheet must still exist.
% If no appropriate spreadsheets are found or the user cancels, the program 
% will quit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initilize output
    xlsOutData = [];
    xlsOutPath = '';
    iStart = 0;
    eventsTable = [];
    
    % set base file names
    xlsName_events = sprintf('%s_%s_RawEvents.xlsx',depName,targetName);
    xlsNameBase_IDs = sprintf('%s_%s_Validated',depName,targetName);

    % get list of loadable spreadsheets
    [~, fNames] = TWD_Common.Utilities.listFiles(dirPath_results, 'xlsx', 'Recursive',false);
    isFile_events = contains(fNames,xlsName_events) & ~contains(fNames,'~$');
    isFile_IDs = contains(fNames,xlsNameBase_IDs) & ~contains(fNames,'~$');
    validNames = isFile_events | isFile_IDs;
    fNames = fNames(validNames);
    if isempty(fNames)
        warning('No valid input spreadsheets found for deployment %s, target %s',depName,targetName)
        validInput = false;
        return
    end
    
    % prompt user to select spreadsheet
    prompt = {...
        'Select input file.';...
        '';...
        'Choose the "RawEvents" spreadsheet to start a new set of';...
        'validations.';...
        '';...
        'Choose a "Validated" spreadsheet to continue or edit an';...
        'existing set (NOTE: the corresponding "RawEvents"';...
        'spreadsheet must still exist).';
        ''};
    [opt,validInput] = listdlg(...
        'ListString',fNames,...
        'PromptString',prompt,...
        'SelectionMode','single',...
        'ListSize',[300,200],...
        'Name','Input Selection');
    if ~validInput
        disp('No file selected')
        return
    end
    inFileName = fNames{opt};
    inFilePath = fullfile(dirPath_results,inFileName);
    
    % extract deployment name, target, and spreadsheet type from filename
    if contains(inFileName,'RawEvents')
        sheetType = 'RawEvents';
    else
        sheetType = 'Validated';
    end
    %nameParts = strsplit(inFileName,{'_','.'});
    %try
        %depName = nameParts{1};
        %targetName = nameParts{2};
    %    sheetType = nameParts{3};
    %catch
    %    warning('Filename format incorrect.\nCould not read deployment name, target, and/or spreadsheet type.')
    %    validInput = false;
    %    return
    %end
    
    % read input spreadsheet
    try
        inTable = readtable(inFilePath);
    catch
        warning('%s could not be opened',inFileName);
        validInput = false;
        return
    end
    
    % process input file and set output according to spreadsheet type
    switch sheetType
        case 'Validated'
            xlsOutData = inTable;
            xlsOutPath = inFilePath;
            iStart = find(isnan(xlsOutData.Species),1);
            
            % make sure comments are cell arrays of strings
            if isnumeric(xlsOutData.Comment)
                xlsOutData.Comment = num2cell(xlsOutData.Comment);
                for ii = 1:height(xlsOutData)
                    if isnan(xlsOutData.Comment{ii})
                        xlsOutData.Comment{ii} = '';
                    else
                        xlsOutData.Comment{ii} = num2str(xlsOutData.Comment{ii});
                    end
                end
            end
            
            % find and load events spreadsheet too
            xlsPath_events = fullfile(dirPath_results,xlsName_events);
            try
                eventsTable = readtable(xlsPath_events);
            catch
                warning('RawEvents spreadsheet not found or could not be opened');
                validInput = false;
                return
            end
            disp('Resuming validations')
            
        case 'RawEvents'
            eventsTable = inTable;
            
            % create new ID table (do not assign output path yet)
            xlsOutData = createIDSpreadsheet(eventsTable,dtFormat,eventMergeOpt,eventMergeVal,caller);
            xlsOutPath = '';
            iStart = 1;
            disp('Starting new validations')
            
        otherwise
            error('Encountered bug: spreadsheet type not recognized')
    end
    
    % all files read successfully - mark input as valid
    validInput = true;
end

% createIDSpreadsheet [SUBFUNCTION] .......................................
function xlsOutData = createIDSpreadsheet(eventsTableRaw,dtFormat,eventMergeOpt,eventMergeVal,caller)
% Creates a spreadsheet to be filled out with species presence info.
% If applicable, event merging is also performed at this stage (so the
% resulting spreadsheet may have fewer rows than the event sheet).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % merge events if necessary
    [eventsTableNew,targetTime] = mergeEvents(eventsTableRaw,dtFormat,eventMergeOpt,eventMergeVal,caller);
    
    % create output table
    nEvents = height(eventsTableNew);
    targetTimeCellstr = strrep(cellstr(targetTime),':','_');
    switch caller
        case 'BWD'
            %%% time
            timeAppendVars = {'TimeWithBeaked'};
            timeTable = table(targetTimeCellstr,'VariableNames',timeAppendVars);
            %%% ID
            IDAppendVars = {'NClicksBeaked_Validation','Species','Comment'};
            IDAppendTable = table(NaN(nEvents,1),NaN(nEvents,1),repmat({''},nEvents,1),'VariableNames',IDAppendVars);
        case 'SWD'
            %%% time
            timeAppendVars = {'TimeWithSperms'};
            timeTable = table(targetTimeCellstr,'VariableNames',timeAppendVars);
            %%% ID
            IDAppendVars = {'Species','Comment'};
            IDAppendTable = table(NaN(nEvents,1),repmat({''},nEvents,1),'VariableNames',IDAppendVars);
        otherwise
            error('Unknown caller "%s"',caller)
    end
    xlsOutData = [eventsTableNew(:,1:2),timeTable,eventsTableNew(:,3:end),IDAppendTable];
end

% mergeEvents [SUBFUNCTION] ...............................................
function [eventsTableNew,targetTime] = mergeEvents(eventsTableRaw,dtFormat,eventMergeOpt,eventMergeVal,caller)
% Merges adjacent events based on specified criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize
    eStartsRaw = datetime(eventsTableRaw.StartTime,'Format',dtFormat);
    eEndsRaw = datetime(eventsTableRaw.EndTime,'Format',dtFormat);
    eventDur = eEndsRaw - eStartsRaw;
    switch caller
        case 'BWD'
            nClicksRaw = [eventsTableRaw.NClicksAll,eventsTableRaw.NClicksBeaked_EventDet];
        case 'SWD'
            nClicksRaw = eventsTableRaw.NClicks;
    end
    nCol = size(nClicksRaw,2);

    switch eventMergeOpt
        case 'none' % no merging
            eStartsNew = eStartsRaw;
            eEndsNew = eEndsRaw;
            nClicksNew = nClicksRaw;
            targetTime = eventDur;
            
        case 'timegap' % merge events within a certain time limit of one another
            % get gaps between events
            gaps = eStartsRaw(2:end) - eEndsRaw(1:(end-1));
            
            % assess gap length, and merge events that are too close
            longGap = gaps >= eventMergeVal;
            eStartsNew = eStartsRaw([true;longGap]);
            eEndsNew = eEndsRaw([longGap;true]);
            nEventsNew = numel(eStartsNew);
            
            % now sum data for event segments that were merged
            originalStart = ismember(eStartsRaw,eStartsNew);
            iSum = cumsum(originalStart);
            nClicksNew = zeros(nEventsNew,nCol);
            targetTime = duration(zeros(nEventsNew,3));
            for ii = 1:nEventsNew
                doSumii = iSum == ii;
                nClicksNew(ii,:) = sum(nClicksRaw(doSumii,:),1);
                targetTime(ii) = sum(eventDur(doSumii));
            end
            
        case 'calendar' % merge all events that occur within a certain calendar period, e.g. every day or hour
            %%% For this to work, eventMergeVal must be a string which is
            %%% one of the supported "dateshift" units. Most relevant are:
            %%% 'hour', 'day', 'week', 'month', and 'year'
            
            % get adjusted start times for each event
            eStartsInterim = dateshift(eStartsRaw,'start',eventMergeVal,'current');
            
            % find redundant starts and merge them
            [eStartsNew,~,iRawNew] = unique(eStartsInterim);
            nEventsNew = numel(eStartsNew);
            nClicksNew = zeros(nEventsNew,nCol);
            targetTime = duration(zeros(nEventsNew,3));
            for ii = 1:nEventsNew
                rawii = iRawNew == ii;
                nClicksNew(ii,:) = sum(nClicksRaw(rawii,:),1);
                targetTime(ii) = sum(eEndsRaw(rawii) - eStartsRaw(rawii));
            end
            
            % get end times
            eEndsNew = dateshift(eStartsNew,'start',eventMergeVal,'next');
            
        case 'timebin' % merge events that occur within a designated time bin
            binSize = eventMergeVal;
            
            % determine unit from which to start binning 
            %%% This is the next-greatest unit after the greatest unit in 
            %%% the time bin (so for example, 10-min bins start on 1st 
            %%% hour, 10-hour bins start on 1st day, etc.). Note that weeks
            %%% and months are not valid units, to use those it would be
            %%% better to use the "calendarDuration" datatype which has its
            %%% own drawbacks.
            timeUnitVec = {'year','day','hour','minute'};
            timeFcns = {@days,@hours,@minutes};
            nUnits = numel(timeUnitVec);
            unitIdx = 1;
            needUnit = true;
            while needUnit && unitIdx<nUnits
                currentUnitVal = feval(timeFcns{unitIdx},binSize);
                if floor(currentUnitVal) > 0
                    % unit found
                    needUnit = false;
                else
                    % evaluate next unit
                    unitIdx = unitIdx + 1;
                end
            end
            startUnit = timeUnitVec{unitIdx};

            % create array of bins
            startBinning = dateshift(min(eStartsRaw),'start',startUnit,'current');
            binStarts = (startBinning:binSize:max(eEndsRaw))';
            binEnds = binStarts + binSize;
            %nBins = numel(binStarts);
            
            % determine which bins to keep, type 1:
            % find which bins each recording fits into
            %%% use this approach if/when segment splitting is enabled
            %{
            nEventsRaw = numel(eStartsRaw);
            binHasData = false(nBins,1);
            for ii = 1:nEventsRaw
                binHasDataii = binEnds > eStartsRaw(ii) & binStarts < eEndsRaw(ii);
                binHasData = binHasData | binHasDataii;
            end
            eStartsNew = binStarts(binHasData);
            eEndsNew = binEnds(binHasData);
            %}
            
            % determine which bins to keep, type 2:
            % assign bins to each segment and accumulate 
            %iBins = zeros(nEventsRaw,1);
            %for ii = 1:nEventsRaw
            %    iBinii = find(binEnds > eStartsRaw(ii), 1, 'start');    
            %    iBins(ii) = iBinii;
            %end
            %iBins = unique(iBins);
            %eStartsNew = binStarts(iBins);
            %eEndsNew = binEnds(iBins);
            %%%%%
            nEventsRaw = numel(eStartsRaw);
            eStartsInterim = eStartsRaw;
            for ii = 1:nEventsRaw
                iBinii = find(binEnds > eStartsRaw(ii), 1, 'first');    
                eStartsInterim(ii) = binStarts(iBinii);
            end
            
             % find redundant starts and merge them
            [eStartsNew,~,iRawNew] = unique(eStartsInterim);
            nEventsNew = numel(eStartsNew);
            nClicksNew = zeros(nEventsNew,nCol);
            targetTime = duration(zeros(nEventsNew,3));
            for ii = 1:nEventsNew
                rawii = iRawNew == ii;
                nClicksNew(ii,:) = sum(nClicksRaw(rawii,:),1);
                targetTime(ii) = sum(eEndsRaw(rawii) - eStartsRaw(rawii));
            end
            
            % get end times
            eEndsNew = eStartsNew + binSize;
            
        otherwise
            error('Unrecognized event merging option')
    end
    
    % convert nClicksNew to cell array, where each column is one cell
    %%% This will help deal with different possible numbers of columns
    %%% between BWD and SWD in a single line
    nClicksNewCell = mat2cell(nClicksNew',ones(1,nCol))';
    nClicksNewCell = cellfun(@transpose,nClicksNewCell,'UniformOutput',false);
    
    % redo table
    eventsTableNew = table(cellstr(eStartsNew),cellstr(eEndsNew),nClicksNewCell{:},...
        'VariableNames',eventsTableRaw.Properties.VariableNames);
end


%% createMainWindows ------------------------------------------------------
function hf = createMainWindows(caller,clickDiscrim,monID)
% Creates and positions the main figure windows relative to the main
% display resolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% Original positions for reference (not OuterPosition):
%   hfSummary = figure('Position',[10,415,870,580]);
%   hfOverlay = figure('Position',[10,25,1000,350]);
%   hfClicks = figure('Position',[895,595,600,400]);
%   hfClickParams = figure('Position',[1510,595,400,400]);
% Best approach may be to have "screen" dimensions and "saving" dimensions
    
    % define left-right-bottom figure OuterPosition buffer, and Windows 10 
    % toolbar offset (note that these are invariant to screen resolution)
    buff = 7;
    yShift = 40;

    % get normalized window dimensions
    [nposSummary,nposOverlay,nposClick,nposClickSeq,nposClickParams] = defineWinDim(yShift);
    
    % get actual screen resolution
    mp = get(0,'MonitorPositions');
    try
        scrRes = mp(monID,:);
    catch
        warning('Cannot find monitor %d, using monitor 1 instead',monID)
        scrRes = mp(1,:);
    end
    scrXOffset = scrRes(1) - 1;
    scrYOffset = scrRes(2) - 1;
    scrWidth = scrRes(3);
    scrHeightActual = scrRes(4);
    scrHeightEffective = scrHeightActual - yShift;
    
    % Transformations need to be done as follows:
    %   xf = 1 + xfNorm.*w - buff
    %   yf = 1 + yfNorm.*hE - buff + yShift
    %   wf = wfNorm.*w + 2*buff
    %   hf = hfNorm.*hA + buff
    % create translation and scale vectors
    %addVec = [1-buff, 1+yShift-buff, 2*buff, buff];
    addVec = [1-buff+scrXOffset, 1+yShift-buff+scrYOffset, 2*buff, buff];
    multVec = [scrWidth, scrHeightEffective, scrWidth, scrHeightEffective];
    
    % get figure positions in pixels
    pos = struct();
    pos.Summary = round(nposSummary.*multVec + addVec);
    pos.Overlay = round(nposOverlay.*multVec + addVec);
    pos.Click = round(nposClick.*multVec + addVec);
    pos.ClickSeq = round(nposClickSeq.*multVec + addVec);
    pos.ClickParams = round(nposClickParams.*multVec + addVec);
    
    % create figures according to detector type
    figNamesAll = fieldnames(pos);
    switch caller
        case 'BWD'
            figNames = figNamesAll;
        case 'SWD'
            figNames = figNamesAll(~ismember(figNamesAll,'Overlay'));
        otherwise
            error('Unknown caller "%s", caller')
    end
    nFigs = numel(figNames);
    hf = struct();
    for ii = 1:nFigs
        figii = figNames{ii};
        hf.(figii) = figure('OuterPosition',pos.(figii));
    end
    
    % create and add data to Click Params figure
    % (stored within UserData and may also apply figure resize)
    createClickParamsFigData(hf.ClickParams,clickDiscrim);
end

% defineWinDim [SUBFUNCTION] ..............................................
function [npos_Summary,npos_Overlay,npos_Click,npos_ClickSeq,npos_ClickParams] = defineWinDim(yShift)
% Defines normalize OuterPosition values for each window, based on
% available screen space and aspect ratios of the original SBP figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define reference effective screen size and aspect ratio (x/y)
    scrResRef = [1920, 1080-yShift+14];
    scrAR = scrResRef(1)/scrResRef(2);

    % define 1-to-1 OuterPosition aspect ratios (x/y) for Summary, Overlay, 
    % and Click windows (based on SBP)
    oprSummary = (886/673)/scrAR;
    oprOverlay = (1016/443)/scrAR;
    oprClick = (616/493)/scrAR;
    
    % derive normalized OuterPosition values for each window
    %%% Summary
    nposH_Summary = 0.601;
    nposW_Summary = nposH_Summary*oprSummary;
    nposY_Summary = 1-nposH_Summary;
    nposX_Summary = 0;
    %%% Overlay
    nposH_Overlay = 1-nposH_Summary;
    nposW_Overlay = nposH_Overlay*oprOverlay;
    nposY_Overlay = 0;
    nposX_Overlay = 0;
    %%% Click
    nposH_Click = 0.474;
    nposW_Click = nposH_Click*oprClick;
    nposY_Click = 1-nposH_Click;
    nposX_Click = nposW_Summary;
    %%% ClickSequence
    nposH_ClickSeq = nposH_Click;
    nposW_ClickSeq = 1 - (nposX_Click + nposW_Click);
    nposY_ClickSeq = nposY_Click;
    nposX_ClickSeq = nposX_Click + nposW_Click;
    %%% ClickParams
    nposH_ClickParams = 1 - nposH_Click; %0.340;
    nposW_ClickParams = 1 - nposW_Overlay;
    nposY_ClickParams = 0; %1 - (nposH_Click + nposH_ClickParams);
    nposX_ClickParams = nposW_Overlay;
    
    % create output position vectors
    npos_Summary = [nposX_Summary, nposY_Summary, nposW_Summary, nposH_Summary];
    npos_Overlay = [nposX_Overlay, nposY_Overlay, nposW_Overlay, nposH_Overlay];
    npos_Click = [nposX_Click, nposY_Click, nposW_Click, nposH_Click];
    npos_ClickSeq = [nposX_ClickSeq, nposY_ClickSeq, nposW_ClickSeq, nposH_ClickSeq];
    npos_ClickParams = [nposX_ClickParams, nposY_ClickParams, nposW_ClickParams, nposH_ClickParams];
end

% createClickParamsFigData [SUBFUNCTION] ..................................
function createClickParamsFigData(hfClickParams,clickDiscrim)
% Defines various data useful for displaying click parameter and
% discrimination info, and stores within figure UserData. May also resize
% the figure to optimize dimensions based on amount of content.
% Assumes 3 tables positioned as follows:
%   - Occurrence (top left)
%   - Discrimination (bottom left)
%   - Parameters (right)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define colours for UITable cells
    colActiveHex = struct(...
        'pass','#bfffbf',... % Hex for [0.75,1,0.75]; pale green
        'fail','#ffbfbf',... % Hex for [1,0.75,0.75]; pale red 
        'unassessed','#ffffbf'); % Hex for [1,1,0.75]; pale yellow
    colInactiveHex = struct(...
        'pass','#80bf80',... % Hex for [0.5,0.75,0.5]; dark green
        'fail','#bf8080',... % Hex for [0.75,0.5,0.5]; dark red 
        'unassessed','#bfbf80'); % Hex for [0.75,0.75,0.5]; dark yellow
    colActiveDec = struct(...
        'pass',[0.75,1,0.75],... 
        'fail',[1,0.75,0.75],... 
        'unassessed',[1,1,0.75]); 
    colInactiveDec = struct(...
        'pass',[0.5,0.75,0.5],... 
        'fail',[0.75,0.5,0.5],... 
        'unassessed',[0.75,0.75,0.5]);
    
    % get number of discriminators
    nDiscrim = numel(clickDiscrim);
    
    % get number of features
    %%% ClickDiscriminator property "CritTable" lists all available
    %%% features, but Params table should show only those used. So, need
    %%% to loop through each discriminator and accumulate those features
    %%% that are used. (can just be a sum of UseCategory > 0)
    nFeaturesTotal = height(clickDiscrim(1).critTable);
    featureUseList = false(nFeaturesTotal,1);
    for ii = 1:nDiscrim
        featureUseListii = clickDiscrim(ii).critTable.UseCategory > 0;
        featureUseList = featureUseList | featureUseListii;
    end
    nFeatures = sum(featureUseList);
    
    % DETERMINE UI COMPONENT POSITIONS
    
    % define static dimensions
    %%% padding
    borderPad = 20;
    intPad = 25;
    txtYPadCorrect = 10;
    %%% UITable parts
    tableHeadHeight = 20;
    tableRowHeight = 18; % for adaptive table heights, always add +2 to cumulative height (including headers and horz scrollbar)
    tableVertScrollWidth = 18;
    tableHorzScrollHeight = 17;
    %%% text objects
    txtHeight = 20;
    txtFontSize = 10;
    
    % get current figure dimensions
    figPosOriginal = hfClickParams.Position;
    figOutPosOriginal = hfClickParams.OuterPosition;
    figPos = figPosOriginal;
    
    % define table dimensions
    %%% Occurance table
    occurTableColWidths = {50, 80};
    occurTableWidth = sum([occurTableColWidths{:},2]);
    occurTableHeight = tableRowHeight*2 + 2;
    %%% Discrimination table
    discrimTableColWidths = {50, 80-tableVertScrollWidth};
    discrimTableWidth = sum([discrimTableColWidths{:},tableVertScrollWidth,2]);
    discrimTableHeightMin = round(tableRowHeight*1.5) + 2; % limiting factor for figure height (will resize figure if too small)
    discrimTableHeightBest = tableRowHeight*nDiscrim + tableHorzScrollHeight + 2; % used for shrinking figure if possible
    %%% Parameters table
    paramTableColWidths = [{100, 70},repmat({50},1,nDiscrim)];
    paramTableWidthBest = sum([paramTableColWidths{:},tableVertScrollWidth,2]);
    paramTableWidthMin = 100; % limiting factor for figure width (will resize figure if too small)
    paramTableHeightBest = tableHeadHeight + tableRowHeight*nFeatures + tableHorzScrollHeight + 2; % used for shrinking figure if possible
    
    % define anonymous functions for determining X or Y positions 
    % based on other objects
    xRelativeTo = @(pos) sum(pos([1,3]));
    yRelativeTo = @(pos) sum(pos([2,4]));
    
    % define UI component positions (top-down, i.e. assuming Y is measured
    % from top to bottom - this will be corrected later).
    %%% Occurance
    pos_txtOccur = [borderPad, borderPad-txtYPadCorrect, occurTableWidth, txtHeight];
    pos_tableOccur = [borderPad, yRelativeTo(pos_txtOccur), occurTableWidth, occurTableHeight];
    %%% Discrimination
    pos_txtDiscrim = [borderPad, yRelativeTo(pos_tableOccur)+intPad-txtYPadCorrect, discrimTableWidth, txtHeight];
    pos_tableDiscrim = [borderPad, yRelativeTo(pos_txtDiscrim), discrimTableWidth, hfClickParams.Position(4)-borderPad-yRelativeTo(pos_txtDiscrim)];
    %%% Parameters
    pos_txtParam = [xRelativeTo(pos_tableOccur)+intPad, borderPad-txtYPadCorrect, hfClickParams.Position(3)-borderPad-(xRelativeTo(pos_tableOccur)+intPad), txtHeight];
    pos_tableParam = [pos_txtParam(1), yRelativeTo(pos_txtParam), pos_txtParam(3), hfClickParams.Position(4)-borderPad-yRelativeTo(pos_txtParam)]; 
    
    % adjust figure dimensions
    
    % check if figure too short
    figTooShort = pos_tableDiscrim(4) < discrimTableHeightMin;
    if figTooShort
        % too short; increase figure height
        %warning('ClickParams figure height increased')
        hfClickParams.Position(4) = pos_tableDiscrim(2) + discrimTableHeightMin + borderPad;
        drawnow
        pos_tableDiscrim(4) = discrimTableHeightMin;
        pos_tableParam(4) = hfClickParams.Position(4) - borderPad - pos_tableParam(2);
    else
        % check if figure too tall
        figHeightBest1 = pos_tableDiscrim(2) + discrimTableHeightBest + borderPad;
        figHeightBest2 = pos_tableParam(2) + paramTableHeightBest + borderPad;
        figHeightBest = max([figHeightBest1,figHeightBest2]);
        figTooTall = figPos(4) > figHeightBest;
        if figTooTall
            % too tall; reduce figure height
            %warning('ClickParams figure height reduced')
            hfClickParams.Position(4) = figHeightBest;
            drawnow
            pos_tableDiscrim(4) = hfClickParams.Position(4) - borderPad - pos_tableDiscrim(2);
            pos_tableParam(4) = hfClickParams.Position(4) - borderPad - pos_tableParam(2);
            figOutHeightDiff = figOutPosOriginal(4) - hfClickParams.OuterPosition(4);
            hfClickParams.OuterPosition(2) = hfClickParams.OuterPosition(2) + figOutHeightDiff;
            drawnow
        end
    end
    
    % check if figure too narrow
    figTooNarrow = pos_tableParam(3) < paramTableWidthMin;
    if figTooNarrow
        % too narrow; increase figure width
        %warning('ClickParams figure width increased')
        hfClickParams.Position(3) = pos_tableParam(1) + paramTableWidthMin + borderPad;
        drawnow
        pos_tableParam(3) = paramTableWidthMin;
        pos_txtParam(3) = paramTableWidthMin;
    else
        % check if figure too long
        figWidthBest = pos_tableParam(1) + paramTableWidthBest + borderPad;
        if hfClickParams.Position(3) > figWidthBest
            % figure too long; reduce figure width
            %warning('ClickParams figure width reduced')
            hfClickParams.Position(3) = figWidthBest;
            drawnow
            pos_tableParam(3) = hfClickParams.Position(3) - borderPad - pos_tableParam(1);
            pos_txtParam(3) = pos_tableParam(3);
            figOutWidthDiff = figOutPosOriginal(3) - hfClickParams.OuterPosition(3);
            hfClickParams.OuterPosition(1) = hfClickParams.OuterPosition(1) + figOutWidthDiff;
            drawnow
        end
    end
    
    % flip top-down UI positions to bottom-up
    top2bottom = @(pos) [pos(1), hfClickParams.Position(4)-(pos(2)+pos(4)), pos(3), pos(4)];
    pos_txtOccur = top2bottom(pos_txtOccur);
    pos_tableOccur = top2bottom(pos_tableOccur);
    pos_txtDiscrim = top2bottom(pos_txtDiscrim);
    pos_tableDiscrim = top2bottom(pos_tableDiscrim);
    pos_txtParam = top2bottom(pos_txtParam);
    pos_tableParam = top2bottom(pos_tableParam);
    
    % save data to figure UserData
    hfClickParams.UserData = struct();
    hfClickParams.UserData.colActiveHex = colActiveHex;
    hfClickParams.UserData.colInactiveHex = colInactiveHex;
    hfClickParams.UserData.colActiveDec = colActiveDec;
    hfClickParams.UserData.colInactiveDec = colInactiveDec;
    hfClickParams.UserData.txtFontSize = txtFontSize;
    hfClickParams.UserData.pos_txtOccur = pos_txtOccur;
    hfClickParams.UserData.pos_tableOccur = pos_tableOccur;
    hfClickParams.UserData.pos_txtDiscrim = pos_txtDiscrim;
    hfClickParams.UserData.pos_tableDiscrim = pos_tableDiscrim;
    hfClickParams.UserData.pos_txtParam = pos_txtParam;
    hfClickParams.UserData.pos_tableParam = pos_tableParam;
    hfClickParams.UserData.occurTableColWidths = occurTableColWidths;
    hfClickParams.UserData.discrimTableColWidths = discrimTableColWidths;
    hfClickParams.UserData.paramTableColWidths = paramTableColWidths;
end


%% inspectAllEventClicks --------------------------------------------------
function [clickStarts,clickMatIdx,clickPassedDiscrim,clickParamVals,clickParamPass,clickQ] = inspectAllEventClicks(eStart,eEnd,eSegStarts,eSegEnds,matStarts,matFileNames,dirPath_mat,clickDiscrim)
% Loads all clicks in an event, returns their start times, MAT file 
% indices, and a measure of quality, and also runs all discriminators 
% against them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES: 
% - click quality is currently measured as SNR, but might be useful to try 
%   PP also.

    % initialize
    nDiscrim = numel(clickDiscrim);
    clickStarts = datetime.empty(0,1);
    clickMatIdx = double.empty(0,1);
    clickPassedDiscrim = logical.empty(0,nDiscrim);
    clickParamVals = table();
    clickParamPass = cell(nDiscrim,1);
    clickQ = double.empty(0,1);
    clickFields = {...
        'pos';...
        'peakFr';...
        'F0';...
        'dur';...
        'slope';...
        'nSamples';...
        'bw3db';...
        'bw10db';...
        'rmsSignal';...
        'rmsNoise';...
        'snr';...
        'ppSignal';...
        'yFilt';...
        'yNFilt';...
        'specClick';...
        'specNoise'};
    nClickFields = numel(clickFields);

    % get pseudo-end times for MAT files
    %%% Corresponds to onset of next MAT file.
    %%% Much faster than loading each MAT to get duration.
    [recStarts,~,iMatRec] = unique(matStarts);
    recEndsPseudo = [recStarts(2:end);datetime(Inf,Inf,Inf)];
    matEndsPseudo = recEndsPseudo(iMatRec);
    
    % get event segments within event range
    %%% Small detail: only segments with a start time within the wanted
    %%% event range are considered part of that event; end times don't
    %%% count.
    isMemberSegment = (eSegStarts >= eStart) & (eSegStarts < eEnd);
    eSegStartsSub = eSegStarts(isMemberSegment);
    eSegEndsSub = eSegEnds(isMemberSegment);
    
    % loop through each MAT file, load the ones that contain segments of 
    % interest, and inspect the relevant clicks
    nMat = numel(matStarts);
    for ii = 1:nMat
        % look for segments in range
        matStartii = matStarts(ii);
        matEndPseudoii = matEndsPseudo(ii);
        iSegsii = find((eSegStartsSub >= matStartii) & (eSegEndsSub <= matEndPseudoii));
        nSegmentsii = numel(iSegsii);
        
        % if there are segments, load MAT file
        if nSegmentsii > 0
            matPathii = fullfile(dirPath_mat,matFileNames{ii});
            dataii = load(matPathii);
            clickStartsii = matStartii + seconds(dataii.pos(:,1));
            clickQii = dataii.snr;
            nClicksii = numel(clickStartsii);
            
            % if the end of the last segment of interest corresponds to the 
            % end of the recording, it's possible that there may be a few 
            % clicks of interest that will be left out, because segment end 
            % is a rounded estimate of recording end. If that's the case, 
            % then edit the end of the last segment so that it corresponds 
            % to the true end.
            recEndii = matStartii + seconds(sum(dataii.rawDur));
            iSegLastii = iSegsii(end);
            if dateshift(recEndii,'start','second') == eSegEndsSub(iSegLastii)
                eSegEndsSub(iSegLastii) = recEndii;
            end
            
            % determine which clicks are within the segments of interest
            isClickInRangeii = false(nClicksii,1);
            for jj = 1:nSegmentsii
                iSegjj = iSegsii(jj);
                eSegStartjj = eSegStartsSub(iSegjj);
                eSegEndjj = eSegEndsSub(iSegjj);
                isClickInRangejj = (clickStartsii >= eSegStartjj) & (clickStartsii < eSegEndjj);
                
                isClickInRangeii = isClickInRangeii | isClickInRangejj;
            end
            
            nClicksInRangeii = sum(isClickInRangeii);
            if nClicksInRangeii > 0
                % retain data for the clicks in range only
                for jj = 1:nClickFields
                    fieldjj = clickFields{jj};
                    dataii.(fieldjj) = dataii.(fieldjj)(isClickInRangeii,:);
                end
                clickStartsii = clickStartsii(isClickInRangeii);
                clickQii = clickQii(isClickInRangeii);
                
                % run all discriminators against all clicks
                clickPassedDiscrimii = false(nClicksInRangeii,nDiscrim);
                clickParamValsii = table();
                for jj = 1:nDiscrim
                    discrimjj = clickDiscrim(jj);
                    [passedjj,clickParamPassjj,clickParamValsjj] = discrimjj.evalTargetClicks(dataii);
                    clickPassedDiscrimii(:,jj) = passedjj;
                    
                    % add any new variables to the click parameter table
                    clickParamsii = clickParamValsii.Properties.VariableNames;
                    clickParamsjj = clickParamValsjj.Properties.VariableNames;
                    isNewParamjj = ~ismember(clickParamsjj,clickParamsii);
                    clickParamValsii = [clickParamValsii,clickParamValsjj(:,isNewParamjj)];
                    
                    % save pass status for each feature
                    clickParamPass{jj} = [clickParamPass{jj};clickParamPassjj];
                end
                
                % store info for each click
                %%% ignore warnings
                clickStarts = [clickStarts;clickStartsii];
                clickMatIdx = [clickMatIdx;repmat(ii,nClicksInRangeii,1)];
                clickPassedDiscrim = [clickPassedDiscrim;clickPassedDiscrimii];
                clickParamVals = [clickParamVals;clickParamValsii];
                clickQ = [clickQ;clickQii];
            end
        end
    end
end


%% limitNumClicks ---------------------------------------------------------
function subset = limitNumClicks(clickStarts,eSegStarts,eSegEnds,isTarget,clickQ,nClicksMax,caller)
% Returns logical indices of which clicks to keep, based on the maximum
% allowed number of clicks. Clicks to keep are prioritized based on the
% quality of their segments, where quality is measured in one of multiple
% ways. For BWD, segment quality is the proportion of target clicks in the
% segment. For SWD (where all clicks are target clicks), segment quality is 
% the average quality of clicks within the segment, where click quality can
% be any single-number metric where higher is better (e.g. SNR or
% peak-to-peak amplitude). Segments are added until the cumulative 
% total number of clicks exceeds the allowed limit. After that, excess 
% clicks in the last segment are culled to make sure that the limit is 
% respected.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize
    nClicksAll = numel(clickStarts);
    nClicksTarget = sum(isTarget);

    % check if there are too many target clicks
    if nClicksTarget > nClicksMax
        
        % sort segments as appropriate
        switch caller
            case 'BWD'
                % sort by proportion of target clicks
                limitTypeStr = 'proportion of target clicks';
                [iSegSort,clickSegIdx,nClicksTargetSeg] = sortSeg_Prop(clickStarts,eSegStarts,eSegEnds,isTarget);
            case 'SWD'
                % sort by click "quality"
                limitTypeStr = 'mean quality of target clicks';
                [iSegSort,clickSegIdx,nClicksTargetSeg] = sortSeg_clickQ(clickStarts,eSegStarts,eSegEnds,isTarget,clickQ);
        end
         warning('%d clicks passed, limiting to %d;\nsegments prioritized based on %s',nClicksTarget,nClicksMax, limitTypeStr)
        
        % get number of target clicks in each segment (sorted)
        nClicksTargetSegSort = nClicksTargetSeg(iSegSort);
        
        % add segments in order of greatest proportion until maximum click
        % limit is reached
        segTargetClickSum = cumsum(nClicksTargetSegSort);
        iSegKeepSort = 1:find(segTargetClickSum > nClicksMax,1,'first');
        
        % keep track of lowest-quality segment retained, for which excess
        % clicks will be culled
        iSegCullSort = iSegKeepSort(end);
        
        % get number of excess clicks
        nClicksTargetExcess = segTargetClickSum(iSegCullSort) - nClicksMax;
        
        % revert segments to original order
        iSegKeep = iSegSort(iSegKeepSort);
        iSegCull = iSegSort(iSegCullSort);
        
        % cull clicks in the lowest-quality segment to keep number of
        % clicks under maximum limit
        cullClick = false(nClicksAll,1);
        if nClicksTargetExcess > 0
            isTargetClickInTruncSeg = ismember(clickSegIdx,iSegCull) & isTarget;
            iTargetClickInTruncSeg = find(isTargetClickInTruncSeg);
            iCull = iTargetClickInTruncSeg((end-nClicksTargetExcess+1):end);
            cullClick(iCull) = true;
        end
        
        % get subset of clicks
        isClickInGoodSeg = ismember(clickSegIdx,iSegKeep);
        subset = isTarget & isClickInGoodSeg & ~cullClick;
    else 
        subset = isTarget;
    end
end

% sortSeg_Prop [SUBFUNCTION] ..............................................
function [iSegSort,clickSegIdx,nClicksTargetSeg] = sortSeg_Prop(clickStarts,eSegStarts,eSegEnds,isTarget)
% Sort segments by proportion of target clicks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize
    nClicksAll = numel(clickStarts); %%% total number of clicks in event
    nSeg = numel(eSegStarts); %%% number of segments in event
    clickSegIdx = zeros(nClicksAll,1); %%% container to hold segment membership of each click
    nClicksTargetSeg = zeros(nSeg,1); %%% container to hold number of target clicks in each segment
    segPropTarget = zeros(nSeg,1); %%% container to hold proportion of target clicks in each segment
    
    % get number and proportion of target clicks in each segment
    for ii = 1:nSeg
        segStartii = eSegStarts(ii);
        segEndii = eSegEnds(ii);

        % isolate clicks in segment
        isClickInSegii = (clickStarts >= segStartii) & (clickStarts < segEndii);
        nClicksSegii = sum(isClickInSegii);
        clickSegIdx(isClickInSegii) = ii;

        % get proportion of target clicks
        isTargetii = isTarget & isClickInSegii;
        nClicksTargetii = sum(isTargetii);
        segPropTargetii = nClicksTargetii/nClicksSegii;

        % store results
        nClicksTargetSeg(ii) = nClicksTargetii;
        segPropTarget(ii) = segPropTargetii;
    end

    % change NaNs to -Inf in proportion vector, because NaNs cause
    % problems when sorting. NaN means that a segment has no clicks 
    % and therefore is useless, but they are treated as high values by 
    % the "sort" function. -Inf will give these segments lowest 
    % priority instead.
    segPropTarget(isnan(segPropTarget)) = -Inf;

    % sort segments by proportion of target clicks
    [~,iSegSort] = sort(segPropTarget,'descend');
end

% sortSeg_clickQ [SUBFUNCTION] ...............................................
function [iSegSort,clickSegIdx,nClicksTargetSeg] = sortSeg_clickQ(clickStarts,eSegStarts,eSegEnds,isTarget,clickQ)
% sort segments by mean "quality" of target clicks (variable clickQ).
% clickQ may be any single-number metric where higher is better, such as 
% SNR or peak-to-peak amplitude.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize
    nClicksAll = numel(clickStarts); %%% total number of clicks in event
    nSeg = numel(eSegStarts); %%% number of segments in event
    clickSegIdx = zeros(nClicksAll,1); %%% container to hold segment membership of each click
    nClicksTargetSeg = zeros(nSeg,1); %%% container to hold number of target clicks in each segment
    segMeanQ = zeros(nSeg,1); %%% container to hold mean quality of clicks in each segment
    
    % get number and average quality of clicks in each segment
    for ii = 1:nSeg
        segStartii = eSegStarts(ii);
        segEndii = eSegEnds(ii);

        % isolate clicks in segment
        isClickInSegii = (clickStarts >= segStartii) & (clickStarts < segEndii);
        clickSegIdx(isClickInSegii) = ii;

        % get average quality of target clicks in segment
        isTargetii = isTarget & isClickInSegii;
        nClicksTargetii = sum(isTargetii);
        clickQii = clickQ(isTargetii);
        segMeanQii = mean(clickQii);

        % store results
        nClicksTargetSeg(ii) = nClicksTargetii;
        segMeanQ(ii) = segMeanQii;
    end
    
    % change NaNs to -Inf in proportion vector, because NaNs cause
    % problems when sorting. NaN means that a segment has no clicks 
    % and therefore is useless, but they are treated as high values by 
    % the "sort" function. -Inf will give these segments lowest 
    % priority instead.
    segMeanQ(isnan(segMeanQ)) = -Inf;
    
    % sort segments by quality target clicks
    [~,iSegSort] = sort(segMeanQ,'descend');
end


%% getClicks --------------------------------------------------------------
function data = getClicks(clickStartTarget,clickMatIdxTarget,matStarts,matNames,dirPath_mat)
% Retrieves data struct containing clicks within the time range specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - watch out for hard-coded indices of click vs. non-click variables.
% Might be best to have var names extracted ahead of time and passed
% here...

    % initialize output
    %%% removed rawStart and rawDur. Field indices have been adjusted.
    data = struct(...
        'fs',[],...
        'offset',[],...
        'pos',[],...
        'peakFr',[],...
        'F0',[],...
        'dur',[],...
        'slope',[],...
        'nSamples',[],...
        'bw3db',[],...
        'bw10db',[],...
        'rmsSignal',[],...
        'rmsNoise',[],...
        'snr',[],...
        'ppSignal',[],...
        'yFilt',[],...
        'yNFilt',[],...
        'specClick',[],...
        'specNoise',[]);
    dataFields = fieldnames(data);
    nonClickFields = dataFields(1:2);
    clickFields = dataFields(3:end);
    nNonClickFields = numel(nonClickFields);
    nClickFields = numel(clickFields);  
    
    % isolate MAT files of interest
    iMatKeep = unique(clickMatIdxTarget);
    matNamesSub = matNames(iMatKeep);
    matStartsSub = matStarts(iMatKeep);
    nMatSub = numel(matNamesSub);
    
    % loop through each MAT file and load click data
    for ii = 1:nMatSub
        matNameii = matNamesSub{ii};
        matPathii = fullfile(dirPath_mat,matNameii);
        matStartii = matStartsSub(ii);
        dataii = load(matPathii);
        
        % get start times of clicks in current MAT
        clickStartsii = matStartii + seconds(dataii.pos(:,1));
        
        % check which clicks are in range
        isClickInRangeii = ismember(clickStartsii,clickStartTarget);
        
        % isolate data of interest
        for jj = 1:nClickFields
            fieldjj = clickFields{jj};
            dataVarjj = dataii.(fieldjj)(isClickInRangeii,:);
            data.(fieldjj) = [data.(fieldjj);dataVarjj];
        end
    end
    
    % add non-click info to struct
    if nMatSub > 0
        dataSample = dataii;
    else
        % if there were no files, then load the 1st one to get data. Also
        % use this to get proper dimensions for click variables.
        dataSample = load(fullfile(dirPath_mat,matNames{1}));
        for ii = 1:nClickFields
            fieldii = clickFields{ii};
            dataVarii = dataSample.(fieldii)([],:);
            data.(fieldii) = dataVarii;
        end
    end
    for ii = 1:nNonClickFields
        fieldii = nonClickFields{ii};
        data.(fieldii) = dataSample.(fieldii);
    end
end


%% plotSummaryFigs --------------------------------------------------------
function plotSummaryFigs(hf,data,depName,tStart,fRange,caller)
% Adapted from "bw_combined_analysis_verify"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % isolate data of interest
    ICITarget = diff(data.pos(:,1))*1000;
    FPeakTarget = data.peakFr;
    bw10dbTarget = data.bw10db;
    bw3dbTarget = data.bw3db;
    specClickTarget = data.specClick;
    specNoiseTarget = data.specNoise;
    nTarget = numel(FPeakTarget);
    Fs = data.fs;
    
    %SBP: delete ici not relevant to real ici
    ICITargetAdj = ICITarget(ICITarget>50 & ICITarget<1000);
    
    %SBP: find ici histogram bin with highest count (most likely to be true ici)
    ICIBins = 0:10:1000;
    [ICICounts,ICIEdges] = histcounts(ICITargetAdj,ICIBins);
    [~,iMaxICI] = max(ICICounts); % N used to be N(:), not sure why. Should be a vector regardless.
    ICIModeBin = [ICIEdges(iMaxICI),ICIEdges(iMaxICI+1)];
    
    %SBP: find peak frequency histogram bin with highest count
    %FPeakBins = 0:1:100;
    FPeakBins = floor(fRange(1)):1:ceil(fRange(2));
    [FPeakCounts,FPeakEdges] = histcounts(FPeakTarget,FPeakBins);
    [~,iMaxFPeak] = max(FPeakCounts);
    FPeakModeBin = [FPeakEdges(iMaxFPeak),FPeakEdges(iMaxFPeak+1)];
    
    %SBP: find median peak frequency, median -10dB and -3dB endpoints (for plots)
    peakMed = median(FPeakTarget);
    med10db = median(bw10dbTarget,1);
    med3db = median(bw3dbTarget,1);
    
    % plot graphs
    %SBP: sort spectras by peak frequency and prepare for plotting concetanated spectrogram
    [~,iSortFPeak] = sort(FPeakTarget);
    specSorted = specClickTarget(iSortFPeak,:);
    specSorted = specSorted.'; % Not sure why this is, watch out for it
    nfft = size(specSorted,1)*2;
    f = 0:(Fs/2000)/(nfft/2-1):Fs/2000; % This calculation is rather confusing, consider simplifying
    
    % calculate mean spectra for clicks and noise
    meanSpecClickTarget = TWD_Common.computeMeanSpectrum(specClickTarget);
    meanSpecNoiseTarget = TWD_Common.computeMeanSpectrum(specNoiseTarget);
    
    % set current figure and edit 
    set(0,'CurrentFigure',hf);
    clf
    hf.Name = sprintf('%s_%s',depName, char(tStart));
    
    % subplot 1: peak frequencies
    ha1 = subplot(2,2,1);
    histogram(ha1,FPeakTarget,FPeakBins,'EdgeColor','none','FaceColor','k','FaceAlpha',1)
    %xlim(ha1,[0,125])
    xlim(ha1,fRange)
    xlabel(ha1,'Peak Frequency (kHz)')
    ylabel(ha1,'Count')
    text(ha1,0.6,0.9,sprintf('F_{Peak} = %d-%d kHz',FPeakModeBin(1),FPeakModeBin(2)),'Unit','normalized')
    ha1.Tag = 'PeakFreqDistAxes';
    
    % subplot 2: ICI
    ha2 = subplot(2,2,2);
    histogram(ha2,ICITargetAdj,ICIBins,'EdgeColor','none','FaceColor','k','FaceAlpha',1)
    xlim(ha2,[0,1000])
    xlabel(ha2,'Inter-click interval (ms)')
    ylabel(ha2,'Count')
    text(ha2,0.6,0.9,sprintf('ICI = %d-%d ms',ICIModeBin(1),ICIModeBin(2)),'Unit','Normalized');
    ha2.Tag = 'ICIDistAxes';
    
    % subplot 3: power spectrum
    ha3 = subplot(2,2,3);
    ha3.NextPlot = 'add';
    plot(ha3,f,meanSpecClickTarget,'LineWidth',2)
    plot(ha3,f,meanSpecNoiseTarget,':k','LineWidth',2) 
    % ylim([-100 -20])
    %xlim(ha3,[5,95])
    switch caller
        case 'BWD'
            xlim(ha3,[5,95])
        case 'SWD'
            xlim(ha3,[2,95])
    end
    a3YLim = ha3.YLim;
    xlim(ha3,fRange)
    ylim(ha3,a3YLim)
    line(ha3, [med10db(1) med10db(1)],a3YLim, 'Color',[0.65 0.65 0.65], 'LineWidth',1); 
    line(ha3, [med10db(2) med10db(2)],a3YLim, 'Color',[0.65 0.65 0.65], 'LineWidth',1); 
    line(ha3, [med3db(1) med3db(1)],a3YLim, 'Color',[0.8 0.8 0.8], 'LineWidth',1); 
    line(ha3, [med3db(2) med3db(2)],a3YLim, 'Color',[0.8 0.8 0.8], 'LineWidth',1); 
    line(ha3, [peakMed peakMed],a3YLim, 'Color',[1 0.2 0.2], 'LineWidth',1);
    xlabel(ha3,'Frequency (kHz)') 
    ylabel(ha3,'Normalized amplitude (dB)') % what makes it normalized?
    title(ha3,sprintf('Mean click spectra, n=%d',nTarget),'FontWeight','bold')
    text(ha3,0.5,0.9,sprintf('-10dB lower = %d kHz',round(med10db(1))),'Unit','Normalized');
    ha3.Tag = 'PowSpecAxes';
    
    % subplot 4: concatenated spectrogram
    ha4 = subplot(2,2,4);
    conspec = imagesc(ha4,1:nTarget, f, specSorted); axis xy; colormap(gray);
    conspec.Tag = 'ConcatenatedSpectrogram';
    %ylim(ha4,[0,100]);
    ylim(ha4,fRange)
    xlabel(ha4,'Click number')
    ylabel(ha4,'Frequency (kHz)')
    title(ha4,'Clicks sorted by peak frequency','FontWeight','bold')
    ha4.NextPlot = 'add';
    ha4.Tag = 'ConcatenatedSpectrogramAxes';
end


%% plotSpectrumOverlay ----------------------------------------------------
function plotSpectrumOverlay(hf,data,refSpecData,ii,nEvents,beakTime,nBeakedOriginal)
% Does just that. The code for this part has not changed much, I mostly
% just wrapped it into a function. The reference spectra data have been
% structified to avoid passing a huge number of input arguments (or 
% repeatedly loading the spectra MAT files).
% I think the reference spectra system and the overlay plot itself would 
% benefit from a good cleanup, but that's for another time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % extract variables
    f = refSpecData.f;
    refSpec = refSpecData.spec;

    % get mean spectrum
    meanSpec = TWD_Common.computeMeanSpectrum(data.specClick);
    
    %normalize mean spectrum for detection event
    [~, fbin1] = min(abs(f.fAMAR-10)); %%% index of frequency in f.fAMAR that is closest to 10
    [~, fbin2] = min(abs(f.fAMAR-35));
    [~, fbin3] = min(abs(f.fAMAR-20));
    [~, fbin4] = min(abs(f.fAMAR-80));
    meanSpec = meanSpec - min(meanSpec(fbin1:fbin2)); %%% translate down so that 0 = min between 10 and 35 kHz
    meanSpec = meanSpec/max(meanSpec(fbin3:fbin4)); %%% normalize to max within 20-80 kHz
    %%% I guess this is done to stay consistent with how refSpecs were computed?

    % prepare overlay figure window
    set(0,'CurrentFigure',hf)
    clf
    hf.Name = sprintf('UBW %s/%s',num2str(ii),num2str(nEvents));
    axes(hf,'NextPlot','add');

    plot(f.fHARP,refSpec.meanSpecMe,'b','LineWidth',2)
    plot(f.fHARP,refSpec.meanSpecZc,'g','LineWidth',2)
    plot(f.fHARP,refSpec.meanSpecMd,'r','LineWidth',2)
    plot(f.fHARP,refSpec.meanSpecBW31,'m','LineWidth',2)
    plot(f.fHARP,refSpec.meanSpecBW38,'k','LineWidth',2)
    plot(f.fAMARref,refSpec.meanSpecHa,'y','LineWidth',2)
    plot(f.fARRAY,refSpec.meanSpecMm,'--b','LineWidth',2)
    plot(f.fARRAY,refSpec.meanSpecMb,'c','LineWidth',2)

    plot(f.fAMAR,meanSpec,'k','Linewidth',4)

    % annotate plot with duration of beaked whale detections
    text(0.75,0.9,sprintf('TimeWithBeaked = %s',char(beakTime)),'Unit','normalized')
    
    % annotate plot with number of clicks
    %%% If clicks have been culled, the original number is displayed in
    %%% brackets
    nClicksBeaked = numel(data.dur);
    nClicksStr = sprintf('# clicks = %d',nClicksBeaked);
    if nBeakedOriginal > nClicksBeaked
        nClicksStr = [nClicksStr,sprintf(' (%d)',nBeakedOriginal)];
    end
    text(0.8,0.8,nClicksStr,'Unit','normalized')

    ylim([-0.2 1.5])
    xlim([5 95])% JS?: "changed from 80"

    legend('Me','Zc','Md','BW31','BW38','Ha','Mm','Mb','Location','EastOutside')
end


%% plotClick --------------------------------------------------------------
function plotClick(hfClick,hfClickSeq,data,clickStarts,sortOpt,iSort,img,Fs,fRange,figNameBase,caller)
% Plots the waveform and spectrogram of a single click, and also displays 
% the occurence and maximum amplitude of nearby clicks as a stem plot.
% Partially based on based on 
% "beaked_automatic_overlay_combined_det_AMAR250"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% winStartEndTime = 2-element vector, in milliseconds:
%   - 1st element = buffer before start of click (formerly "s")
%   - 2nd element = end of window (formerly "e")
% for s = 101, winStartEndTime(1) = 0.4;
% for e = 600 (BWD original), winStartEndTime(2) = 2.396;
% for e = 800 (SWD original), winStartEndTime(2) = 3.196;


    % define, extract, or compute common variables
    iClickTarget = iSort(img); % index of target click
    %tClickTarget = data.pos(iClickTarget,1);
    dtClickTarget = clickStarts(iClickTarget);
    nClicksTotal = numel(data.dur);
    localWinLengthHalf = 2; % half length of search window around target click, in seconds (for finding neighbouring clicks)
    offset = data.offset;
    
    % identify clicks around target
    %tClicksLocal_Relative = data.pos(:,1) - tClickTarget;
    tClicksLocal_Relative = seconds(clickStarts - dtClickTarget);
    iClicksLocal = find(tClicksLocal_Relative >= -localWinLengthHalf & tClicksLocal_Relative <= localWinLengthHalf);
    iClicksLocal = [iClicksLocal(1)-1; iClicksLocal; iClicksLocal(end)+1];
    iClicksLocal = iClicksLocal(iClicksLocal > 0 & iClicksLocal <= nClicksTotal);
    tClicksLocal_Relative = tClicksLocal_Relative(iClicksLocal);
    isTarget_Relative = tClicksLocal_Relative == 0;

    % extract target and local clicks
    xClickTarget = data.yFilt(iClickTarget,:);
    durClickTarget = data.dur(iClickTarget);
    xClicksLocal = data.yFilt(iClicksLocal,:);
    xNoiseLocal = data.yNFilt(iClicksLocal,:);
    
    % set window range and click plotting range depending on caller
    switch caller
        case 'BWD'
            winStartEndTime = [0.4, 2.396];
            xRange = [0.1, 1.9];
        case 'SWD'
            winStartEndTime = [0.4, 3.196];
            xRange = [0.1, 2.7];
    end
    
    % set plotting parameters (adapted from SBP's code)
    %%% waveform
    s = round((winStartEndTime(1)/1000)*Fs) + 1;
    e = round((winStartEndTime(2)/1000)*Fs) + 1;
    dd=e-s+1;
    t=0:(dd/(Fs/1000))/(dd-1):dd/(Fs/1000);
    %%% spectra
    winlength = 0.2; %based on HARP fs 200k, 40 DFT
    olPerc = 39/40; %based on 40 DFT, 98% overlap
    df = round(winlength/1000*Fs); %EDITED - round to integer
    ol = min([round(df*olPerc), df-1]);
    [~,F,T,P] = spectrogram(xClickTarget(s:e),df,ol,df,Fs);
    T = T*1000;
    F = F/1000;
    %FViewMax = min([100,Fs/2000]); % max frequency to display in plot (i.e. yMax)
    
    % set up click range patch
    %%% determine click time range within window
    sc = offset - s + 1;
    sct = (sc/Fs)*1000;
    ect = sct + durClickTarget;
    %%% click patch properties
    patchT = [sct;sct;ect;ect];
    patchCol = [1,0.9,0.75];
    patchEdgeCol = [0.9,0.8,0.6];
    
    % translate all times
    t = t - sct;
    T = T - sct;
    patchT = patchT - sct;
    xRange = xRange - sct;
    
    % find peaks and compute RMS noise for all local clicks
    xPeaks = max(abs(xClicksLocal),[],2);
    xPeakTarget = xPeaks(isTarget_Relative);
    xNoiseRMS = rms(xNoiseLocal,2);
    %xRMSLocal = 10.^(data.rmsSignal(iClicksLocal)./20);
    %xRMSTarget = 10.^(data.rmsSignal(iClickTarget)./20);
    %xRMSNoise = 10.^(data.rmsNoise(iClicksLocal)./20);

    
    % setup figure window
    set(0,'CurrentFigure',hfClick)
    clf
    hfClick.Name = [figNameBase,'_click'];

    % plot waveform
    ha_wav = subplot(2,1,1);
    plot(t,xClickTarget(s:e),'k');
    yMin = ha_wav.YLim(1);
    yMax = ha_wav.YLim(2);
    hold on
    hp = patch(patchT,[yMin,yMax,yMax,yMin],patchCol,'EdgeColor',patchEdgeCol);
    uistack(hp,'bottom')
    ha_wav.Layer = 'top';
    xlim(xRange);
    ylim([yMin,yMax]);
    ylabel(ha_wav,'Amplitude [counts]','fontsize',10,'fontweight','b')
    title(ha_wav,sprintf('click detail # %d/%d (%s)',img,nClicksTotal,sortOpt))

    % plot spectrogram
    ha_spec = subplot(2,1,2);
    surf(T,F,10*log10(abs(P)),'EdgeColor','none');
    axis xy; colormap(TWD_Common.blue); view(0,90); %ylim([0 FViewMax]); 
    xlim(xRange);
    ylim(fRange);
    xlabel(ha_spec,'Time [ms]','fontsize',10,'fontweight','b')
    ylabel(ha_spec,'Frequency [kHz]','fontsize',10,'fontweight','b');
    
    % link waveform and spectrogram
    linkaxes([ha_wav,ha_spec],'x');
    
    
    % setup click sequence figure window
    set(0,'CurrentFigure',hfClickSeq)
    clf
    hfClickSeq.Name = [figNameBase,'_click_sequence'];
    
    % plot local click sequence
    ha_local = axes();
    stem(ha_local,tClicksLocal_Relative,xPeaks)
    %stem(ha_local,tClicksLocal_Relative,xRMSLocal)
    hold on
    stem(ha_local,0,xPeakTarget,'r');
    %stem(ha_local,0,xRMSTarget,'r');
    noiseCol = [0.85, 0.325, 0.098]; % from 'lines' map #2
    area(ha_local,tClicksLocal_Relative,xNoiseRMS,'FaceColor',noiseCol,'EdgeColor',noiseCol,'FaceAlpha',0.5);
    %area(ha_local,tClicksLocal_Relative,xRMSNoise,'FaceColor',noiseCol,'EdgeColor',noiseCol,'FaceAlpha',0.5);
    xlim(ha_local,[-localWinLengthHalf,localWinLengthHalf]);
    
    box('on');
    grid('on');
    xlabel(ha_local,'Time [s]','fontsize',10,'fontweight','b')
    ylabel(ha_local,'Amplitude','fontsize',10,'fontweight','b')
    title(ha_local,'Local Click Sequence')
end


%% showClickParams --------------------------------------------------------
function showClickParams(hf,valTableTarget,paramPassTablesTarget,discrimPassMatTarget,discrimNames,isDiscrimActive,tClickStartTarget,iSort,img,showParamPass,figNameBase)
% Creates a UITable displaying click parameters, and whether or not the
% click passed discrimination criteria for each discriminator.
% Assumes certain data are stored within UserData property of figure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get colours for UITable cells, and text object font size
    colActiveHex = hf.UserData.colActiveHex;
    colInactiveHex = hf.UserData.colInactiveHex;
    colActiveDec = hf.UserData.colActiveDec;
    colInactiveDec = hf.UserData.colInactiveDec;
    txtFontSize = hf.UserData.txtFontSize;

    % isolate data
    iClick = iSort(img);
    %%% features
    paramNames = valTableTarget.Properties.VariableNames';
    paramVals = table2cell(valTableTarget(iClick,:))';
    nParams = numel(paramVals);
    %%% discrimination
    nDiscrim = numel(discrimNames);
    discrimPass = discrimPassMatTarget(iClick,:)';
    discrimPassStr = repmat({'FAIL'},nDiscrim,1);
    discrimPassStr(discrimPass) = {'PASS'};
    
    % if pass status for each discriminator should be displayed, get that
    % information
    if showParamPass
    
        % define anonymous function for creating coloured cells in a UITable
        % - works by exploiting UITable's undocumented support for HTML code; 
        % creates a coloured, single-cell HTML table to be inserted within a
        % cell of the UI table (thanks Friedrich from MATLAB Answers)
        htmlCell = @(txt,col) ['<html><table border=0 width=400 bgcolor=',col,'><TR><TD>',txt,'</TD></TR> </table></html>'];

        % process pass status for individual criteria
        % ultimate result is a table of pass/fail/NA with rows = features and
        % columns = discriminators
        paramPassCell = cell(nParams,nDiscrim);
        for ii = 1:nDiscrim
            paramPassTableii = paramPassTablesTarget{ii};
            if isDiscrimActive(ii)
                cellCols = colActiveHex;
            else
                cellCols = colInactiveHex;
            end

            for jj = 1:nParams
                paramjj = paramNames{jj};

                % check if current discriminator assesses current parameter
                if ismember(paramjj,paramPassTableii.Properties.VariableNames)
                    if paramPassTableii.(paramjj)(iClick)
                        strij = 'PASS';
                        colij = cellCols.pass; 
                    else
                        strij = 'FAIL';
                        colij = cellCols.fail; 
                    end
                else
                    strij = 'N/A';
                    colij = cellCols.unassessed;
                end

                % store loop results
                paramPassCell{jj,ii} = htmlCell(strij,colij);
            end
        end
        
        % set up other info for UI table
        paramPassColNames = discrimNames';
        paramTableColWidths = hf.UserData.paramTableColWidths;
        
    else
        % don't show param pass status
        paramPassCell = cell(1,0);
        paramPassColNames = cell(1,0);
        paramTableColWidths = hf.UserData.paramTableColWidths(1:2);
    end
        
    % compile data
    paramData = [paramNames,paramVals,paramPassCell];
    discrimData = [discrimNames,discrimPassStr];
    
    % prepare background colours for discrim table
    discrimPassActive = discrimPass & isDiscrimActive;
    discrimFailActive = ~discrimPass & isDiscrimActive;
    discrimPassInactive = discrimPass & ~isDiscrimActive;
    discrimFailInactive = ~discrimPass & ~isDiscrimActive;
    
    discrimColMat_PassActive = repmat(colActiveDec.pass,sum(discrimPassActive),1);
    discrimColMat_FailActive = repmat(colActiveDec.fail,sum(discrimFailActive),1);
    discrimColMat_PassInactive = repmat(colInactiveDec.pass,sum(discrimPassInactive),1);
    discrimColMat_FailInactive = repmat(colInactiveDec.fail,sum(discrimFailInactive),1);
    
    discrimCol = zeros(nDiscrim,3);
    discrimCol(discrimPassActive,:) = discrimColMat_PassActive;
    discrimCol(discrimFailActive,:) = discrimColMat_FailActive;
    discrimCol(discrimPassInactive,:) = discrimColMat_PassInactive;
    discrimCol(discrimFailInactive,:) = discrimColMat_FailInactive;
    
    % determine click occurance information
    tClickStart = tClickStartTarget(iClick);
    strDate = char(tClickStart,'yyyy-MM-dd');
    strTime = char(tClickStart,'HH:mm:ss.SSS');
    %strRecTime = char(tClickStartRec,'hh:mm:ss.SSS'); % !event merging invalidates relative time as a useful measure
    occurData = {...
        'Date', strDate;...
        'Time',  strTime};
    
    % update figure window
    set(0,'CurrentFigure',hf)
    clf
    hf.Name = [figNameBase,'_click_params'];
    

    % create UI components
    uiCell = cell(6,1);
    %%% txtOccur
    uiCell{1} = uicontrol('Style','text',...
            'Position',hf.UserData.pos_txtOccur,...
            'String','Occurrence',...
            'FontSize',txtFontSize,...
            'HorizontalAlignment','left');
    %%% tableOccur
    uiCell{2} = uitable(...
        'Position',hf.UserData.pos_tableOccur,...
        'Data', occurData,...
        'RowName',[],...
        'ColumnName',[],...
        'ColumnWidth',hf.UserData.occurTableColWidths,...
        'RowStriping','off');
    %%% txtDiscrim
    uiCell{3} = uicontrol('Style','text',...
            'Position',hf.UserData.pos_txtDiscrim,...
            'String','Discrimination',...
            'FontSize',txtFontSize,...
            'HorizontalAlignment','left');
    %%% tableDiscrim
    uiCell{4} = uitable(...
        'Position',hf.UserData.pos_tableDiscrim,...
        'Data', discrimData,...
        'RowName',[],...
        'ColumnName',[],...
        'ColumnWidth',hf.UserData.discrimTableColWidths,...
        'BackgroundColor',discrimCol);
    %%% txtParam
    uiCell{5} = uicontrol('Style','text',...
            'Position',hf.UserData.pos_txtParam,...
            'String','Parameters',...
            'FontSize',txtFontSize,...
            'HorizontalAlignment','left');
    %%% tableParam
    uiCell{6} = uitable(...
        'Position',hf.UserData.pos_tableParam,...
        'Data', paramData,...
        'RowName',[],...
        'ColumnName',[{'Feature','Value'},paramPassColNames],...
        'ColumnWidth',paramTableColWidths);
    
    % change units of everything to normalized
    %for ii = 1:numel(uiCell)
    %    uiCell{ii}.Units = 'normalized';
    %end
end


%% highlightConSpecClick --------------------------------------------------
function hp = highlightConSpecClick(haConSpec, iClickConSpec, patchCol, patchTag, nClicks)
% Highlight a click on the concatenated spectrogram plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define transparency as a function of number of clicks (using a
    % transformed logistic function). This is done because highlighting 
    % can be harder to see when there are many clicks.
    patchAlphaMin = 0.2;
    patchAlphaMax = 0.8;
    patchAlpha = (1./(1+exp((250-nClicks)./50))).*(patchAlphaMax-patchAlphaMin) + patchAlphaMin;
    
    % define patch properties
    %patchAlpha = 0.25;
    xPatch = iClickConSpec + [-0.5; -0.5; 0.5; 0.5];
    yPatch = [haConSpec.YLim(1); haConSpec.YLim(2); haConSpec.YLim(2); haConSpec.YLim(1)];
    
    % delete existing patch if there is one
    hp = findobj(haConSpec, 'Type','patch', 'Tag',patchTag);
    if ~isempty(hp)
        delete(hp)
    end

    % draw new patch
    hp = patch(haConSpec, xPatch, yPatch, patchCol, 'EdgeColor','none', 'FaceAlpha',patchAlpha, 'Tag',patchTag);
end

%% saveTable --------------------------------------------------------------
function saveTable(xlsOutData,xlsOutPath,wipe)
% Attempts to save species ID spreadsheet.
% If the save fails, user is prompted to retry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    saveComplete = false;
    while ~saveComplete
        try
            if isfile(xlsOutPath) && wipe
                delete(xlsOutPath);
            end
            writetable(xlsOutData,xlsOutPath)
            [~,xlsFileName,xlsExt] = fileparts(xlsOutPath);
            fprintf('Spreadsheet "%s" saved\n',[xlsFileName,xlsExt]);
            %disp('Spreadsheet saved')
            saveComplete = true;
        catch ME
            % if save fails, issue warning dialog and possibility to retry
            OKPressed = false;
            warnStr = sprintf('Could not save spreadsheet:\n\n%s\n\nPress OK to retry or X to cancel.',ME.message);
            hfWarn = warndlg(warnStr,'Warning: save failed','modal');
            
            % overwrite callback function for OK button
            hpbOK = findobj(hfWarn,'Tag','OKButton');
            hpbOK.Callback = @warnOK_Callback;
            
            % wait for user action
            waitfor(hfWarn)
            
            % check if OK button was pressed. If not, cancel the save.
            if ~OKPressed
                saveComplete = true;
            else
                disp('Save canceled')
            end
        end
    end
    
    % nested function: OK button callback .................................
    function warnOK_Callback(src,edata)
        OKPressed = true;
        delete(gcbf)
    end
end


%% saveTableAs ------------------------------------------------------------
function [filePath,saveCanceled] = saveTableAs(xlsOutData,dirPath_results,depName,targetName)
% Ask user to specify a spreadsheet name for saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define base file name
    fileNameBase = sprintf('%s_%s_Validated',depName,targetName);
    
    % user input loop
    prompt = sprintf('Saving new "%s" spreadsheet...\nEnter a unique string to append to filename, or leave blank:',fileNameBase);
    saveState = 'processing';
    while strcmp(saveState,'processing')
        % prompt user to specify filename
        suffCell = inputdlg(prompt,'New Spreadsheet');
        
        if isempty(suffCell)
            % cancel pressed - abort
            saveState = 'cancel';
            continue
        else
            % create file path
            suffStr = suffCell{1};
            if isempty(suffStr)
                fileName = [fileNameBase,'.xlsx'];
            else
                fileName = [fileNameBase,'_',suffStr,'.xlsx'];
            end
            filePath = fullfile(dirPath_results,fileName);

            % check if file already exists, and prompt user for next action 
            % if it does
            if isfile(filePath)
                owPrompt = sprintf('File "%s" already exists. Overwrite?',fileName);
                owOpt = questdlg(owPrompt,'Warning: file exists','Yes','No','Cancel','Yes');
                switch owOpt
                    case 'Yes'
                        % proceed with save
                        saveState = 'dosave';

                    case 'No'
                        % ask for new string
                        continue
                        
                    case {'Cancel',''}
                        % cancel
                        saveState = 'cancel';
                        
                    otherwise
                        error('Encountered bug: option not recognized')
                end
            else
                % file doesn't yet exist - proceed with save
                saveState = 'dosave';
            end
        end
    end
    
    % evaluate state and take appropriate action
    switch saveState
        case 'dosave'
            % save file
            saveTable(xlsOutData,filePath,true);
            saveCanceled = false;
            %fprintf('Saved data to new spreadsheet: "%s"\n',fileName)
            
        case 'cancel'
            filePath = '';
            saveCanceled = true;
            disp('Save canceled')
    end
end


%% CLICK PICKER CALLBACKS =================================================
% FigButtonMotionFcn_summary ..............................................
function FigButtonMotionFcn_summary(src, edata)
% Highlight the spectrum over which the mouse cursor is hovering on the
% concatenated spectrogram plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get cursor position relative to concatenated spectrogram axis
    ha = findobj(src, 'Type','axes', 'Tag','ConcatenatedSpectrogramAxes');
    
    cp = ha.CurrentPoint;
    x = cp(1,1);
    y = cp(1,2);
    
    % find and delete existing highlight patch if it exists
    patchTag = 'SpecHover';
    hp = findobj(ha, 'Type','patch', 'Tag',patchTag);
    if ~isempty(hp)
        delete(hp);
    end
    
    % determine if cursor is in range or not
    withinX = x >= ha.XLim(1) && x <= ha.XLim(2);
    withinY = y >= ha.YLim(1) && y <= ha.YLim(2);
    
    % draw patch to highlight the area over which cursor is hovering
    if withinX && withinY
        patchCol = [1, 0, 0];
        
        highlightConSpecClick(ha, round(x), patchCol, patchTag, floor(ha.XLim(2)));
    end
end


% FigButtonDownFcn_summary ................................................
function FigButtonDownFcn_summary(src, edata)
% Select a click on the concatenated spectrogram plot based on cursor
% location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % remove the motion callback (better to do it here for debugging)
    src.WindowButtonMotionFcn = '';

    % return highlight patch if there is one
    hp = findobj(src, 'Type','patch', 'Tag','SpecHover');
    
    % get index of click being hovered over (if any)
    if ~isempty(hp)
        src.UserData = mean(hp.XData);
        delete(hp);
    else
        src.UserData = [];
    end
end


%% DISCRIMINATOR SWITCHER =================================================
% createDiscrimWindow .....................................................
function [hf,ht] = createDiscrimSwitchWindow(discrimNames,discrimDefaults)
% creates window for toggling discriminators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define variables
    pad = 10;
    pbWidth = 50;
    pbHeight = 25;
    yGap = 10;
    txtHeight = 20;
    
    % OK button position
    posPushbuttonOK = [...
        pad,...
        pad,...
        pbWidth,...
        pbHeight];
    
    % table position
    tColWidths = [60,60,60];
    posTable = [...
        pad,...
        pad + pbHeight + yGap,...
        sum(tColWidths) + 100,...
        140];
    
    % text position
    posText = [...
        pad,...
        posTable(2) + posTable(4) + yGap,...
        posTable(3),...
        txtHeight];

    % define window position
    posWin = [...
        600,...
        200,...
        pad*2 + posText(3),...
        pad*2 + pbHeight + yGap*2 + posTable(4) + txtHeight];
    
    % initialize UserData (determins if action proceeds or not)
    UDInitial = struct(...
        'DefaultSettings',[],...
        'CurrentSettings',[],...
        'Proceed',false);
    
    % set default settings
    UDInitial.DefaultSettings = discrimDefaults;
    
    % create figure
    hf = figure(...
        'Position',posWin,...
        'Name','Click Discriminators',...
        'CloseRequestFcn',@discrimUI_CloseReq,...
        'Visible','off',...
        'MenuBar','none',...
        'UserData',UDInitial);
    
    % create OK button
    hpb = uicontrol(...
        'Tag','pushbuttonOK',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posPushbuttonOK,...
        'Style','pushbutton',...
        'String','OK',...
        'FontSize',10,...
        'Callback',@discrimUI_pushbuttonOK_Callback);
    
    % create UITable
    ht = uitable(...
        'Tag','uitableDiscrim',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posTable,...
        'ColumnName',{'Positive','Negative','Disabled'},...
        'ColumnWidth',num2cell(tColWidths),...
        'ColumnEditable',[true,true,true],...
        'RowName',discrimNames,...
        'CellEditCallback',@discrimTable_Callback);
    
    % create header text
    htxt = uicontrol(...
        'Tag','text_Discrim',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posText,...
        'Style','text',...
        'String','Choose discriminator settings',...
        'FontSize',10);
    
end


% discrimUI_CloseReq ......................................................
function discrimUI_CloseReq(src,edata)
% function to execute when user clicks on X button.
% Does not actually delete the window, just hides it and cancels any
% changes to the discriminator settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get child table
    ht = findobj(src,'Tag','uitableDiscrim');

    % reset settings
    ht.Data = src.UserData.CurrentSettings;
    
    % set "Proceed" to false
    src.UserData.Proceed = false;

    % hide window
    src.Visible = 'off';
    
    % resume input
    uiresume(src);
end


% discrimUI_pushbuttonOK_Callback .........................................
function discrimUI_pushbuttonOK_Callback(src,edata)
% Function to call when OK button is pressed.
% Hides discriminator switcher and signals go-ahead for changing the
% settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get window
    hf = src.Parent;
    
    % get table
    ht = findobj(hf,'Tag','uitableDiscrim');
    
    % check if data changed or not
    oldData_Num = cell2mat(hf.UserData.CurrentSettings);
    newData_Num = cell2mat(ht.Data);
    
    % set "Proceed" to true if data has changed
    if any(oldData_Num(:) ~= newData_Num(:))
        hf.UserData.Proceed = true;
    else
        hf.UserData.Proceed = false;
    end
    
    % hide window
    hf.Visible = 'off';
    
    % resume input
    uiresume(hf);
end


% discrimTable_Callback ...................................................
function discrimTable_Callback(src,edata)
% Function to execute when a selection has been made in the discriminator
% settings table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get indices of change
    iRowChanged = edata.Indices(1);
    iColChanged = edata.Indices(2);
    
    % process user input
    if edata.NewData == true
        % case where new box was checked -> deselect all other boxes in the
        % row
        nCol = size(src.Data,2);
        iColOther = 1:nCol;
        iColOther = iColOther(iColOther ~= iColChanged);
        src.Data(iRowChanged,iColOther) = {false};
    else
        % case where user clicked on already selected box -> cancel
        % operation
        src.Data{iRowChanged,iColChanged} = true;
    end
end


%% EVENT SWITCHER =========================================================
% createEventSwitchWindow .................................................
function [hf,ht] = createEventSwitchWindow()
% Creates window for selecting which event to process.
% Allows backtracking to already completed events.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define variables
    pad = 10;
    pbWidth = 50;
    pbHeight = 25;
    yGap = 10;
    pbGap = 20;
    txtHeight = 20;
 
    % pushbutton positions
    %%% OK
    posPushbuttonOK = [...
        pad,...
        pad,...
        pbWidth,...
        pbHeight];
    %%% Last
    posPushbuttonLast = [...
        pad+pbWidth+pbGap,...
        pad,...
        pbWidth,... 
        pbHeight];
    
    % table position
    tColWidths = [110,60,50];
    posTable = [...
        pad,...
        pad+pbHeight + yGap,...
        sum(tColWidths)+64,...
        500];
    
    % text position
    posText = [...
        pad,...
        posTable(2) + posTable(4) + yGap,...
        max([posTable(3),pbWidth*2+pbGap]),...
        txtHeight];

    % define window position
    posWin = [...
        600,...
        200,...
        pad*2 + posText(3),...
        pad*2 + pbHeight + yGap*2 + posTable(4) + txtHeight];
    
    % initialize UserData (determins if action proceeds or not)
    UDInitial = struct('Proceed',false);
    
    % create figure
    hf = figure(...
        'Position',posWin,...
        'Name','Event Switcher',...
        'CloseRequestFcn',@eventSwitchUI_CloseReq,...
        'Visible','off',...
        'MenuBar','none',...
        'UserData',UDInitial);
    
    % create OK button
    hbpOK = uicontrol(...
        'Tag','pushbuttonOK_Event',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posPushbuttonOK,...
        'Style','pushbutton',...
        'String','OK',...
        'FontSize',10,...
        'Callback',@eventSwitchUI_pushbuttonOK_Callback);
    
    % create Last button
    hbpLast = uicontrol(...
        'Tag','pushbuttonLast_Event',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posPushbuttonLast,...
        'Style','pushbutton',...
        'String','Last',...
        'FontSize',10,...
        'Callback',@eventSwitchUI_pushbuttonLast_Callback);
    
    % create event switcher table
    ht = uitable(...
        'Tag','uitableEventSwitch',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posTable,...
        'ColumnName',{'EventStart','Species','Active'},...
        'ColumnWidth',num2cell(tColWidths),...
        'ColumnEditable',[false,false,true],...
        'CellEditCallback',@eventSwitchTable_Callback);
    
    % create header text
    htxt = uicontrol(...
        'Tag','text_Event',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posText,...
        'Style','text',...
        'String','Select Event',...
        'FontSize',10);
end


% eventSwitchUI_CloseReq ..................................................
function eventSwitchUI_CloseReq(src,edata)
% function to execute when user clicks on X button of event switcher 
% window. Does not actually delete the window, just hides it and cancels 
% any action.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set "Proceed" to false
    src.UserData.Proceed = false;

    % hide window
    src.Visible = 'off';
    
    % resume input
    uiresume(src);
end


% eventSwitchUI_pushbuttonOK_Callback .....................................
function eventSwitchUI_pushbuttonOK_Callback(src,edata)
% Function to call when OK button for event switcher is pressed.
% Hides event switcher and signals go-ahead for changing the settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get window
    hf = src.Parent;
    
    % signal go-ahead
    hf.UserData.Proceed = true;
    
    % hide window
    hf.Visible = 'off';
    
    % resume input
    uiresume(hf);
end


% eventSwitchUI_pushbuttonLast_Callback ...................................
function eventSwitchUI_pushbuttonLast_Callback(src,edata)
% Function to call when "Last" button for event switcher is pressed.
% Hides event switcher and signals go-ahead for changing the settings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get window
    hf = src.Parent;
    
    % get event table
    ht = findobj('Parent',hf,'Type','uitable');
    
    % find index of first incomplete event
    sp = cell2mat(ht.Data(:,2));
    n = numel(sp);
    iLeftOff = min([n,find(isnan(sp),1,'first')]); % note this is not necessarily the latest entry
    
    % set incomplete event as the active one
    ht.Data(:,3) = {false};
    ht.Data{iLeftOff,3} = true;
    
    % signal go-ahead
    hf.UserData.Proceed = true;
    
    % hide window
    hf.Visible = 'off';
    
    % resume input
    uiresume(hf);
end


% eventSwitchTable_Callback ...............................................
function eventSwitchTable_Callback(src,edata)
% Function to execute when a selection has been made in the event switcher
% table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get index of target column
    iCol = edata.Indices(2);
    
    % get index of row that changed
    iChanged = edata.Indices(1);
    
    % process user input
    if edata.NewData == true
        % case where new row was selected -> deselect all other rows
        n = size(src.Data,1);
        iOther = 1:n;
        iOther = iOther(iOther ~= iChanged);
        src.Data(iOther,iCol) = {false};
    else
        % case where user clicked on already selected row -> cancel
        % operation
        src.Data{iChanged,iCol} = true;
    end
end


% updateEventSwitchTable ..................................................
function updateEventSwitchTable(ht,xlsOutData,iEvent)
% updates the table in the event switcher window. Call this every time the
% window is open.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % remember, table columns are:
    % EventStart, Species, Active
    
    % initialize
    n = height(xlsOutData);
    sp = xlsOutData.Species;
    
    % get completed events
    hasID = ~isnan(sp);
    
    % determine logical index of current event
    isActive = false(n,1);
    isActive(iEvent) = true;
    
    % create table data
    eStartCell = cellstr(xlsOutData.StartTime);
    spCell = num2cell(sp);
    activeCell = num2cell(isActive);
    dataCell = [eStartCell,spCell,activeCell];
    
    % create colour matrix (incomplete events will be red)
    colID = [1,1,1];
    colNoID = [1,0.75,0.75];
    col = zeros(n,3);
    col(hasID,:) = repmat(colID,sum(hasID),1);
    col(~hasID,:) = repmat(colNoID,sum(~hasID),1);
    
    % update table
    ht.Data = dataCell;
    ht.BackgroundColor = col;
    
end