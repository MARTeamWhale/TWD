%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "editFilterParams"
%   Written by Wilfried Beslin
%   Last updated May. 24, 2023, using MATLAB R2018b
%
%   Description:
%   Use this script to edit the list of cutoff frequencies and filter
%   orders used when filtering clicks during processing. Filters are
%   implemented as digital Butterworth bandpass filters.
%
%   Cutoff frequencies are dependent on sampling rate. Cutoffs that work 
%   well for one sampling rate may not work as well (or at all) for a 
%   different rate. Therefore, cutoff frequencies should be specified 
%   explicitly for each sampling rate that may be encountered. If working 
%   with new data that was sampled at a previously unused sampling rate, 
%   use this script to add cutoffs for the new rate. Click compilation 
%   cannot proceed on the new data otherwise.
%
%   Note that filter parameters will be sorted by sampling rate once saved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% Known Issues:
% - First initialization following the startup of a MATLAB session may take
% several seconds. Profiling has revealed that the bottleneck is with the
% 'lp2bp' function, which is one of the internal functions called by
% 'butter'. Nothing can be done about it. 
% - Likewise, fvtool can sometimes take a few seconds before showing up
% after clicking a VIEW cell. Popup messages can also have slightly longer-
% than-usual delays from time to time. These seems more unpredictable.

function editFilterParams()

    % determine path to code directory
    dirPath_root = mfilename('fullpath');
    [dirPath_root,~,~] = fileparts(dirPath_root);
    [dirPath_root,~,~] = fileparts(dirPath_root);
    
    % look for module folders
    [~, modFolders] = TWD_Common.Utilities.listFiles(dirPath_root, '[folders]', 'Recursive',false, 'NameMustContain','^_');
    modNames = strrep(modFolders, '_', '');

    % prompt user to specify the module being used
    [iMod, madeSelection] = listdlg(...
        'ListString',modNames,...
        'PromptString','Edit filters for which module?',...
        'SelectionMode','single',...
        'ListSize',[150,60],...
        'Name','Module Selection');
    if ~madeSelection
        return
    end
    module = modNames{iMod};
    modFolder = modFolders{iMod};
    
    % read filter data file
    filtFilePath = fullfile(dirPath_root,modFolder,'BandpassFilterParams.mat');
    load(filtFilePath, 'filtdata')
    
    % create UI for editing filters and initalize it
    disp('Initializing...')
    [hf, ht] = createFilterEditWindow(module, filtdata);
    initializeFiltEditTable(ht, filtdata, module, filtFilePath);
    disp('Done')
    
end


%% createFilterEditWindow -------------------------------------------------
function [hf, ht] = createFilterEditWindow(module, filtdata)
% Creates window for editing bandpass filter cutoff frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define variables
    pad = 20;
    pbWidth = 150;
    pbHeight = 50;
    pbGap = 30;
    pbFontSize = 10;
    txtFontSize = 9;
    titFontSize = 11;
    titHeight = 30;
    figWidth = 1200;
    figHeight = 500;
    
    % POSITIONS ...........................................................
    
    % table position
    tColWidths = [100,80,80,60,60,40,250];
    posTable = [...
        pad,...
        pad*2 + pbHeight,...
        sum(tColWidths)+50,...
        figHeight - (pad*3 + pbHeight)];
    
    % pushbutton positions
    %%% Save
    posPushbuttonSave = [...
        pad*2,...
        pad,...
        pbWidth,...
        pbHeight];
    %%% Quit
    posPushbuttonQuit = [...
        posPushbuttonSave(1) + posPushbuttonSave(3) + pbGap,...
        pad,...
        pbWidth,... 
        pbHeight];
    
    % text position
    %%% Title
    posTextTitle = [...
        posTable(1) + posTable(3) + pad,...
        figHeight - (pad/2 + titHeight),...
        figWidth - (posTable(1) + posTable(3) + pad*2),...
        titHeight];
    %%% Description
    posTextDescription = [...
        posTextTitle(1),...
        pad,...
        posTextTitle(3),...
        figHeight - (pad*2 + posTextTitle(4))];

    % define window position
    posWin = [...
        600,...
        200,...
        figWidth,...
        figHeight];
    
    % CREATE OBJECTS ......................................................
    
    % create figure
    hf = figure(...
        'Position',posWin,...
        'Name','Bandpass Filter Editor',...
        'NumberTitle','off',...
        'CloseRequestFcn',@filtEditUI_CloseReq,...
        'Visible','on',...
        'MenuBar','none');
    
    % create table
    ht = uitable(...
        'Tag','uitableFiltEdit',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posTable,...
        'ColumnName',[filtdata.Properties.VariableNames, {'IsStable','VIEW','Error'}],...
        'ColumnWidth',num2cell(tColWidths),...
        'ColumnEditable',[true,true,true,true,false,true,false],...
        'CellEditCallback',@filtEditTable_Callback);
    
    % create Save button
    hbpSave = uicontrol(...
        'Tag','pushbuttonSave',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posPushbuttonSave,...
        'Style','pushbutton',...
        'String','Save',...
        'FontSize',pbFontSize,...
        'Callback',@pushbuttonSave_Callback);
    
    % create Quit button
    hbpQuit = uicontrol(...
        'Tag','pushbuttonQuit',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posPushbuttonQuit,...
        'Style','pushbutton',...
        'String','Quit',...
        'FontSize',pbFontSize,...
        'Callback',@pushbuttonQuit_Callback);
    
    % create title text
    htxtTitle = uicontrol(...
        'Tag','textTitle',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posTextTitle,...
        'Style','text',...
        'String',sprintf('Bandpass Filter Settings (%s)',module),...
        'FontSize',titFontSize,...
        'FontWeight','bold');
    
    % create description text
    descriptionString = [...
        {['Use this tool to add or edit digital filters that will ',...
        'be applied to audio data when running ',module,'. A filter ',...
        'must be specified for every potential sampling rate. All ',...
        'filters are designed as Butterworth bandpass filters.']};...
        {''};...
        {''};...
        {'Directions:'};...
        {''};...
        {['- Edit the “SamplingRate”, “CutoffFreq1”, "CutoffFreq2”, ',...
        'and "Order" fields directly as needed']};...
        {''};...
        {['- Enter values in the last row to add a new filter for a ',...
        'new sampling rate']};...
        {''};...
        {'- Delete a row by removing all data from it'};...
        {''};...
        {['- Rows that are incomplete or contain invalid entries will ',...
        'be highlighted in red with an error message']};...
        {''};...
        {['- Click a cell under the "VIEW" column to view properties ',...
        'of the filter that will be constructed using the settings ',...
        'specified in that row. Note that filters for multiple rows ',...
        'can be viewed at once as different tabs, but the tab ',...
        'sequence may not necessarily match the order of the rows in ',...
        'the table!']};...
        {''};...
        {['- The "IsStable" column indicates if the filter is ',...
        'reliable or not. Unstable filters may result in significant ',...
        'artefacts if they are applied to real data and should not ',...
        'be used.']};...
        {''};...
        {['- Press the "Save" button to save changes. Note that this ',...
        'will not work if there are any rows with error messages.']}];
    htxtDescription = uicontrol(...
        'Tag','textDescription',...
        'Parent',hf,...
        'Units','pixels',...
        'Position',posTextDescription,...
        'Style','text',...
        'String',descriptionString,...
        'FontSize',txtFontSize,...
        'HorizontalAlignment','left');
    
end


%% initializeFiltEditTable ------------------------------------------------
function initializeFiltEditTable(ht, filtdata, module, filtFilePath)
% Sets up the FiltEdit table based on available data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numFilts = height(filtdata);

    % initialize UserData
    rowCols = struct(...
        'odd', [1, 1, 1],...
        'even', [0.94, 0.94, 0.94],...
        'bad', [1, 0.75, 0.75],...
        'last',[1, 1, 0.75]);
    cellColsHex = struct(...
        'good', '#4DB34D',... % Hex for [0.3, 0.7, 0.3]; green
        'bad', '#FF4D4D'); % Hex for [1.0, 0.3, 0.3]; red
    
    ht.UserData = struct(...
        'module', module,...
        'filtfile', filtFilePath,...
        'filtdata', {filtdata},...
        'filters', {cell(numFilts,1)},...
        'fvtools', {cell(numFilts,1)},...
        'rowcols', rowCols,...
        'stablecellcols', cellColsHex);

    % fill table
    dataCell = [table2cell(filtdata), repmat({'', false, ''},numFilts,1)];    
    ht.Data = dataCell;
    
    % initialize row colours
    ht.BackgroundColor = ones(numFilts,3);
    
    % create digital filters, evaluate filter stability, and set colours
    for ii = 1:numFilts
        processRow(ht, ii);
    end
    
    % add blank row at the bottom
    addBlankRow(ht);
end


%% updateFiltEditTable ----------------------------------------------------
function updateFiltEditTable(ht, iRow, iCol, oldVal)
% Process a change to a filter parameter in the table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % remember, table columns are:
    % SamplingRate, CutoffFreq1, CutoffFreq2, Order, IsStable, VIEW, Error
    
    userVal = ht.Data{iRow, iCol};
    
    % ensure the changed cell is numeric
    if ischar(userVal)
        ht.Data{iRow, iCol} = str2double(userVal);
    elseif isnan(userVal)
        ht.Data{iRow, iCol} = [];
    end
    newVal = ht.Data{iRow, iCol};
    
    % if Fs was modified, make sure the same value does not already exist;
    % if it does, cancel the operation.
    if iCol == 1
        numEntries = size(ht.Data,1) - 1;
        otherRows = true(numEntries, 1);
        otherRows(iRow) = false;

        fsExists = ismember(newVal,[ht.Data{otherRows, 1}]);
        if fsExists
            warndlg('Cannot have more than one entry with the same SamplingRate', 'Warning', 'modal');
            ht.Data{iRow, iCol} = oldVal;
            return
        end
    end
    
    % check if last row was edited or not. It it was, add new entry.
    lastRowEdited = iRow == size(ht.Data, 1);
    if lastRowEdited
        ht.UserData.filters = [ht.UserData.filters; {{}}];
        ht.UserData.fvtools = [ht.UserData.fvtools; {[]}];
        addBlankRow(ht);
        ht.Data{iRow,6} = false;
    end
    
    % process row
    processRow(ht, iRow)
end


%% processRow -------------------------------------------------------------
function processRow(ht, iRow)
% Process entries in a certain row. Includes creation and assessment of
% filters, setting of row colours and error messages, or deletion of row if
% it is empty.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check if row is empty
    emptyRow = cellfun('isempty', ht.Data(iRow,1:4));
    if emptyRow
        % delete row and associated filter
        ht.Data(iRow,:) = [];
        ht.BackgroundColor(end-1,:) = [];
        ht.UserData.filters(iRow) = [];
        ht.UserData.fvtools(iRow) = [];
        recolourRows(ht);
        
        % update the index within fvtools for all those that came after the
        % deleted row
        numFiltsNew = size(ht, 1) - 1;
        for ii = iRow:numFiltsNew
            toolii = ht.UserData.fvtools{ii};
            if ~isempty(toolii)
                toolii.UserData.idx = toolii.UserData.idx - 1;
            end
        end
    else
        % if fvtool is active, deactivate it
        activeFVTool = ht.Data{iRow, 6};
        if activeFVTool
            close(ht.UserData.fvtools{iRow})
        end
        
        % create and evaluate filter for the modified entry
        errMsg = setFilter(ht, iRow);
        ht.Data{iRow,end} = errMsg;
        
        % recolour rows
        recolourRows(ht);
        
        % restore fvtool if needed
        if activeFVTool
            ht.Data{iRow, 6} = true;
            processFVTool(ht, iRow);
        end
    end
end


%% addBlankRow ------------------------------------------------------------
function addBlankRow(ht)
% Add a new empty row at the bottom of the table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ht.Data = [ht.Data; {[], [], [], [], '', '', ''}];
    
    % recolour
    recolourRows(ht)
end


%% recolourRows -----------------------------------------------------------
function recolourRows(ht)
% Recolour all rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numFilts = size(ht.Data, 1) - 1;

    % get row indices
    oddRows = 1:2:numFilts;
    evenRows = 2:2:numFilts;
    badRows = ~cellfun('isempty', ht.Data(1:numFilts,end));
    
    % create colour matrix
    rowColMat = ones(numFilts, 3);
    rowColMat(oddRows,:) = repmat(ht.UserData.rowcols.odd, numel(oddRows), 1);
    rowColMat(evenRows,:) = repmat(ht.UserData.rowcols.even, numel(evenRows), 1);
    rowColMat(badRows,:) = repmat(ht.UserData.rowcols.bad, sum(badRows), 1);
    rowColMat = [rowColMat; ht.UserData.rowcols.last];
    
    % assign matrix
    ht.BackgroundColor = rowColMat;
end


%% setFilter --------------------------------------------------------------
function errMsg = setFilter(ht, idx)
% Create a digital filter, evaluate its stability, and update the relevant
% area of the UI table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % extract variables
    Fs = ht.Data{idx,1};
    Fc1 = ht.Data{idx,2};
    Fc2 = ht.Data{idx,3};
    ord = ht.Data{idx,4};
    
    % validate values
    errMsg = '';
    checkInt = @(val) validateattributes(val, {'numeric'}, {'positive','integer','scalar'});
    try
        checkInt(Fs);
        checkInt(Fc1);
        checkInt(Fc2);
        checkInt(ord);
        try
            assert(Fc2 < Fs/2, 'CutoffFreq2 must be less than SamplingRate/2');
            assert(Fc1 < Fc2, 'CutoffFreq1 must be less than CutoffFreq2');
            assert(mod(ord,2) == 0, 'Order must be an even number')
        catch ME
            errMsg = ME.message;
        end
    catch
        errMsg = 'Missing or invalid value(s)';
    end

    if isempty(errMsg)
        % build filter based on module
        switch ht.UserData.module
            case 'BWD'
                makeFilt = @makeFilt_butter;
            case 'SWD'
                makeFilt = @makeFilt_designfilt;
            otherwise
                error('Unknown module')
        end

        try
            f = makeFilt(Fs, Fc1, Fc2, ord);
            errMsg = '';
        catch
            f = {};
            errMsg = 'Error creating filter';
        end
    else
        f = {};
    end
    
    % set filter data to UserData
    ht.UserData.filters{idx} = f;
    
    % evaluate filter stability, if there is one
    if isempty(errMsg)
        s = isstable(f{:});
        if s
            stableStr = 'YES';
            stableCol = ht.UserData.stablecellcols.good;
        else
            stableStr = 'NO';
            stableCol = ht.UserData.stablecellcols.bad;
        end
        
        % set cell text and colour by exploiting UITable's undocumented
        % support for HTML code; creates a coloured, single-cell HTML table
        % to be inserted within a cell of the UI table (thanks Friedrich
        % from MATLAB Answers)
        ht.Data{idx,5} = ['<html><table border=0 width=400 bgcolor=',stableCol,'><TR><TD>',stableStr,'</TD></TR> </table></html>'];
    else
        ht.Data{idx,5} = '';
    end
    
    
    % NESTED FUNCTIONS ....................................................
    
    % Create a filter using the 'butter' function
    function f = makeFilt_butter(Fs, Fc1, Fc2, ord)
        [B,A] = butter(ord/2, [Fc1,Fc2]/(Fs/2));
        f = {B, A};
    end

    % Create a filter using the 'designfilt' function
    function f = makeFilt_designfilt(Fs, Fc1, Fc2, ord)
        f = {designfilt('bandpassiir',...
                'DesignMethod','butter',...
                'SampleRate',Fs,...
                'FilterOrder',ord,...
                'HalfPowerFrequency1',Fc1,...
                'HalfPowerFrequency2',Fc2)};
    end
end


%% processFVTool ----------------------------------------------------------
function processFVTool(ht, idx)
% Create or delete fvtool instances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check if VIEW box is ticked or not
    activateTool = ht.Data{idx,6};
    
    if activateTool
        % enable the tool if possible
        f = ht.UserData.filters{idx};
        if ~isempty(f)
            % create tool for this filter
            hfvt = fvtool(f{:});
            hfvt.UserData = struct(...
                'ht', ht,...
                'idx', idx);
            hfvt.CloseRequestFcn = @fvtool_CloseReq;
            
            % set Fs manually if needed
            if strcmp(ht.UserData.module, 'BWD')
                hfvt.Fs = ht.Data{idx,1};
            end
            
            % add tool to container
            ht.UserData.fvtools{idx} = hfvt;
        else
            % filter does not exist, so cancel operation
            ht.Data{idx,6} = false;
        end
    else
        % check if there is an active tool or not
        hfvt = ht.UserData.fvtools{idx};
        if isempty(hfvt)
            % no active tool; this can happen if the user clicks too fast.
            % In this case, restore checkmark.
            ht.Data{idx,6} = true;
        else
            % delete active tool
            close(ht.UserData.fvtools{idx})
            ht.UserData.fvtools{idx} = [];
        end
    end
end


%% createFiltDataTable ----------------------------------------------------
function [filtdata, errMsg] = createFiltDataTable(ht)
% Create a filtdata table from the UITable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    errMsg = '';
    try
        % NOTE: it is better to use cell2mat followed by array2table,
        % rather than a single call to cell2table. This method ensures that
        % the resulting table is completely numeric.
        filtdataMat = cell2mat(ht.Data(1:(end-1),1:4));
        filtdata = array2table(filtdataMat, 'VariableNames',ht.UserData.filtdata.Properties.VariableNames);
        filtdata = sortrows(filtdata,1);
    catch ME
        filtdata = table.empty(0,0);
        errMsg = ME.message;
    end
end


%% saveFiltData -----------------------------------------------------------
function success = saveFiltData(ht)
% Save filtdata to file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check that table has no issues
    dataOK = all(cellfun('isempty',ht.Data(1:(end-1),end)));
    if dataOK
        % check if there any unstable filters
        %unstableFilters = any(strcmp(ht.Data(1:(end-1),5), 'NO')); % this won't work with HTML strings
        filts = ht.UserData.filters;
        numFilt = numel(filts);
        idx = 1;
        unstableFilters = false;
        while ~unstableFilters && idx <= numFilt
            currentFilt = filts{idx};
            unstableFilters = ~isstable(currentFilt{:});
            idx = idx + 1;
        end
        if unstableFilters
            promptStr = sprintf('WARNING: there are unstable filters.\nUsing these filters can result in unreliable data.\nAre you sure you want to save anyway?');
            opt = questdlg(promptStr, 'Unstable Filters Detected', 'Yes', 'No', 'No');
            if any(strcmp(opt, {'No',''}))
                success = false;
                return
            end
        end
        
        % do save
        try
            % create new filtdata
            [filtdata, errMsg] = createFiltDataTable(ht);
            assert(isempty(errMsg), errMsg);
            
            % assign filtdata to UITable
            ht.UserData.filtdata = filtdata;
            
            % save filtdata
            filePath = ht.UserData.filtfile;
            save(filePath, 'filtdata')
            
            % show message to user
            msgbox('Filter data saved successfully', 'Data Saved', 'help', 'modal')
            success = true;
        catch ME
            errMsg = sprintf('Failed to save filter data:\n%s', ME.message, 'modal');
            msgbox(errMsg, 'Save Failed', 'error', 'modal')
            success = false;
        end
    else
        % return warning
        msgbox(sprintf([...
            'There are unresolved issues with the filter settings.\n',...
            'These issues need to be fixed before data can be saved.']),...
            'Cannot Save', 'warn', 'modal');
        success = false;
    end
end


%% CALLBACKS --------------------------------------------------------------

% filtEditTable_Callback ..................................................
function filtEditTable_Callback(src, edata)
% Function to execute when a selection has been made in the event switcher
% table.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get index of row and column that changed
    iRow = edata.Indices(1);
    iCol = edata.Indices(2);
    
    % check if modified column was VIEW
    if strcmp(src.ColumnName{iCol}, 'VIEW')
        if iRow < size(src.Data, 1)
            % process fvtool (de)activation
            processFVTool(src, iRow);
        else
            % cancel operation if last row was modified
            src.Data{iRow, iCol} = edata.PreviousData;
        end
    else
        % process edited column
        updateFiltEditTable(src, iRow, iCol, edata.PreviousData);
    end
end


% pushbuttonSave_Callback .................................................
function pushbuttonSave_Callback(src,edata)
% Process pushes of the 'Save' button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get window and UITable handles
    hf = src.Parent;
    ht = findobj(hf, 'Type', 'UITable');
    
    % do save
    saveFiltData(ht);
end


% pushbuttonQuit_Callback .................................................
function pushbuttonQuit_Callback(src,edata)
% Process pushes of the 'Quit' button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get window handle
    hf = src.Parent;
    
    % close window
    close(hf)
end


%% fvtool_CloseReq ........................................................
function fvtool_CloseReq(src, edata)
% Close fvtool window and remove checkmark from VIEW column
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ht = src.UserData.ht;
    idx = src.UserData.idx;
    
    % find index of the tool being deleted
    %%% doesn't work; src is not the right type of object
    %idxAllTools = find(~cellfun('isempty',ht.UserData.fvtools));
    %idx = idxAllTools(src == [ht.UserData.fvtools{idxAllTools}]);
    
    % remove checkmark and tool handle
    ht.Data{idx, 6} = false;
    ht.UserData.fvtools{idx} = [];

    % delete tool
    delete(src);
end


%% filtEditUI_CloseReq ....................................................
function filtEditUI_CloseReq(src, edata)
% Close the main window and any associated objects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % check if filtdata has been updated
    ht = findobj(src, 'Type', 'UITable');
    [fd, tableErrMsg] = createFiltDataTable(ht);
    if isempty(tableErrMsg)
        unsavedChanges = ~isequal(fd, sortrows(ht.UserData.filtdata,1));
        if unsavedChanges
            promptStr = sprintf('There are unsaved changes.\nWould you like to save before quitting?');
            opt = questdlg(promptStr, 'Unsaved Changes', 'Yes', 'No', 'Cancel', 'Yes');
            switch opt
                case 'Yes'
                    saved = saveFiltData(ht);
                    if ~saved
                        % cancel if there is a problem with saving
                        return
                    end
                case 'No'
                    % do nothing
                case {'Cancel',''}
                    % cancel quit
                    return
                otherwise
                    error('Unknown option')
            end
        end
    else
        % special processing for case where presence or absence of unsaved
        % changes cannot be determined
        promptStr = sprintf('WARNING: Changes were detected, but they cannot be saved as they are.\n\nQUITTING NOW WILL RESULT IN THE LOSS OF ANY UNSAVED CHANGES.\n\nAre you sure you want to quit anyway?');
        opt = questdlg(promptStr, 'Incomplete Changes', 'Yes', 'No', 'No');
        if any(strcmp(opt, {'No',''}))
            % cancel quit
            return
        end
    end

    % find active fvtools and close them
    fvtools = ht.UserData.fvtools;
    activeFVTools = find(~cellfun('isempty', fvtools))';
    for ii = activeFVTools
        close(fvtools{ii});
    end

    % delete main window
    delete(src);
end