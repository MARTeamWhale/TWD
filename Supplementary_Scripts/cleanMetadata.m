%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "cleanMetadata" 
%   Written by Wilfried Beslin
%   Last Updated Nov. 7, 2019, using MATLAB R2018b
%
%   Description:
%   Removes unwanted Triton detector output files, leaving only .c and .cTg
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function cleanMetadata()

    % INPUT - CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% full path to metadata folder
    dirPath = 'C:\Users\BeslinW\Documents\BeakedWhales_Analyses\BWD_TEST_subdirs\metadata';
    %%% decide if subfolders should also be cleaned (all levels)
    includeSub = true;
    %%% list of unwanted file types (case insensitive)
    extRemove = {'gtg','us','w','ccc'};
    % END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % prompt user to confirm operation
    if includeSub
        fStr = 'folder and all its subfolders';
    else
        fStr = 'folder';
    end
    promptStr = sprintf('%s\n\nAre you sure you want to delete the following file types from this %s?\n\n%s',dirPath,fStr,strjoin(extRemove,', '));
    opt = questdlg(promptStr,'Warning','Yes','No','Yes');
    
    % clean folder(s)
    if strcmp(opt,'Yes')
        doClean(dirPath,extRemove,includeSub);
        disp('done')
    end
end

%% ------------------------------------------------------------------------
function doClean(dirPath,extRemove,includeSub)
    fprintf('Cleaning %s\n',dirPath);
    nExt = numel(extRemove);
    for ii = 1:nExt
        extii = extRemove{ii};
        delStr = [dirPath,filesep,'*.',extii];
        delete(delStr)
    end
    
    % check subfolders
    if includeSub
        % get subdirectories
        dirInfoFull = dir(dirPath);
        subdirs = {dirInfoFull.name}';
        subdirs = subdirs([dirInfoFull.isdir]);
        dotExpr = '^(\.+)$'; % regular expression to find a just series of dots
        dotRegExOut = regexp(subdirs,dotExpr,'match');
        notDots = cellfun('isempty',dotRegExOut);
        subdirs = subdirs(notDots);
        subdirPaths = fullfile(dirPath,subdirs);

        % clean files in subdirectories
        for ii = 1:numel(subdirPaths)
            doClean(subdirPaths{ii},extRemove,includeSub);
        end 
    end
end