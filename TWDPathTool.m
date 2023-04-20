%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "TWDPathTool"
%   Written by WB
%   Last updated Sep. 27, 2019, using MATLAB R2018b
%
%   Description:
%   Adds or removes a specified TWD/BWD/SWD folder to the MATLAB path. If
%   added, the folder paths will be moved to the top.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% - This is very simple for now, ideally I'd like it to have features like:
%   -- prompt user to select a TWD version, starting in m-file folder
%       --- ensure user selects a TWD folder, otherwise throw error
%   -- prompt user to select modules to include (if multiple versions)
%       --- perhaps with UITable Or listbox 
%   -- maybe check if:
%       --- current version is overshadowed
%       --- current version is already at top 
%       --- current version is not yet at top
%   -- option to save changes permanently

function TWDPathTool()

    % determine path to current TWD directory
    dirPath_root = mfilename('fullpath');
    [dirPath_root,fcnName,~] = fileparts(dirPath_root);
    [~,versionName,~] = fileparts(dirPath_root);
    
    % prompt user for action
    noteStr1 = '- Changes will persist only for the current MATLAB session unless saved';
    noteStr2 = '- If multiple versions of BWD or SWD exist within this folder, they will each be added to the path in alphabetical order';
    promptStr = sprintf('%s\n%s\nHow would you like to edit the MATLAB path?\n\nNOTES:\n%s\n%s',versionName,dirPath_root,noteStr1,noteStr2);
    opt = questdlg(promptStr,fcnName,'Add to Top','Remove from Path','Add to Top');
    
    if isempty(opt)
        return
    else
        % get all subfolders within this TWD
        disp('Retrieving TWD folder paths...')
        TWDPaths = genpath(dirPath_root);
        TWDPaths = strsplit(TWDPaths,';')';
        
        % process user option
        switch opt
            case 'Add to Top'
                disp('Adding folders to top of search path...')
                addpath(TWDPaths{:},'-begin');
            case 'Remove from Path'
                disp('Removing folders from search path...')
                rmpath(TWDPaths{:});
            otherwise
                error('Unrecognized option')
        end
        disp('Done')
    end
end