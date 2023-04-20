%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "cleanMetadata" 
%   Written by Wilfried Beslin
%   Last Updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   Removes unwanted Triton detector output files from "metadata" folders, 
%   such as gTg, us, w, and ccc. Note that this code will only search for
%   files within folders called "metadata".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUT - CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full path to metadata folder
dirPath = 'D:\BWD_Results\metadata';

% decide if subfolders should also be cleaned (all levels)
includeSub = true;

% list of unwanted file types (case insensitive)
%%% default is {'gtg','us','w','ccc'}
extRemove = {'gtg','us','w','ccc'}; 

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


do_cleanMetadata(dirPath, includeSub, extRemove)


%%  -----------------------------------------------------------------------
function do_cleanMetadata(dirPath, includeSub, extRemove)

    % prompt user to confirm operation
    if includeSub
        fStr = 'folder and all its subfolders';
    else
        fStr = 'folder';
    end
    promptStr = sprintf('%s\n\nAre you sure you want to delete the following file types from this %s?\n\n%s',dirPath,fStr,strjoin(extRemove,', '));
    opt = questdlg(promptStr,'Warning','Yes','No','Yes');
    
    % clean folder(s)
    if strcmp(opt, 'Yes')
        filesToDelete = TWD_Common.Utilities.listFiles(dirPath, extRemove, includeSub, 'MustContain','(metadata[/\\][^/\\]+)$');
        if ~isempty(filesToDelete)
            fprint('Removing the following files:\n    %s', strjoin(filesToDelete, '\n    '))
            delete(filesToDelete)
            disp('done')
        else
            disp('Found nothing to delete.')
        end
    end
end