%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "cleanMetadata" 
%   Written by Wilfried Beslin
%   Last Updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   Removes unwanted Triton detector output files, such as gTg, us, w, and
%   ccc. This script can also clean multiple folders in a root path at
%   once, but in this case, it will only search for folders called
%   "metadata".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUT - CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full path to the folder or root path to be cleaned
dirPath = 'D:\TWD_analysis\DEP_yyyy_mm\metadata';

% decide if subfolders should be cleaned. 
% If "true", this option will only look for folders called "metadata".
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
        fStr = 'all "metadata" folders in this path';
    else
        fStr = 'this folder';
    end
    promptStr = sprintf('%s\n\nAre you sure you want to delete the following file types from %s?\n\n%s', dirPath, fStr, strjoin(extRemove,', '));
    opt = questdlg(promptStr,'Warning','Yes','No','Yes');
    
    % clean folder(s)
    if strcmp(opt, 'Yes')
        if includeSub
            filesToDelete = TWD_Common.Utilities.listFiles(dirPath, extRemove, 'Recursive',true, 'MustContain','([/\\]metadata[/\\][^/\\]+)$');
        else
            filesToDelete = TWD_Common.Utilities.listFiles(dirPath, extRemove, 'Recursive',false);
        end
        numFiles = numel(filesToDelete);
        if numFiles > 0
            for ii = 1:numFiles
                file_ii = filesToDelete{ii};
                fprintf('Removing file "%s"\n', file_ii)
                try
                    delete(file_ii)
                catch ME
                    warning('Failed to remove file:\n%s', ME.message)
                end
            end
            disp('done')
        else
            disp('Found nothing to delete.')
        end
    end
end