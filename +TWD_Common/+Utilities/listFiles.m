function [filepaths, filenames] = listFiles(dirpath, varargin)
%
% Utility function for retrieving the paths and names of specific files
% within a directory (including other folders, if specified). Nested
% directories can also be searched.
%
% SYNTAX:
%   [filepaths, filenames] = listFiles(dirpath)
%   [__] = listFiles(dirpath, fileTypes)
%   [__] = listFiles(__, Name,Value)
%
% INPUT ARGUMENTS:
%   Required
%   .......................................................................
%   "dirpath" - string specifying full or relative path of directory of 
%       interest
%   .......................................................................
%
%   Optional
%   .......................................................................
%   "fileTypes" - string or cell array of strings specifying the type(s) of
%       files to include. These are usually file extensions such as 'txt'
%       or 'csv' (no preceeding dot is necessary), or empty strings for 
%       files with no extensions. The following special identifiers are
%       also supported:
%
%           '[folders]' - returns folders. This may be specified alone or
%               in addition to other extensions.
%
%           '**' - return all files, but no folders. Only works when
%               specified alone.
%
%           '*' - returns everything, including files and folders. This is 
%               the default. Only works when specified alone.
%   .......................................................................
%
%   Name-Value Pairs
%   .......................................................................
%   "Recursive" - logical specifying whether or not subdirectories should 
%       be searched too. Default is false.
%   .......................................................................
%   "MustContain" - regular expression specifying elements that a path must 
%       have for it to be included
%   .......................................................................
%   "NameMustContain" - like "MustContain", but only applied to file name
%       rather than the full path
%   .......................................................................
%   "MustNotContain" - regular expression specifying elements that a path 
%       should not have for it to be included
%   .......................................................................
%   "NameMustNotContain" - like "MustNotContain", but only applied to file 
%       name rather than the full path
%   .......................................................................
%   
% OUTPUT:
%   .......................................................................
%   "filepaths" - cell array of strings containg the path of each file
%   .......................................................................
%   "filenames" - cell array of strings containg the names of each file
%       with the specified extension(s). Extensions are included in the
%       names too.
%   .......................................................................
%
% DEPENDENCIES:
%   <none>
%
%
% Written by Wilfried Beslin
% Last updated 2023/04/13 using MATLAB R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - The 'export_fig' package is an excellent folder on which to test this
% function, as it contains a great variety of file types, nested folders,
% and some complex file names.

    % input validation
    p = inputParser;
    validExpr = @(arg) validateattributes(arg,{'char'},{'row'});
    defaultLimitExpr = '.*';
    defaultExclusionExpr = '';
    
    % dirpath
    addRequired(p, 'dirpath', @isdir)
    % fileTypes [OPTIONAL]
    defaultTypes = '*';
    validType = @(arg) validateattributes(arg,{'char','cell'},{});
    addOptional(p, 'fileTypes', defaultTypes, validType)
    % Recursive [PARAMETER]
    defaultRecursive = false;
    validRecursive = @(arg) validateattributes(arg,{'logical'},{'scalar'});
    addParameter(p, 'Recursive', defaultRecursive, validRecursive)
    % MustContain [PARAMETER]
    addParameter(p, 'MustContain', defaultLimitExpr, validExpr);
    % NameMustContain [PARAMETER]
    addParameter(p, 'NameMustContain', defaultLimitExpr, validExpr);
    % MustNotContain [PARAMETER]
    addParameter(p, 'MustNotContain', defaultExclusionExpr, validExpr);
    % NameMustNotContain [PARAMETER]
    addParameter(p, 'NameMustNotContain', defaultExclusionExpr, validExpr);
    
    parse(p, dirpath, varargin{:})
    fileTypes = p.Results.fileTypes;
    recursive = p.Results.Recursive;
    pathLimitExpr = p.Results.MustContain;
    nameLimitExpr = p.Results.NameMustContain;
    pathExclusionExpr = p.Results.MustNotContain;
    nameExclusionExpr = p.Results.NameMustNotContain;
    % end input parsing
    
    % process file extension specification
    if ischar(fileTypes)
        if strcmp(fileTypes, '*')
            % all files and folders
            includeFolders = true;
        elseif strcmp(fileTypes, '**')
            % all files but no folders
            includeFolders = false;
        else
            % specific file type; process as 1-element cell array
            fileTypes = {fileTypes};
        end
    end
    
    fileTypeExpr = '.*'; % initialize extension regex to search for absolutely anything
    if iscellstr(fileTypes)
        % check for the [folders] tag
        folderTag = strcmpi(fileTypes, '[folders]');
        includeFolders = any(folderTag);
        
        % process remaining file types, if any
        fileTypes = fileTypes(~folderTag);
        if isempty(fileTypes)
            fileTypeExpr = '';
        else
            % first, determine if files with or without extensions are to
            % be included (it could be one or the other, or both)
            noExt = ismember(fileTypes, '');
            processFilesWithoutExts = any(noExt);
            processFilesWithExts = ~all(noExt);
            
            % create regular expression for finding files without an
            % extension, in case it is needed later. This expression should
            % match any file name that does not contain any dots (except 
            % perhaps at the very beginning). It also should not accept any
            % slashes - this is to be used on file names only, not file
            % paths.
            noExtProcessExpr = '^(\.?[^./\\]+)$';

            % create file match expression
            if processFilesWithExts
                % for file types with an extension, create regular
                % expression for finding all file names containing the
                % appropriate extensions (start by using the extensions in
                % the fileTypes variable, making sure to clean them of any
                % dots that the user may have added). This regex must be
                % able to support file names with additional dots within
                % them, and it must reject any slashes - it is only to be
                % used on file names.
                postDotExpr = '(?<=(^\.?))[^\.].*';
                cleanExtMatch = regexp(fileTypes(~noExt), postDotExpr, 'match');
                cleanExt = [cleanExtMatch{:}];
                extExpr = strcat('\.', cleanExt);

                fileTypeExpr = ['^[^/\\]+(', strjoin(extExpr,'|'), ')$'];
                
                % add the expression for files without extensions, if
                % needed
                if processFilesWithoutExts
                    fileTypeExpr = ['(',fileTypeExpr,')|(',noExtProcessExpr,')'];
                end
                
            elseif processFilesWithoutExts
                % if only files without extensions are to be returned, then
                % set the file matching expression for that case only
                fileTypeExpr = noExtProcessExpr;
            end
        end
    end
    
    % read all files from folder (and subfolders, if specified)
    if recursive
        dirFcnArg = [dirpath, filesep, '**'];
    else
        dirFcnArg = dirpath;
    end
    dirFileInfo = dir(dirFcnArg);
    
    filenamesAll = {dirFileInfo.name}';
    filepathsAll = strcat({dirFileInfo.folder}', filesep, filenamesAll);
    isFolder = [dirFileInfo.isdir]';
    
    % do a preliminary trimming of paths to keep
    if includeFolders
        dotExpr = '^(\.+)$'; % regular expression to find '.' or '..'
        keep1 = ~meetsRegEx(filenamesAll, dotExpr, true);
    else
        keep1 = ~isFolder;
    end
    filepaths = filepathsAll(keep1);
    filenames = filenamesAll(keep1);
    
    % remove stuff that doesnt match specified extensions
    keep2 = isFolder(keep1) | meetsRegEx(filenames, fileTypeExpr, false);
    filepaths = filepaths(keep2);
    filenames = filenames(keep2);
    
    % remove stuff that doesn't meet inclusinon expression, or does meet 
    % exclusion expressions
    keep3 =...
        meetsRegEx(filepaths, pathLimitExpr, true) &...
        meetsRegEx(filenames, nameLimitExpr, true) &...
        ~meetsRegEx(filepaths, pathExclusionExpr, true) &...
        ~meetsRegEx(filenames, nameExclusionExpr, true);
    filepaths = filepaths(keep3);
    filenames = filenames(keep3);
end


%% ------------------------------------------------------------------------
function m = meetsRegEx(str, expr, cs)
    if cs
        regexfcn = @regexp;
    else
        regexfcn = @regexpi;
    end
    exprOut = feval(regexfcn, str, expr);
    m = ~cellfun('isempty',exprOut);
end