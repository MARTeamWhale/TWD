%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "readDateTime"
%   Written by WB
%   Last updated Dec. 3, 2019, using MATLAB R2018b
%
%   Returns a datetime object corresponding to the timestamp in a file 
%   name (or list of file names). Full paths are also accepted.
%   IMPORTANT: this only works if the date/time is encoded somewhere in the 
%   filename as "yyyyMMdd.HHmmss", where "." is any character that is not a
%   digit or a slash.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dt = readDateTime(fileName)
    
    % check if input is a string or cellstr and process accordingly
    if ischar(fileName)
        % extract datetime from filename
        dt = doDTRead(fileName);
    elseif iscellstr(fileName)
        % extract datetime from each filename
        dt = NaT(size(fileName)); % initialize datetime array
        for ii = 1:numel(dt)
            fileNameii = fileName{ii};
            dt(ii) = doDTRead(fileNameii);
        end
    else
        % input argument not recognized
        error('Bad input');
    end
end

%% doDTRead ---------------------------------------------------------------
function dt = doDTRead(fileName)
% Extract datetime from char string
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % define regex for reading datetime
    expr = '\d{8}[^\d\\/]\d{6}'; % expression for finding strings of the form "DDDDDDDD.DDDDDD", where "." is any character that is not a digit or a slash
    iSep = 9; % index of separating character

    % extract datetime info
    dtStr = regexp(fileName,expr,'match');
    try
        dtStr = dtStr{end}; % "end" ensures that only the match for filename is used, not other potential matches in folder path
        dtStr = dtStr([1:(iSep-1),(iSep+1):end]);
        dt = datetime(dtStr,'InputFormat','yyyyMMddHHmmss');
    catch
        warning('Could not read timestamp for file "%s"',fileName)
        dt = NaT;
    end
end