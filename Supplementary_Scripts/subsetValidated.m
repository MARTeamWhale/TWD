%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "subsetValidated"
%   Written by Wilfried Beslin
%   Last updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   Isolates a subset of events in a "validated" spreadsheet based on
%   species ID code. The output is saved as a new spreadsheet. Each output
%   sheet will have a unique filename, so it's possible to create new
%   subsets without overwriting old ones.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full path to validation spreadsheet to subset
filePath = 'D:\TWD_analysis\DEP_yyyy_mm\results01\DEP_Target_Validated.xlsx';

% species codes to include
%%% - must be a column vector, use ";" to separate multiple digits
%%% - note that any occurence of these digits within an ID code will be
%%%     considered a match. So for example, if "0" is specified, any code
%%%     with the digit zero somewhere (like 20) will be included.
codeNum = [0;8];

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


valSubset = do_subsetValidated(filePath, codeNum);


%% ------------------------------------------------------------------------
function outData = do_subsetValidated(filePath, codeNum)
    % load spreadsheet
    inData = readtable(filePath);

    % convert codes to strings
    codeStr = cellstr(num2str(codeNum));
    SIDStr = cellstr(num2str(inData.Species));

    % get rows containing species ID codes of interest
    keep = contains(SIDStr,codeStr);

    % create subset of table
    outData = inData(keep,:);

    % save output to file (get unique filename to avoid overwriting
    % existing files)
    [dirPath, inFileName, fileExt] = fileparts(filePath);
    outFileSuffix = '_Sub';
    outFileNum = 1;
    fileSaved = false;
    while ~fileSaved
        outFileName = [inFileName,outFileSuffix,num2str(outFileNum)];
        outFilePath = fullfile(dirPath,[outFileName,fileExt]);
        if logical(exist(outFilePath,'file'))
            outFileNum = outFileNum + 1;
        else
            writetable(outData,outFilePath)
            fprintf('Saved subset validations table to:\n"%s"\n', outFilePath)
            fileSaved = true;
        end
    end
end