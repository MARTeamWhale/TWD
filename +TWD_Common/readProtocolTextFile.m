%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "readrotocolTextFile"
%   Written by WB
%   Last updated Oct. 13, 2022 using MATLAB R2018b
%
%   Description:
%   reads the detection protocol contained within a 'DetectionProtocol.txt'
%   file..
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function detProtocol = readProtocolTextFile(dirPath)

    detFilePath = fullfile(dirPath,'DetectionProtocol.txt');
    [detFileID,detFileErr] = fopen(detFilePath,'r');
    if detFileID ~= -1
        detProtocol = fscanf(detFileID,'%s');
        fclose(detFileID);
    else
        error('Could not determine detection protocol from file:\n"%s"\n\n%s',detFilePath,detFileErr)
    end
end