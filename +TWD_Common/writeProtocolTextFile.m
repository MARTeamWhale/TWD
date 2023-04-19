%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "writeProtocolTextFile"
%   Written by WB
%   Last updated Oct. 13, 2022 using MATLAB R2018b
%
%   Description:
%   Creates a simple text file denoting the detection protocol that was
%   used during a TWD run.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeProtocolTextFile(dirPath_out,detProtocol)

    txtFileName = 'DetectionProtocol.txt';
    txtFilePath = fullfile(dirPath_out,txtFileName);
    [txtFileID,txtFileError] = fopen(txtFilePath,'w');
    if txtFileID ~= -1
        fprintf(txtFileID,'%s',detProtocol);
        fclose(txtFileID);
    else
        warning('Could not write detection protocol file (%s):\n%s',txtFileName,txtFileError)
    end
end