%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "readTritonFile"
%   Written by WB
%   Last updated Dec. 14, 2018, using MATLAB R2016b
%
%   Reads in a Triton .c or .cTg file and returns the results as a cell
%   array.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = readTritonFile(fPath)

    fID = fopen(fPath);
    data = textscan(fID,'%f %f %s');
    fclose(fID);
end