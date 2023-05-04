%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "BWD.run"
%   Written by Wilfried Beslin
%   Last updated May. 4, 2023, using MATLAB R2018b
%
%   Description:
%   Main function for running WB's version of the beaked whale click
%   detector. Takes Triton click file paths and detection parameters as
%   input, and runs through all steps to produce MAT files and "RawEvents"
%   spreadsheets.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run(...
    dirPath_root,...
    depName,...
    dirPath_analysis,...
    dirPath_audio,...
    dirName_detResults,...
    detProtocol,...
    segDur,...
    nfft,...
    nMATClicksMax)

    % check which steps to run
    doClickComp = ~isempty(dirPath_audio);
    doBeakedDet = ~isempty(dirName_detResults);
    if ~doClickComp && ~doBeakedDet
        warning('Nothing to do!')
        return
    end

    % run click compilation
    if doClickComp
        BWD.compileClicks(dirPath_root, dirPath_analysis, dirPath_audio, depName, nfft, nMATClicksMax, segDur);
    end
    
    % run event detection
    if doBeakedDet
        BWD.detectBeaked(dirPath_root, depName, dirPath_analysis, detProtocol, dirName_detResults);
    end
    
    % finished
    disp('Done')
end