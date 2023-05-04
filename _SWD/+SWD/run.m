%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "SWD.run"
%   Written by Wilfried Beslin
%   Last updated May. 4, 2023, using MATLAB R2018b
%
%   Description:
%   Main function for running sperm whale click detector. Takes detection 
%   parameters as input, and runs through all steps to produce MAT files 
%   and "RawEvents" spreadsheets.
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
    doClickDet = ~isempty(dirPath_audio);
    doEventDet = ~isempty(dirName_detResults);
    if ~doClickDet && ~doEventDet
        warning('Nothing to do!')
        return
    end

    % run click detection
    if doClickDet
        SWD.detectSpermClicks(dirPath_root, dirPath_analysis, dirPath_audio, detProtocol, depName, nfft, nMATClicksMax, segDur);
    else
        % check the detection protocol that was used
        detProtocol_mat = TWD_Common.readProtocolTextFile(fullfile(dirPath_analysis,'mat'));
        if ~strcmp(detProtocol, detProtocol_mat)
            % if the protocol that was specified in the Master file is
            % different from the protocol that was previously used to
            % compile clicks, ask the user if they are fine to ignore the
            % detProtocol parameter or if they would rather abort
            prompt = sprintf('The detection protocol specified in the master file was "%s", but the protocol used to compile clicks was "%s". "%s" cannot be used to analyze these clicks. Would you like to continue using "%s" instead? (NOTE: choosing ''No'' will cancel the run)', detProtocol, detProtocol_mat, detProtocol, detProtocol_mat);
            protocol_opt = questdlg(prompt, 'WARNING: Detection Protocol Mismatch', 'Yes', 'No', 'Yes');
            if strcmpi(protocol_opt, 'Yes')
                detProtocol = detProtocol_mat;
            else
                return
            end
        end
    end
    
    % run event detection
    if doEventDet
        SWD.detectSpermEvents(dirPath_root, dirPath_analysis, detProtocol, dirName_detResults, depName);
    end
    
    % finished
    disp('Done')
end