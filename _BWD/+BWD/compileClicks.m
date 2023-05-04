%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "compileClicks"
%   Written by WB, based on "analyze_HARP_AMAR250" by SBP and JS
%   Last updated May. 4, 2023, using MATLAB R2018b
%
%   Description:
%   Extracts click data from WAV files based on Triton detections (cTg 
%   files), and compiles them into MAT files. At least one MAT per WAV,
%   depending on how many clicks were detected.
%   If the process is interrupted, the function will automatically resume
%   where it left off by searching for existing MAT files.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

function compileClicks(dirPath_root, dirPath_analysis, dirPath_audio, depName, nfft, nClicksMax, segDur)
    
    % define parameters for click extraction and spectral analysis
    clickOffset_sec = [0.001, 0.004]; % length of buffer before and after click, in seconds (previously hardcoded to 250 and 999 samples)
    noiseOffset_sec = 0.005; % offset before start of click from which noise samples are computed, in seconds (previously hardcoded to 1250 samples)  
    
    % get/set metadata and output MAT folder paths
    dirPath_metadata = fullfile(dirPath_analysis,'metadata');
    dirPath_out = fullfile(dirPath_analysis,'mat');
    
    % get cTg and wav file names
    [filePaths_ctg, fileNames_ctg] = TWD_Common.Utilities.listFiles(dirPath_metadata, 'cTg', 'Recursive',true);
    [filePaths_wav, fileNames_wav] = TWD_Common.Utilities.listFiles(dirPath_audio, 'wav', 'Recursive',true);
    nFiles = numel(fileNames_wav);
    
    % load filter cutoff frequencies
    load(fullfile(dirPath_root,'BandpassFilterParams.mat'), 'filtdata')
    
    % define output directory and check if there is existing output
    if isfolder(dirPath_out)
        % check for existing MAT files, get timestamp of the last one, and 
        % compare that with WAV file names to determine where to resume
        [~, matNames] = TWD_Common.Utilities.listFiles(dirPath_out, 'mat');
        if ~isempty(matNames)
            lastTime = TWD_Common.Utilities.readDateTime(matNames{end});
            lastTime.Format = 'yyyyMMdd_HHmmss';
            timeStrParts = strsplit(char(lastTime),'_');
            leftOff = contains(fileNames_wav,timeStrParts{1}) & contains(fileNames_wav,timeStrParts{2});
            iStart = find(leftOff) + 1;
        else
            iStart = 1;
        end
        if iStart <= nFiles
            fprintf('Resuming at recording #%d/%d ("%s")\n',iStart,nFiles,fileNames_wav{iStart})
        else
            return
        end
    else
        % create new output folder and start from beginning
        mkdir(dirPath_out)
        iStart = 1;
        fprintf('Starting new extractions for deployment "%s"\n',depName) 
    end
    
    % loop through each file and compute click data
    for ii = iStart:nFiles
        
        % make sure timestamps of WAV and cTg match
         DTii_wav = TWD_Common.Utilities.readDateTime(fileNames_wav{ii});
         DTii_ctg = TWD_Common.Utilities.readDateTime(fileNames_ctg{ii});
         assert(DTii_wav == DTii_ctg,'WAV and Triton cTg files do not match')
        
        % get file paths
        filePathii_wav = filePaths_wav{ii};
        filePathii_ctg = filePaths_ctg{ii};
        
        % read cTg file
        ctgDataii = BWD.readTritonFile(filePathii_ctg);
        clickStartEndTimesii = [ctgDataii{1},ctgDataii{2}];
        nClicksii = size(clickStartEndTimesii,1);
        
        if nClicksii > 0
            % split by number of clicks
            iClickPartStartsii = (1:nClicksMax:nClicksii)';
            iClickPartEndsii = unique([nClicksMax:nClicksMax:nClicksii,nClicksii])';
            iClickPartRangesii = [iClickPartStartsii,iClickPartEndsii];
            nPartsii = size(iClickPartRangesii,1);

            % loop through each part and compile clicks 
            for jj = 1:nPartsii
                iClicksjj = iClickPartRangesii(jj,1):iClickPartRangesii(jj,2);
                clickStartEndTimesjj = clickStartEndTimesii(iClicksjj,:);
                nClicksjj = iClicksjj(end) - iClicksjj(1) + 1;
                fprintf('Recording %d/%d, part %d/%d: extracting %d clicks\n',ii,nFiles,jj,nPartsii,nClicksjj)
                
                [clickDatajj,recStartjj] = computeClickData(filePathii_wav, segDur, clickStartEndTimesjj, clickOffset_sec, noiseOffset_sec, nfft, filtdata);
                
                % patch for very rare case where all clicks detected by
                % Triton were out-of-bounds, resulting in empty data
                % struct which causes problems later. Rare, but I've had it 
                % happen
                if size(clickDatajj.yFilt,1) == 0
                    disp('     All clicks out of bounds!')
                    continue
                end

                % save click data
                recStartjj.Format = 'yyyyMMdd_HHmmss';
                outFileNameii = sprintf('%s_%s_%d',depName,char(recStartjj),jj);
                outPathii = fullfile(dirPath_out,outFileNameii);
                save(outPathii,'-struct','clickDatajj')
            end
        else
            fprintf('Recording %d/%d: no clicks detected\n',ii,nFiles)
        end
    end
end

%% computeClickData -------------------------------------------------------
function [clickData,recStartDT] = computeClickData(filePath_wav, segDur, clickStartEndTimes, clickTimeOffset, noiseTimeOffset, nfft, filtdata)
% Returns struct of click parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % read audio file info
    wavInfo = audioinfo(filePath_wav);
    Fs = wavInfo.SampleRate;
    totSamples = wavInfo.TotalSamples; % using "totSamples" to avoid confusion with the legacy parameter "nSamples"
    recDur = wavInfo.Duration;
    recDurDR = seconds(recDur);
    
    % read datetime and split recording into segments
    recStartDT = TWD_Common.Utilities.readDateTime(filePath_wav);
    recEndDT = recStartDT + recDurDR;
    segStartsDT = (recStartDT:segDur:recEndDT)';
    if segStartsDT(end) == recEndDT
        segStartsDT(end) = [];
    end
    segStartMat = datevec(segStartsDT);
    segDursDR = diff([segStartsDT;recEndDT]);
    segDurs = seconds(segDursDR);

    % convert time-based variables to samples
    clickOffset_samples = round(clickTimeOffset*Fs);
    noiseOffset_samples = round(noiseTimeOffset*Fs);
    
    % initialize click and noise position structs
    clickPos = struct(...
        'StartTimes',[],...
        'EndTimes',[],...
        'StartSamples',[],...
        'EndSamples',[],...
        'StartSamplesExtended',[],...
        'EndSamplesExtended',[]);
    noisePos = struct(...
        'StartSamples',[],...
        'EndSamples',[]);

    % get number of clicks and store start/end times
    nClicks = size(clickStartEndTimes,1);
    clickPos.StartTimes = clickStartEndTimes(:,1);
    clickPos.EndTimes = clickStartEndTimes(:,2);

    % begin click extraction
        
    % get start/end samples for each click
    clickPos.StartSamples = round(clickPos.StartTimes*Fs);
    clickPos.EndSamples = round(clickPos.EndTimes*Fs);
    clickPos.StartSamplesExtended = clickPos.StartSamples - clickOffset_samples(1) - 1;
    clickPos.EndSamplesExtended = clickPos.StartSamplesExtended + clickOffset_samples(2) - 1;
    % the "extended" samples result in a fixed duration of length 
    % clickOffset(2); somewhat dangerous because it ignores the 
    % TK-based duration, which could in theory be longer. Granted, this 
    % situation is unlikely, but it would be nice to have a more robust 
    % method.

    % get range of noise samples for each click
    noisePos.StartSamples = clickPos.StartSamples - noiseOffset_samples;
    noisePos.EndSamples = clickPos.StartSamples;
    % I think this should all be shifted by -1, but that's not how
    % it was in the original code.

    % remove clicks that are too close to the start or end of
    % the timeseries (including noise samples)
    isOut = ...
        noisePos.StartSamples < 1 ...
        | clickPos.EndSamplesExtended > totSamples;
    if any(isOut)
        % edit nClicks, clickPos and noisePos.
        % For clickPos and noisePos, use a loop to avoid having to type
        % out every field. Slightly less efficient, but less
        % error-prone if the code changes.
        nClicks = sum(~isOut);
        posFields = unique([fieldnames(clickPos);fieldnames(noisePos)]);
        for ii = 1:numel(posFields)
            fieldii = posFields{ii};
            if isfield(clickPos,fieldii)
                clickPos.(fieldii) = clickPos.(fieldii)(~isOut);
            end
            if isfield(noisePos,fieldii)
                noisePos.(fieldii) = noisePos.(fieldii)(~isOut);
            end
        end
    end

    % extract and filter waveforms of each click (and corresponding noise)
    nClickSamplesExtended = unique(clickPos.EndSamplesExtended - clickPos.StartSamplesExtended) + 1; % should be scalar
    nNoiseSamples = unique(noisePos.EndSamples - noisePos.StartSamples) + 1; % should be scalar
    xClickFilt = NaN(nClicks,nClickSamplesExtended);
    xNoiseFilt = NaN(nClicks,nNoiseSamples);
    for ii = 1:nClicks
        [xClickFiltii, xNoiseFiltii] = BWD.extractFilterClick(filePath_wav,...
            [clickPos.StartSamplesExtended(ii),clickPos.EndSamplesExtended(ii)],...
            [noisePos.StartSamples(ii),noisePos.EndSamples(ii)],...
            filtdata);
        xClickFilt(ii,:) = xClickFiltii;
        xNoiseFilt(ii,:) = xNoiseFiltii;
    end

    % compute click duration (ms)
    clickDur_ms = (clickPos.EndTimes - clickPos.StartTimes)*1000;

    % compute ICI (ms)
    ICI_ms = (clickPos.StartTimes(2:end) - clickPos.StartTimes(1:(end-1)))*1000;

    % compute click features
    [specClick,specNoise,peakFr,F0,bw3db,bw10db,ppSignal,rmsSignal,rmsNoise,snr,slope,nSamples] =...
        TWD_Common.computeParameters(xClickFilt,xNoiseFilt,clickDur_ms,Fs,nfft,clickOffset_samples(1));
    
    % save all features to struct.
    % Designed to be consistent with SBP's code.
    clickData = struct(...
        'bw10db',bw10db,...
        'bw3db',bw3db,...
        'dur',clickDur_ms,...
        'F0',F0,...
        'fs',Fs,...
        'ici',ICI_ms,...
        'nSamples',nSamples,...
        'offset',clickOffset_samples(1),...
        'peakFr',peakFr,...
        'pos',[clickPos.StartTimes,clickPos.EndTimes],...
        'ppSignal',ppSignal,...
        'rawDur',segDurs,...
        'rawStart',segStartMat,...
        'rmsNoise',rmsNoise,...
        'rmsSignal',rmsSignal,...
        'slope',slope,...
        'snr',snr,...
        'specClick',specClick,...
        'specNoise',specNoise,...
        'yFilt',xClickFilt,...
        'yNFilt',xNoiseFilt);
end