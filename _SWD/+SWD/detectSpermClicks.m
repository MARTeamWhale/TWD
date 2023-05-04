%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "detectSpermClicks"
%   Written by Wilfried Beslin
%   Last updated May. 4, 2023, using MATLAB R2018b
%
%   Description:
%   Runs basic sperm whale click detection (CABLE detector) followed by 
%   click discrimination. Results are saved in MAT files. 
%   At least one MAT per WAV, depending on how many clicks were detected.
%   If the process is interrupted, the function will automatically resume
%   where it left off by searching for existing MAT files.
%
%   2019/09/04 update: allows setting onThresh as optional input param
%   2019/09/26 update: allows setting offThresh as optional input param too
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES 
%{
% 2019/10/07
%   - MAT file contents has been edited to be identical to BWD, sacrificing
% efficiency for compatibility

%!!! Be careful with segmentation vs. number of clicks.
% With Triton, clicks were detected from whole files, and MATs were split
% by nClicks only.
% CABLE is less efficient with long files and should segment before
% detection - so essentially, the WAV files need to be split.
% In other words, using CABLE requires that segmentation be applied during 
% MAT file compilation, unlike BWD where segmentation was only considered 
% during event detection. This means that unlike BWD, SWD MAT file 
% timestamps may not always match with WAV file timestamps. The nMaxClicks 
% factor can still be applied for SWD though.

% Maybe add waitbar too
%}

function detectSpermClicks(dirPath_root, dirPath_analysis, dirPath_audio, detProtocol, depName, nfft, nClicksMax, segDur, varargin)

    % get path to detection criteria directory
    dirPath_SpermDetCrit = fullfile(dirPath_root,'DetectionCriteria',detProtocol);
    
    % get list of WAV files
    [wavFilePaths, wavFileNames] = TWD_Common.Utilities.listFiles(dirPath_audio, 'wav', 'Recursive',true);
    nFiles = numel(wavFileNames);
    
    % peek at the first audio file to get Fs
    wavInfo1 = audioinfo(wavFilePaths{1});
    FsAMAR = wavInfo1.SampleRate;
    
    % load filter cutoff frequencies
    load(fullfile(dirPath_root,'BandpassFilterParams.mat'), 'filtdata');
    
    % determine standard filter cutoff frequencies based on dataset sampling rate
    iFiltCutoff = filtdata.SamplingRate == FsAMAR;
    switch sum(iFiltCutoff)
        case 1
            Fc1 = filtdata.CutoffFreq1(iFiltCutoff);
            Fc2 = filtdata.CutoffFreq2(iFiltCutoff);
            ord = filtdata.Order(iFiltCutoff);
        case 0
            error('No bandpass filter cutoff frequencies have been specified for Fs=%d Hz. They must be added using "editFilterParams" before click compilation can be run on this dataset.', Fs)
        otherwise
            error('More than one bandpass filter cutoff frequency specification exists for Fs=%d Hz; cannot determine which one to use.', Fs)
    end
    %{
    switch FsAMAR
        case {250000, 256000}
            Fc1 = 2000;
            Fc2 = 12000;
        otherwise
            error('Unsupported sampling rate')
    end
    %}
    
    % process onThresh and offThresh
    if numel(varargin) == 1
        % varargin{1} expected to contain [onThresh, offThresh]
        threshOn = varargin{1}(1);
        threshOff = varargin{1}(2);
    elseif numel(varargin) == 0
        threshOn = 20; % CABLE value = 10
        threshOff = 2; % CABLE value = 1
    else
        error('Too many input arguments')
    end

    % Define CABLE detection parameters
    % These are CABLE defaults, except for threshOn, threshOff, and IPIRange(1)
    clickDetParams = struct(...
        'threshOn',threshOn,...
        'threshOff',threshOff,...
        'alphaSignal',0.2,...
        'alphaNoiseOn',0.000002,...
        'alphaNoiseOff',0.0002,...
        'minEchoProp',0.6,...
        'IPIRange',[0,9],... % milliseconds; CABLE lower limit was 2
        'maxClickDuration',40,... % milliseconds
        'minClickSep',0,... % seconds
        'minPulseDuration',0.05); % milliseconds
    FsDetect = 48000; % As with CABLE, will downsample to this; faster to do it at this rate
    channelIndex = 1;

    % create digital filter objects
    dFilt = struct();
    %%% click detection; note that this one deviates from the CABLE filter
    dFilt.Detect = designfilt('highpassiir',...
        'DesignMethod','butter',...
        'SampleRate',FsDetect,...
        'FilterOrder',10,...
        'HalfPowerFrequency',2000);
    %%% "standard" wideband filter for discrimination 
    dFilt.Standard = designfilt('bandpassiir',...
        'DesignMethod','butter',...
        'SampleRate',FsAMAR,...
        'FilterOrder',ord,...
        'HalfPowerFrequency1',Fc1,...
        'HalfPowerFrequency2',Fc2);
    %%% Narrowband filter to help see clicks during validation; similar to
    %%% click detection filter, but for native sampling rate
    %dFilt.View = designfilt('bandpassiir',...
    %    'DesignMethod','butter',...
    %    'SampleRate',FsAMAR,...
    %    'FilterOrder',10,...
    %    'HalfPowerFrequency1',2000,...
    %    'HalfPowerFrequency2',24000);
    
    % create ClickDiscriminator object
    spermCritFilePath = fullfile(dirPath_SpermDetCrit,'ClickDiscrimParams.xlsx');
    discrimObj = TWD_Common.ClickDiscriminator(spermCritFilePath);
    
    % create mat directory if it doesn't exist yet. If it does, inspect
    % contents to see where the routine left off.
    dirPath_mat = fullfile(dirPath_analysis,'mat');
    if isfolder(dirPath_mat)
        % check existing MAT files
        iStart = checkExistingFiles(wavFileNames,dirPath_mat);
    else
        % create new folder and start from scratch
        mkdir(dirPath_mat)
        iStart = 1;
        fprintf('Starting new SW click detections for deployment "%s"\n',depName) 
    end
    
    % loop through each file
    for ii = iStart:nFiles
        wavPathii = wavFilePaths{ii};
        wavInfoii = audioinfo(wavPathii);
        
        % verify sampling rate
        assert(wavInfoii.SampleRate == FsAMAR, 'Encountered audio file with unexpected sampling rate. All files in a dataset must have the same rate.');
        
        % extract segment start times and durations
        [segStartMatii,segDursii,segSamplesii] = getSegTimes(wavInfoii,segDur);
        nSegii = numel(segDursii);
        
        % loop through each segment
        for jj = 1:nSegii
            segStartTimejj = segStartMatii(jj,:);
            segDurjj = segDursii(jj);
            segSamplesjj = segSamplesii(jj,:);
            fprintf('Recording %d/%d, part %d/%d: ',ii,nFiles,jj,nSegii)
            
            % detect clicks
            [clickTimesjj,detClickPeaksjj,detClickNoisejj] = detectClicks(wavPathii,segSamplesjj,channelIndex,FsDetect,dFilt.Detect,clickDetParams);
            
            % compile information on clicks
            clickDatajj = extractClicks(wavPathii,segSamplesjj,channelIndex,clickTimesjj,dFilt.Standard,nfft);
            %%% add CABLE detection info
            %clickDatajj.detClickPeaks = detClickPeaksjj';
            %clickDatajj.detClickNoise = detClickNoisejj';
            
            % run click discrimination
            isSpermjj = discrimObj.evalTargetClicks(clickDatajj);
            fprintf('%d sperm whale clicks detected\n',sum(isSpermjj))
            
            % create MAT file(s)
            compileMAT(dirPath_mat,depName,segStartTimejj,segDurjj,nClicksMax,clickDatajj,isSpermjj,wavInfoii.SampleRate); 
            
            % DEBUG
            %{
            xDB = audioread(wavPathii,segSamplesjj);
            xDB = filtfilt(dFilt.Standard,xDB);
            x2DB = filtfilt(dFilt.View,xDB);
            tDB = (1:numel(xDB))/wavInfoii.SampleRate;
            figure
            plot(tDB,xDB,':')
            hold on
            plot(tDB,x2DB)
            yMaxDB = max(abs(xDB))*1.05;
            patchX1DB = [clickDatajj.pos(:,1),clickDatajj.pos(:,1),clickDatajj.pos(:,2),clickDatajj.pos(:,2)]';
            patchX2DB = patchX1DB(:,isSpermjj);
            patchY1DB = repmat([-yMaxDB;yMaxDB;yMaxDB;-yMaxDB],1,numel(isSpermjj));
            patchY2DB = repmat([-yMaxDB;yMaxDB;yMaxDB;-yMaxDB],1,sum(isSpermjj));
            patchCol1DB = [1,1,0.75];
            patchEdgeCol1DB = [0.9,0.9,0.6];
            patchCol2DB = [1,0.75,0.75];
            patchEdgeCol2DB = [0.9,0.6,0.6];
            h_clicks1DB = patch(patchX1DB,patchY1DB,patchCol1DB,'EdgeColor',patchEdgeCol1DB);
            h_clicks2DB = patch(patchX2DB,patchY2DB,patchCol2DB,'EdgeColor',patchEdgeCol2DB);
            uistack(h_clicks2DB,'bottom')
            uistack(h_clicks1DB,'bottom')
            keyboard
            %}
        end
    end
    
    % add a plain text file marking the name of the detection protocol
    TWD_Common.writeProtocolTextFile(dirPath_mat, detProtocol)
end

%% checkExistingFiles -----------------------------------------------------
function iStart = checkExistingFiles(wavFileNames,dirPath_mat)
% checks MAT directory for existing files and determines where to resume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: need to account for possibility that recording was partially
% processed, where it was split into multiple MATs but only a subset of
% MATs were saved to disk. Easiest way to deal with this is to redo last
% file.
    nFiles = numel(wavFileNames);
    [~, matNames] = TWD_Common.Utilities.listFiles(dirPath_mat, 'mat');
    if ~isempty(matNames)
        % MAT files exist - check where left off
        lastTime = TWD_Common.Utilities.readDateTime(matNames{end});
        lastTime.Format = 'yyyyMMdd_HHmmss';
        timeStrParts = strsplit(char(lastTime),'_');
        leftOff = contains(wavFileNames,timeStrParts{1}) & contains(wavFileNames,timeStrParts{2});
        iStart = find(leftOff) + 1;
    else
        % no MAT files yet - start at file 1
        iStart = 1;
    end
    if iStart <= nFiles
        fprintf('Resuming at recording #%d/%d ("%s")\n',iStart,nFiles,wavFileNames{iStart})
    else
        disp('All recordings already processed')
    end
end

%% getSegTimes ------------------------------------------------------------
function [segStartMat,segDurs,segSamples] = getSegTimes(wavInfo,segDur)
% Reads timestamp of a WAV file and returns the start times and durations 
% of its segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get recording duration (duration object)
    recDurDR = seconds(wavInfo.Duration);

    % get recording start and end times (datetime objects
    recStartDT = TWD_Common.Utilities.readDateTime(wavInfo.Filename);
    recEndDT = recStartDT + recDurDR;
    
    % get segment start times (datetime objects)
    segStartsDT = (recStartDT:segDur:recEndDT)';
    if segStartsDT(end) == recEndDT
        segStartsDT(end) = [];
    end
    
    % get segment durations (duration objects)
    segDursDR = diff([segStartsDT;recEndDT]);
    
    % convert datetimes/durations to double
    segStartMat = datevec(segStartsDT);
    segDurs = seconds(segDursDR);
    
    % translate to range of samples for each segment
    segDurs_samples = segDurs*wavInfo.SampleRate;
    segEndSamples = [cumsum(segDurs_samples(1:end-1));wavInfo.TotalSamples];
    segStartSamples = [1;cumsum(segDurs_samples(1:end-1))+1];
    segSamples = [segStartSamples,segEndSamples];
end

%% detectClicks -----------------------------------------------------------
function [clickTimes,clickPeaks,clickNoise] = detectClicks(wavPath,sampleIndex,channelIndex,FsDetect,dFiltDetect,clickDetParams)
% detects clicks using CABLE detector for one file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % load and prepare recording samples for detection
    [x,tCutStart,tCutEnd] = importAudio_Detection(wavPath,sampleIndex,channelIndex,FsDetect,dFiltDetect);
    
    % compute envelope
    xEnv = abs(hilbert(x));
    
    % run CABLE click detector
    %%% Note that output is a 2-by-N matrix!
    [clickSamples,pageNoise] = SWD.CABLE.detectClicks(...
        xEnv,...
        FsDetect,...
        clickDetParams.threshOn,...
        clickDetParams.threshOff,...
        clickDetParams.alphaSignal,...
        clickDetParams.alphaNoiseOn,...
        clickDetParams.alphaNoiseOff,...
        clickDetParams.minEchoProp,...
        clickDetParams.IPIRange,...
        clickDetParams.maxClickDuration,...
        clickDetParams.minClickSep,...
        clickDetParams.minPulseDuration);
    
    % remove 1st and last detections as they are usually artefacts in this
    % case
    clickSamples = clickSamples(:,2:(end-1));
    
    % get click peaks and noise values
    nClicks = size(clickSamples,2);
    clickPeaks = zeros(1,nClicks);
    clickNoise = zeros(1,nClicks);
    for ii = 1:nClicks
        iStartii = clickSamples(1,ii);
        iEndii = clickSamples(2,ii);
        
        clickPowii = xEnv(iStartii:iEndii).^2;
        noiseii = pageNoise(iStartii:iEndii);
        clickPeaks(ii) = max(clickPowii);
        clickNoise(ii) = mean(noiseii);
    end
    
    % convert samples to times
    clickTimes = clickSamples/FsDetect;
    clickTimes = clickTimes + tCutStart;
    
    %%% DEBUG
    %{
    t = (1:numel(x))/FsDetect + tCutStart;
    figure
    plot(t,x)
    hold on
    plot(t,xEnv)
    yMax = max(xEnv);
    patchX = [clickTimes(1,:);clickTimes(1,:);clickTimes(2,:);clickTimes(2,:)];
    patchY = repmat([-yMax;yMax;yMax;-yMax],1,nClicks);
    patchCol = [1,0.9,0.75];
    patchEdgeCol = [0.9,0.8,0.6];
    h_clicks = patch(patchX,patchY,patchCol,'EdgeColor',patchEdgeCol);
    uistack(h_clicks,'bottom')
    keyboard
    %}
end

%% importAudio_Detection --------------------------------------------------------------
function [x,tCutStart,tCutEnd] = importAudio_Detection(wavPath,sampleIndex,channelIndex,FsDetect,dFilt)
% Reads in a WAV file and processes it for click detection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1) read WAV file
    [x,FsRaw] = audioread(wavPath,sampleIndex);
    
    % 2) isolate wanted channel
    x = x(:,channelIndex);
    
    % 3) clip quiet periods at start and end.
    %   These periods seem to occur systematically with duty-cycled AMAR 
    %   recordings, likely from delayed activation and premature 
    %   termination of recording phases within each cycle. Unfortunately, 
    %   they mess with CABLE's click detector (due to recursive filters).
    nRaw = numel(x);
    iStartRaw = max([1,find(abs(x) > 0, 1, 'first')]);
    iEndRaw = min([nRaw,find(x ~= 0, 1, 'last')]);
    %x2 = x; % DEBUG
    x = x(iStartRaw:iEndRaw);
    tCutStart = (iStartRaw - 1)/FsRaw;
    tCutEnd = (nRaw - iEndRaw)/FsRaw;
    
    % 4) correct DC shift (ensure baseline is zero)
    x = x - mean(x);
    %x1 = x; %%% DEBUG
    %x2 = x2 - mean(x2); %%% DEBUG
    
    % 5) resample
    n_resamp = 25;
    [p,q] = rat(FsDetect/FsRaw);
    x = resample(x,p,q,n_resamp);
    %x2 = resample(x2,p,q,n_resamp); %%% DEBUG
    
    % 6) apply noise filter
    x = TWD_Common.filtfiltdf(dFilt,x);
    
    %%% DEBUG
    %{
    t1 =  tCutStart + (0:(numel(x)-1))/FsDetect;
    t2 = (0:(numel(x2)-1))/FsDetect;
    figure;
    plot(t1,x)
    hold on
    plot(t2,filtfilt(dFilt,x2));
    plot(t2,x2,'--');
    keyboard
    %}
    %%% DEBUG
    %{
    t1 = (1:numel(x1))/FsRaw;
    t2 = (1:numel(x))/FsDetect;
    figure;
    plot(t1,x1)
    hold on
    plot(t2,x)
    keyboard
    %}
end

%% extractClicks ----------------------------------------------------------
function clickData = extractClicks(wavPath,sampleIndex,channelIndex,clickTimes,dFilt,nfft)
% Computes and compiles click features for later discrimination.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - Much of this is replication of BWD code
% (compileClicks>computeClickData). Might be a good idea to use one
% common function eventually...
% - Since this is designed for partial compatibility with BWD (notably
% ClickDiscriminator), it is a little inefficient

    % load raw audio
    [x,Fs] = audioread(wavPath,sampleIndex,'native');
    x = x(:,channelIndex);
    totSamples = numel(x); % using "totSamples" to avoid confusion with the legacy parameter "nSamples"

    % define parameters for click extraction (BWD original)
    clickOffset_sec = [0.001, 0.004]; % length of buffer before and after click, in seconds (previously hardcoded to 250 and 999 samples)
    noiseOffset_sec = 0.005; % offset before start of click from which noise samples are computed, in seconds (previously hardcoded to 1250 samples)  
    
    % convert time-based variables to samples
    clickOffset_samples = round(clickOffset_sec*Fs);
    noiseOffset_samples = round(noiseOffset_sec*Fs);
    
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
    nClicks = size(clickTimes,2);
    clickPos.StartTimes = clickTimes(1,:);
    clickPos.EndTimes = clickTimes(2,:);
    
    % get start/end samples for each click
    clickPos.StartSamples = round(clickPos.StartTimes*Fs);
    clickPos.EndSamples = round(clickPos.EndTimes*Fs);
    clickPos.StartSamplesExtended = clickPos.StartSamples - clickOffset_samples(1) - 1;
    clickPos.EndSamplesExtended = clickPos.StartSamplesExtended + clickOffset_samples(2) - 1;
    
    % get range of noise samples for each click
    noisePos.StartSamples = clickPos.StartSamples - noiseOffset_samples;
    noisePos.EndSamples = clickPos.StartSamples;
    
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
    
    % isolate waveform and noise samples of each click 
    % (this part breaks from BWD a bit)
    nClickSamplesExtended = unique(clickPos.EndSamplesExtended - clickPos.StartSamplesExtended) + 1; % should be scalar
    nNoiseSamples = unique(noisePos.EndSamples - noisePos.StartSamples) + 1; % should be scalar
    xClicks = zeros(nClicks,nClickSamplesExtended,class(x)); % NaN is only for floating-point data types
    xClicksFilt = zeros(nClicks,nClickSamplesExtended);
    xNoise = zeros(nClicks,nNoiseSamples,class(x));
    xNoiseFilt = zeros(nClicks,nNoiseSamples);
    for ii = 1:nClicks
        xClickii = x(clickPos.StartSamplesExtended(ii):clickPos.EndSamplesExtended(ii));
        xNoiseii = x(noisePos.StartSamples(ii):noisePos.EndSamples(ii));
        
        % save native waveforms
        xClicks(ii,:) = xClickii;
        xNoise(ii,:) = xNoiseii;
        
        % filtered versions
        xClicksFilt(ii,:) = TWD_Common.filtfiltdf(dFilt,double(xClickii));
        xNoiseFilt(ii,:) = TWD_Common.filtfiltdf(dFilt,double(xNoiseii));
        %%% DEBUG
        %figure
        %ha1 = subplot(2,1,1);
        %plot(xClickii);
        %ha2 = subplot(2,1,2);
        %plot(xClicksFilt(ii,:));
        %linkaxes([ha1,ha2],'x');
        %fvtool(dFilt)
        %keyboard
        %%% END DEBUG
    end
    
    % compute click duration (ms)
    clickDur_ms = (clickPos.EndTimes - clickPos.StartTimes)*1000;
    
    % compute click features
    [specClick,specNoise,peakFr,F0,bw3db,bw10db,ppSignal,rmsSignal,rmsNoise,snr,slope,nSamples,ppSignalWB] =...
        TWD_Common.computeParameters(xClicksFilt,xNoiseFilt,clickDur_ms,Fs,nfft,clickOffset_samples(1));
    
    % save all features to struct.
    % Partially consistent with BWD/SBP's code.
    clickData = struct(...
        'bw10db',bw10db,...
        'bw3db',bw3db,...
        'dur',clickDur_ms',...
        'F0',F0,...
        'nSamples',nSamples,...
        'offset',clickOffset_samples(1),...
        'peakFr',peakFr,...
        'pos',[clickPos.StartTimes',clickPos.EndTimes'],...
        'ppSignal',ppSignal,... %'ppSignalWB',ppSignalWB,...
        'rmsNoise',rmsNoise,...
        'rmsSignal',rmsSignal,...
        'slope',slope,...
        'snr',snr,...
        'specClick',specClick,...
        'specNoise',specNoise,...
        'yFilt',xClicksFilt,...  %'xClicks',xClicks,...
        'yNFilt',xNoiseFilt); %'xNoise',xNoise);
end

%% compileMAT -------------------------------------------------------------
function compileMAT(dirPath_mat,depName,segStartTime,segDur,nClicksMax,clickData,isSperm,Fs)
% Creates one or more MAT files for one segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get data fields in clickData struct
    clickDataFields = fieldnames(clickData);
    clickDataFields = setdiff(clickDataFields,'offset');
    nClickFields = numel(clickDataFields);

    % get number of clicks
    nSpermClicks = sum(isSperm);
    
    if nSpermClicks > 0
        % split by number of clicks
        iClickPartStarts = (1:nClicksMax:nSpermClicks)';
        iClickPartEnds = unique([nClicksMax:nClicksMax:nSpermClicks,nSpermClicks])';
        iClickPartRanges = [iClickPartStarts,iClickPartEnds];
        nParts = size(iClickPartRanges,1);
        
        % count cumulative number of sperm whale clicks among all clicks
        clickNum = cumsum(isSperm);
        
        % initialize filename for this segment
        segStartDT = datetime(segStartTime);
        segStartDT.Format = 'yyyyMMdd_HHmmss';
        outFileNameStart = sprintf('%s_Pm_%s',depName,char(segStartDT));  
        
        % loop through each part 
        for ii = 1:nParts
            % isolate wanted clicks within this part
             iClickPartRangeii = iClickPartRanges(ii,:);
            includeii = isSperm & (clickNum >= iClickPartRangeii(1)) & clickNum <= iClickPartRangeii(2);
            
            % retain data for wanted clicks only
            clickDataii = clickData;
            for jj = 1:nClickFields
                fieldjj = clickDataFields{jj};
                clickDataii.(fieldjj) = clickDataii.(fieldjj)(includeii,:);
            end
            
            % add non-click info as fields
            clickDataii.fs = Fs;
            clickDataii.rawDur = segDur;
            clickDataii.rawStart = segStartTime;
            
            % add ICI just for compatibility
            ICI_ms = (clickDataii.pos(2:end,1) - clickDataii.pos(1:(end-1),1))*1000;
            clickDataii.ici = ICI_ms;
            
            % save MAT file
            outFileNameii = sprintf('%s_%d',outFileNameStart,ii);
            outFilePathii = fullfile(dirPath_mat,outFileNameii);
            save(outFilePathii,'-struct','clickDataii');
        end
    end
end

