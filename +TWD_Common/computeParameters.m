%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "computeParameters"
%   Written by WB, based on "parametersAMAR" by SBP and JS
%   Last updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   This is basically "parametersAMAR" rewritten in the style of the new
%   code. Consistency and compatibility were priorities, so the output from 
%   this version is virtually identical to the old one (occasionally there 
%   are tiny differences likely due to precision issues, e.g. on the order 
%   of 1e-15). In the future, I would like to move to new calculations that
%   produces more appropriate output for the new detector.
%
%   Original header:
%   "Take timeseries out of existing file, convert from normalized data to counts
%   1) calculate spectral received levels RL for click and preceding noise: 
%   calculate spectra, account for bin width to reach dB re counts^2/Hz, compute peak frequency and bandwidth
%   2) calculate RLpp at peak frequency: find min & max value of timeseries, convert to dB
%   3) calculate signal to noise ratio
%
%   sb 090814 % modified by JES for DFO AMAR data (Oct 2018)"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - I find the ppSignal calculation to be very odd. Maybe I should try to
% add my own on the side...

function [specClick,specNoise,peakFr,F0,bw3db,bw10db,ppSignal,rmsSignal,rmsNoise,snr,slope,nSamples,ppSignalWB] =...
    computeParameters(yFilt,yNFilt,dur,fs,N,offset)

    % rename variables for consistency with SBP - change this one day
    xClick = yFilt;
    xNoise = yNFilt;
    nfft = N;
    Fs = fs;
    
    % initialize other variables
    nClicks = size(xClick,1);
    specClick = zeros(nClicks,nfft/2);
    specNoise = zeros(nClicks,nfft/2);
    peakFr = zeros(nClicks,1);
    F0 = zeros(nClicks,1);
    bw3db = zeros(nClicks,3);
    bw10db = zeros(nClicks,3);
    ppSignal = zeros(nClicks,1);
    ppSignalWB = zeros(nClicks,1);
    rmsSignal = zeros(nClicks,1);
    rmsNoise = zeros(nClicks,1);
    snr = zeros(nClicks,1);
    slope = zeros(nClicks,2);
    nSamples = zeros(nClicks,1);
    
    %%% define additional offset relative to TK start/end times to get 
    %%% click position. SBP had hardcoded 30 samples for this for 200 kHz 
    %%% HARP, which was not corrected for 250 kHz AMAR. Here I am adding an
    %%% Fs correction, but to keep things consistent with the 
    %%% AMAR analysis already underway, 250 kHz is used as the reference.
    xBuffer = round(30*(Fs/250000));
    
    % define frequency vector in kHz
    FMax = Fs/2;
    deltaF = FMax/((nfft/2)-1);
    F = 0:deltaF:FMax;
    
    % SBP comment:
    % 1) calculate spectral received levels RL: calculate spectra, account 
    % for bin width (750 Hz) to reach dB re counts^2/Hz
    
    % loop through each click
    for ii = 1:nClicks
        
        %SBP: "pick click as defined by teager energy duration + 30 pts on each side"
        %%% number of points is now corrected for Fs
        xClickii = xClick(ii,:);
        xNoiseii = xNoise(ii,:);
        durii = dur(ii);
        iStartii = (offset+1) - xBuffer; % SBP: "timeseries has points picked before start of click (offset)"
        iEndii = round((offset+1) + (durii/1000*Fs+xBuffer));
        
        % Limit the end point to number of FFT samples
        %%% This means any clicks that are too long will be clipped...
        if iEndii - iStartii > nfft
            iEndii = offset - xBuffer + nfft;
        end
        
        % Joy's quick fix for clicks with longer duration than the
        % extracted window... 
        if iEndii > length(xClickii)
            iEndii = length(xClickii);
        end
        
        % compute power spectrum (log)
        [specClickii,specNoiseii] = computeSpectrum(xClickii,xNoiseii,iStartii,iEndii,nfft,Fs);
        specClick(ii,:) = specClickii;
        specNoise(ii,:) = specNoiseii;
        
        % compute peak frequency (kHz)
        FPeakii = compute_FPeak(specClickii,Fs);
        peakFr(ii) = FPeakii;
        
        % compute centroid frequency (kHz)
        F0ii = compute_F0(specClickii,F);
        F0(ii) = F0ii;
        
        % compute -3 dB bandwidth
        bw3dbii = compute_bwdb(specClickii,Fs,-3);
        bw3db(ii,:) = bw3dbii;
        
        % compute -10 dB bandwidth
        bw10dbii = compute_bwdb(specClickii,Fs,-10);
        bw10db(ii,:) = bw10dbii;
        
        % compute peak-to-peak signal level
        ppSignalii = compute_ppSignal(xClickii,bw10dbii);
        ppSignal(ii) = ppSignalii;
        %%% WB's version
        ppSignalWBii = compute_ppSignalWB(xClickii,iStartii,iEndii);
        ppSignalWB(ii) = ppSignalWBii;
        
        % compute RMS level and SNR
        [rmsSignalii,rmsNoiseii,snrii] = compute_rms(xClickii,xNoiseii,durii,Fs);
        rmsSignal(ii) = rmsSignalii;
        rmsNoise(ii) = rmsNoiseii;
        snr(ii) = snrii;
        
        % compute slope and nSamples
        [slopeii,nSlopeSamplesii] = compute_slope(xClickii,durii,offset,Fs);
        slope(ii,:) = slopeii;
        nSamples(ii) = nSlopeSamplesii;
    end
end

%% computeSpectrum --------------------------------------------------------
function [specClick,specNoise] = computeSpectrum(xClick,xNoise,iStart,iEnd,nfft,Fs)
% Returns the log power spectrum of a click and its noise sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % window click
    winClick = hann(iEnd - iStart + 1)';
    xClickWin = xClick(iStart:iEnd).*winClick;
    nZero = nfft - numel(xClickWin);
    zeroPad = zeros(1,nZero);
    xClickWin = [xClickWin,zeroPad];

    % do FFT
    XMagClick = abs(fft(xClickWin,nfft));
    XPowLogClick = 20*log10(XMagClick);
    
    % window noise
    % JES/SBP comment: "EDITED to start noise window exactly 0.005 s before
    % calculated click start"
    %%% This comment is not really of concern here... this was addressed
    %%% during noise extraction
    winNoise = hann(nfft)';
    if numel(xNoise) >= nfft
        xNoiseWin = xNoise(1:nfft).*winNoise;
    else
        xNoiseWin = zeros(1,nfft);
        xNoiseWin(1:numel(xNoise)) = xNoise;
        xNoiseWin = xNoiseWin.*winNoise;
    end
    
    % do FFT
    XMagNoise = abs(fft(xNoiseWin,nfft));
    XPowLogNoise = 20*log10(XMagNoise);
    
    %SBP: "account for bin width"
    %%% I'm not sure what the significance of this is, but will keep as is 
    %%% for consistency.
    specSub = 10*log10(Fs/nfft);
    specClickSub = XPowLogClick - specSub;
    specNoiseSub = XPowLogNoise - specSub;
        
    % reduce to first half of spectra
    specClick = specClickSub(:,1:nfft/2);
    specNoise = specNoiseSub(:,1:nfft/2); 
end

%% compute_FPeak ----------------------------------------------------------
function FPeak = compute_FPeak(spec,Fs)
% Computes peak frequency of a click
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    n = numel(spec);
    [~,iSpecPeak] = max(spec);
    FPeak = iSpecPeak/n*(Fs/(2*1000));
    %FPeak = F(iSpecPeak)/1000; % slightly different results from SBP; 
    % would be more consistent with F0 calculation though
end

%% compute_F0 -------------------------------------------------------------
function F0 = compute_F0(spec,F)
% Computes centroid frequency of a click
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % convert spectrum to linear scale
    specLin = 10.^(spec/20);
    
    % compute centroid frequency
    F0 = (sum(F.*specLin.^2)/sum(specLin.^2))/1000; % Au 1993, equation 10-3
    
end

%% compute_bwdb -----------------------------------------------------------
function bwdb = compute_bwdb(spec,Fs,dBVal)
% Computes a decibel bandwidth (e.g. -3 or -10) for one click
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = numel(spec);

    % get threshold
    [specPeak,iSpecPeak] = max(spec);
    th = specPeak - abs(dBVal);
    
    % find points where spectrum is above or below threshold
    belowThr = spec < th;
    
    % compute intercept points
    % find 1st element to the left of the peak where power 
    % drops below threshold
    %iStartAboveThr = find(belowThr & 1:n < iSpecPeak, 1,'last') + 1;
    iStartAboveThr = find(belowThr & 1:n < iSpecPeak, 1,'last') - 1; % for consistency with SBP
    iStartAboveThr = max([1,iStartAboveThr]);

    % find 1st element to the right of the peak where power 
    % drops below threshold
    %iEndAboveThr = find(belowThr & 1:n > iSpecPeak, 1,'first') - 1;
    iEndAboveThr = find(belowThr & 1:n > iSpecPeak, 1,'first') + 1; % for consistency with SBP
    iEndAboveThr = min([n,iEndAboveThr]);
    
    % compute bandwidth (lower limit, upper limit, and magnitude)
    bwdbLow = (Fs/(2*1000))*iStartAboveThr/n;
    bwdbUpp = (Fs/(2*1000))*iEndAboveThr/n;
    bwdbMag = bwdbUpp - bwdbLow;
    
    %%% Interesting discovery:
    % This operation:
    %   bwdbLow = (Fs/(2*1000))*(iStartAboveThr/numel(spec));
    % Does NOT produce quite the same results as this:
    %   bwdbLow = (Fs/(2*1000))*iStartAboveThr/numel(spec);
    % There must be a precision loss somewhere
    %%%
    
    bwdb = [bwdbLow,bwdbUpp,bwdbMag];
end

%% compute_ppSignal -------------------------------------------------------
function ppSignal = compute_ppSignal(x,bw10db)
% Computes peak-to-peak recieved level of one click
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get min and max magnitudes
    xMax = max(x);
    xMin = min(x);
    ppCount = xMax + abs(xMin); % only correct if xMin is negative!
    
    % convert to dB
    ppSignal = 20*log10(ppCount);
    
    %SBP: "correct sound pressure levels on a plot of ambient noise levels,
    % calculate -10 dB bandwidth and add 10 log (bandwidth) to PppAll"
    %%% I don't think ambient noise correction is done anymore.
    %%% -10dB correction is... but why?
    bw10log = 10*log10(bw10db(3)*1000);
    ppSignal = ppSignal + bw10log;
end

%% compute_ppSignalWB -------------------------------------------------------
function ppSignal = compute_ppSignalWB(x,iStart,iEnd)
% Computes peak-to-peak recieved level of one click
% WB's version - literaly just difference between max and min, but also
% takes into account click range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get min and max magnitudes
    xSub = x(iStart:iEnd);
    xMax = max(xSub);
    xMin = min(xSub);
    ppCount = xMax - xMin;
    
    % convert to dB
    ppSignal = 20*log10(ppCount);
end

%% compute_rms ------------------------------------------------------------
function [rmsSignal,rmsNoise,snr] = compute_rms(xClick,xNoise,dur,Fs)
% Computes RMS level of signal and noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % click
    nClickMax = round(355*(Fs/250000));
    nClickActual = round((dur/1000)*Fs);
    nClick = min([nClickActual,nClickMax]);
    
    rmsClickLin = sqrt(sum(xClick(101:101+nClick).^2)/nClick);
    rmsSignal = 20*log10(rmsClickLin);
    
    
    % noise
    nNoise = numel(xNoise);
    
    rmsNoiseLin = sqrt(sum(xNoise.^2)/nNoise);
    rmsNoise = 20*log10(rmsNoiseLin);
    
    
    % SNR
    snr = rmsSignal - rmsNoise;
end

%% compute_slope ----------------------------------------------------------
function [slope,nSlopeSamples] = compute_slope(x,dur,offset,Fs)
% Fits a line to the frequencies with greatest amplitude in the spectrogram
% (-8 dB threshold).
% This is essentially a rewrite of "specFit" from the original code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set thresholds and spectrogram settings
    dbTh = -8;
    winLength_ms = 0.2; % window and NFFT size for spectrogram, in milliseconds
    ovl = 39/40; % 98% overlapping factor
    
    % transform spectrogram parameters
    winLength = round(winLength_ms*(Fs/1000));
    novl = min([round(winLength*ovl), winLength-1]);
    halfL = 2*winLength; %%% not really sure why
    
    % buffer around click
    nMax = numel(x) - offset;
    n = round(dur*(Fs/1000));
    n = min(n,nMax);
    d = floor(halfL-(n/2));
    
    % sample waveform
    iStart = offset - d + 1;
    iEnd = offset + n + d; % equivalent result to the original calculation in specFit
    xSub = x(iStart:iEnd);
    
    % generate spectrogram
    [~,F,T,XPow] = spectrogram(xSub,winLength,novl,winLength,Fs);
    
    % convert time to ms and frequencies to kHz
    T = T*1000;
    F = F/1000;
    nT = numel(T);
    
    % get vectors of maximum frequencies and their magnitudes
    % (in log scale)
    %[XPowPeaks,iFPeaks] = max(XPow,[],1);
    [XPowPeaks,iFPeaks] = max(XPow(3:end,:),[],1); % SBP-consistent version... but why remove 1st 3 frequencies?
    %FPeaks = F(iFPeaks)'; % time series of peak frequencies; not used by SBP-consistent code, but useful in future
    XPLogPeaks = 10*log10(XPowPeaks);
    
    % isolate significant bandwidth
    [XPLogMax,itFMax] = max(XPLogPeaks);
    ETh = XPLogMax - abs(dbTh);
    beforeFMax = 1:nT < itFMax;
    afterFMax = 1:nT > itFMax;
    belowTh = XPLogPeaks < ETh;
    %iFLow = max([1,find(beforeFMax & belowTh,1,'last') + 1]);
    %iFHigh = min([nT,find(afterFMax & belowTh,1,'first') - 1]);
    iFLow = max([1,find(beforeFMax & belowTh,1,'last')]); % SBP-consistent version
    iFHigh = min([nT,find(afterFMax & belowTh,1,'first')]); % SBP-consistent version
    %FPeaksIntense = FPeaks(iFLow:iFHigh);
    %TPeaksIntense = T(iFLow:iFHigh);
    TPeaksIntense = T(1:numel(iFLow:iFHigh)); % SBP-consistent version
    FPeaksIntense = iFPeaks(iFLow:iFHigh)*F(2); % SBP-consistent version (what???)
    
    % compute slope
    fitLine = polyfit(TPeaksIntense,FPeaksIntense,1);
    %slope = fitLine(1);
    slope = fitLine; % mainly to avoid confusion and keep consistent with SBP's format
    
    % compute slope duration
    nSlopeSamples = numel(FPeaksIntense);
    %slopeDur = nSlopeSamples/(Fs/1000); % move to this in future
end