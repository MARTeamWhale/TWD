%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "extractFilterClick"
%   Written by WB, based on "extracFilterTimeseriesAMAR250" by SBP and JS
%   Last updated Apr. 27, 2023, using MATLAB R2018b
%
%   Description:
%   Extracts audio samples from a WAV file for one click, and applies a 
%   butterworth filter in both directions (zero-phase filtering).
%   A key difference with SBP's function is that click and noise sample
%   positions are not recalculated, they must be entered as parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES
% - The double call to audioread is inefficient. If noise samples are
%   always taken immediately before the click, could use just one call.

function [yFiltClick, yNFiltClick] = extractFilterClick(...
    wavFilePath, clickPos, noisePos, filtdata)

    % define channel
    c = 1;
    
    % get click waveform
    [y, Fs] = audioread(wavFilePath, clickPos, 'native');
    y = double(y); 
    %%% original class is a signed int type corresponding to bits-per-
    %%% sample (usually 16 - so get "int16"). Might be worth keeping 
    %%% actually as it would reduce file size (int16 uses 4x fewer bytes 
    %%% than double). Would work well with a class implementation to get 
    %%% double on-demand. Fs could also be stored as uint32 (4 bytes).
    %%% UPDATE: this is not true. int16 class helps reduce memory use, but
    %%% not disk space use.
    y = y(:,c);
    y = y.'; % Other code must work with row vectors

    % get noise waveform
    yN = audioread(wavFilePath,noisePos,'native');
    yN = double(yN);
    yN = yN(:,c);
    yN = yN.';

    % create bandpass filter
    %%% isolate appropriate filter cutoff frequencies from list
    iFiltCutoff = filtdata.SamplingRate == Fs;
    switch sum(iFiltCutoff)
        case 1
            Fc1 = filtdata.Cutoff1;
            Fc2 = filtdata.Cutoff2(iFiltCutoff);
        case 0
            error('No bandpass filter cutoff frequencies have been specified for Fs=%d Hz. They must be added to "FilterCutoffs.mat" before click compilation can be run on this dataset.', Fs)
        otherwise
            error('More than one bandpass filter cutoff frequency specification exists for Fs=%d Hz; cannot determine which one to use.', Fs)
    end
    %{
    %%% filter cutoff frequencies
    Fc1 = 5000; 
    %Fc1 = 4000; % TMP CHANGED FOR SPERMS
    switch Fs
        case 96000
            Fc2 = 47000;
        case 200000
            Fc2 = 95000;
        case 250000
            Fc2 = 120000;
            %Fc2 = 24000; % TMP CHANGED FOR SPERMS
        case 128000
            Fc2 = 62000;
        case 192000
            Fc2 = 95000;
        case 375000
            Fc2 = 180000;
        case 300000 % added by WB to test Brazil recording
            Fc2 = 140000;
        otherwise
            error('Unsupported sampling rate')
    end
    %}
    %%% filter order (full, not half)
    ord = 10;
    %%% filter coefficients
    [B,A] = butter(ord/2, [Fc1,Fc2]/(Fs/2));
    
    % apply filter to click and noise
    yFiltClick = filtfilt(B,A,y);
    yNFiltClick = filtfilt(B,A,yN);
    
end