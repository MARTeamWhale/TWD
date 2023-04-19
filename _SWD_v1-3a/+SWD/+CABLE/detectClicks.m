%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function "detectClicks"
%   Written by Wilfried Beslin
%   Last Updated May 2018, using MATLAB version R2015a
%   Toolbox Dependencies:
%       Signal Processing Toolbox
%
%   Description:
%       Custom click detector optimized to detect multi-pulsed sperm whale
%       clicks for later classification. It is based on the Page test
%       algorithm (Page, 1954) used by PAMGUARD and described by
%       Miller (2010) and Zimmer (2011). Page test detections are validated
%       post-hoc to ensure that clicks are not too short, and that the
%       multi-pulsed structure is not broken up.
%
%   Input:
%       xEnv [n-by-1 double]:
%           Vector of waveform envelope amplitudes
%       Fs [1-by-1 double]:
%           Sampling rate in seconds
%       threshOn [1-by-1 double]:
%           Linear SNR threshold which controls the onset of click
%           detection events
%       threshOff [1-by-1 double]:
%           Linear SNR threshold which dictates the limits of a click range
%       alphaSignal [1-by-1 double]:
%           Smoothing factor for the exponential signal power estimation 
%           filter
%       alphaNoiseOn [1-by-1 double]:
%           Smoothing factor for the exponential noise power estimation
%           filter, used while a click has been detected
%       alphaNoiseOff [1-by-1 double]:
%           Smoothing factor for the exponential noise power estimation 
%           filter, used while only noise is present
%       minEchoProp [1-by-1 double]:
%           Minimum proportion of a click envelope peak beyond which the 
%           next click will be considered an echo
%       IPIRange [1-by-2 double]:
%           The minimum and maximum IPI limits, in milliseconds
%       maxClickDuration [1-by-1 double]:
%           Maximum allowed duration of a click, in milliseconds
%       minClickSep [1-by-1 double]:
%           Minimum expected separation between clicks, in milliseconds
%           Use only if there is just one whale present.
%       minPulseDuration [1-by-1 double]:
%           Minimum expected duration of individual pulses within sperm
%           whale clicks, in milliseconds
%
%   Output:
%       clickRanges [2-by-n double]:
%           Matrix of start and end samples for each click
%       pageNoise [n-by-1 double]:
%           Vector of noise power estimates for each sample
%
%   References:
%       Page, S. E. (1954). “Continuous inspection schemes,” Biometrika 41, 
%           100–115.
%       Miller, B. S. (2010). Acoustically derived growth rates and three-
%           dimensional localisation of sperm whales (Physeter 
%           macrocephalus) in Kaikoura, New Zealand [PhD Thesis] 
%           (University of Otago, Dunedin, New Zealand).
%       Zimmer, W. M. X. (2011). Passive Acoustic Monitoring of Cetaceans
%           (Cambridge University Press, Cambridge, UK).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% - "Separation" here differs from "Interval" in that it describes the 
%   space BETWEEN clicks/pulses. "Interval" refers to the difference in 
%   occurence times, and so includes the range of the first click/pulse.

function [clickRanges,pageNoise] = detectClicks(xEnv,Fs,threshOn,threshOff,alphaSignal,alphaNoiseOn,alphaNoiseOff,minEchoProp,IPIRange,maxClickDuration,minClickSep,minPulseDuration)

    % INITIALIZE VARIABLES
    
    %%% signal power
    xPow = xEnv.^2;
    clear xEnv
    
    %%% convert milliseconds to samples
    minClickSamples = ms2samples(IPIRange(1),Fs);
    maxClickSamples = ms2samples(maxClickDuration,Fs);
    minClickSepSamples = ms2samples(minClickSep,Fs);
    maxIPISamples = ms2samples(IPIRange(2),Fs);
    minPulseSamples = ms2samples(minPulseDuration,Fs);
    % END VARIABLE INITIALIZATION
    
    % 1) Page test
    [clickRanges,pageNoise] = pageTest(xPow,threshOn,threshOff,alphaSignal,alphaNoiseOff,alphaNoiseOn,maxClickSamples);
    
    % 2) Validate Page detections
    clickRanges = postPageTest(xPow,clickRanges,minPulseSamples,alphaSignal,minClickSepSamples,minEchoProp,maxIPISamples,maxClickSamples,minClickSamples);
end

%% pageTest ---------------------------------------------------------------
function [clickRanges,pageNoise] = pageTest(xPow,threshOn,threshOff,alphaSignal,alphaNoiseOff,alphaNoiseOn,maxClickSamples)
% Main click detection function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % INITIALIZE VARIABLES
    
    %%% number of sample
    nsamples = numel(xPow);
       
    %%% loop variables
    ndetects = 0;
    alphaNoiseii = alphaNoiseOff;
    powSignalii = xPow(1);
    powNoiseii = rms(xPow);
    signalPresent = false;
    Vprev = powSignalii/powNoiseii;
    istart = 1;
    
    %%% output containers
    clickRanges = zeros(2,round(nsamples/2));
    pageNoise = zeros(nsamples,1);
    % END VARIABLE INITIALIZATION
    
    % Start Page test
    for ii = 1:nsamples
        % update noise and signal power estimates
        powSignalii = powSignalii*(1-alphaSignal) + xPow(ii)*alphaSignal;
        powNoiseii = powNoiseii*(1-alphaNoiseii) + xPow(ii)*alphaNoiseii;
        pageNoise(ii) = powNoiseii;
        
        % calculate test statistic V
        Vii = powSignalii/powNoiseii;
        
        % determine state trigger decision
        if ~signalPresent % SIGNAL CURRENTLY ABSENT
            % update next detection start sample
            if (Vii > threshOff) && (Vprev <= threshOff)
                istart = ii;
            end
            % check if signal present, if so start new click
            if Vii > threshOn
                ndetects = ndetects + 1;
                clickRanges(1,ndetects) = istart;
                alphaNoiseii = alphaNoiseOn;
                signalPresent = true;
            end
        else % SIGNAL CURRENTLY PRESENT
            % check if signal no longer present or too long. If so, end it.
            if (Vii < threshOff) || (ii-istart == maxClickSamples)
                clickRanges(2,ndetects) = ii-1;
                alphaNoiseii = alphaNoiseOff;
                istart = ii+1;
                signalPresent = false;
            elseif ii == nsamples
                clickRanges(2,ndetects) = ii;
            end
        end
  
        % update V
        Vprev = Vii;
    end
    
    % remove unused spaces in clickRanges 
    clickRanges = clickRanges(:,1:ndetects);
end

%% postPageTest -----------------------------------------------------------
function clickRanges = postPageTest(xPow,clickRanges,minPulseSamples,alphaSignal,minClickSepSamples,minEchoProp,maxIPISamples,maxClickSamples,minClickSamples)
% Checks if Page-detected clicks meet certain criteria, and removes or 
% reprocess those that don't. This is how it works:
%   1) remove all detections with duration < minPulseSamples
%   2) shift the start and end point of all detections to account for 
%       signal filter delay
%   3) Loop through each click and check the following:
%       3.1) If minClickSepSamples > 0: Is click too close to previous?
%           YES: remove --> next click
%           NO: --> 3.2)
%       3.2) Should click be merged with the next one?
%           YES: merge until no longer necessary, then --> 3.3)
%           NO: --> 3.3)
%       3.3) Is click too short to be a full click?
%           YES: remove --> next click
%           NO: --> next click
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % initialize variables
    % NOTE: validationFuncs is an array of nested functions that can share 
    % data. The order in which each function is called is important!!
    if minClickSepSamples == 0
        validationFuncs = {...
            @merge,...
            @validateShortness};
    else
        validationFuncs = {...
            @validateSeparation,...
            @merge,...
            @validateShortness};
    end
    nsamples = numel(xPow);
    ndetects = size(clickRanges,2);
    mergeProp = minEchoProp.^2;
    m = 1;

    % 1) Remove really small clicks (too small to be pulses)
    tooSmall = clickDurations(1:ndetects) < minPulseSamples;
    clickRanges(:,tooSmall) = [];
    ndetects = ndetects - sum(tooSmall);
    %dbmsg(vb,sprintf('Removed %d clicks',sum(tooSmall)))

    % 2) Correct for signal filter delay
    signalDelay = round(max(grpdelay(alphaSignal,[1,alphaSignal-1])));
    clickRanges = clickRanges - signalDelay;
    clickRanges(clickRanges<1) = 1;
    clickRanges(clickRanges+signalDelay == nsamples) = nsamples; % in case last signal may have been cut short
    %dbmsg(vb,sprintf('Click starts shifted by %d samples',signalDelay))
    
    % 3) Loop through and validate each click
    while m <= ndetects
        %dbmsg(vb,sprintf('m=%d, ndetects=%d, iteration%d', m, ndetects, iter))
        
        cont = false;
        for ii = 1:numel(validationFuncs)
            feval(validationFuncs{ii})
            if cont
                break
            end
        end

        % Next iteration
        m = m+1;
    end
    %dbmsg(vb,'Click validation complete.')

    % NESTED VALIDATION FUNCTIONS .........................................

    % 3.1) Remove clicks which occur too soon after the previous one
    function validateSeparation()
        %dbmsg(vb,'   Checking separation...')
        if m > 1
            if clickRanges(1,m) - clickRanges(2,m-1) < minClickSepSamples
                clickRanges(:,m) = [];
                ndetects = ndetects-1;
                m = m-1;
                cont = true;
                %dbmsg(vb,sprintf('     Click removed for being within %g samples of previous',minClickSepSamples))
            end
        end
    end

    % 3.2) If next click is within maxIPI of the current one, and the 
    % interval between the start of the current click and end of the next 
    % click is less than maxClickSamples, and the peak of the next click 
    % isn't a significant fraction of the current peak, then merge them. 
    % Keep doing this until the condition is no longer true.
    function merge()
        %dbmsg(vb,'   Assessing merge conditions...')

        trymerge = true;
        while trymerge
            trymerge = false;
            % if this isn't the last click...
            if m < ndetects
                nextClickSep = clickRanges(1,m+1) - clickRanges(2,m);
                % ...and the interval between this click and the next is 
                % within maxIPISamples, and the length of the click that would 
                % result should they be merged is within maxClickSamples...
                if (nextClickSep < maxIPISamples) &&...
                        (clickDurations(m) + nextClickSep + clickDurations(m+1) <= maxClickSamples)
                    currentPeak = max(xPow(clickRanges(1,m):clickRanges(2,m)));
                    nextPeak = max(xPow(clickRanges(1,m+1):clickRanges(2,m+1)));
                    % ...and the peak of the next click is not a
                    % significant fraction of the peak of the current 
                    % click...
                    if (nextPeak/currentPeak < mergeProp) 
                        % ...then merge them.
                        clickRanges(2,m) = clickRanges(2,m+1);
                        clickRanges(:,m+1) = [];
                        ndetects = ndetects-1;
                        trymerge = true;
                        %dbmsg(vb,['     Click merged with next; ',num2str(ndetects),' detections remain'])
                    end
                end
            end
        end
    end

    % 3.3) remove clicks that are too short to be full clicks
    function validateShortness()
        %dbmsg(vb,'   Checking shortness...')
        if clickDurations(m) < minClickSamples
            clickRanges(:,m) = [];
            ndetects = ndetects-1;
            m = m-1;
            cont = true;
            %dbmsg(vb,sprintf('     Click removed for having < %g samples',minClickLength))
        end
    end

    % NESTED UTILITY FUNCTIONS ............................................

    % compute durations of clicks
    function duration = clickDurations(clickNum)
        duration = diff(clickRanges(:,clickNum))+1;
    end
end

%% ms2samples -------------------------------------------------------------
function sampDur = ms2samples(msDur,Fs)
% Converts time in milliseconds to samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES: thinking of making this a more general utility function

    sampDur = round(Fs*(msDur/1000));
end