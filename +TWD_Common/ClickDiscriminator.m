%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class "ClickDiscriminator"
%   Written by Wilfried Beslin
%   Last update Apr. 12, 2023, using MATLAB R2018b
%
%   Description:
%   Class for controlling all aspects of click discrimination. Major uses 
%   include:
%   - defines the discrimination functions
%   - manages properties of discrimination criteria (i.e. loading from 
%   file, storing, updating, retrieving)
%   - evaluates a dataset against the discrimination functions. In other
%   words, use this class to determine if clicks correspond to target
%   clicks or not.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef ClickDiscriminator < handle
    %% PROPERTIES =========================================================
    properties (SetAccess = private)
        critTable
        missTol % array of miss tolerance numbers
    end
    properties (Access = private, Constant)
        % This property associates criteria names with their evaluation
        % functions. MAKE SURE IT IS CONSISTENT WITH THE EVALUATION
        % FUNCTIONS DEFINED AT THE END ("METHODS - PUBLIC STATIC").
        % The function handles will be used by 'feval' in a loop to assess 
        % each criterion one-by-one. Less cumbersome than calling all the 
        % evaluation functions explicitly in sequence, especially if not 
        % all criteria are used. It would also be possible (and easier) to 
        % use 'feval' on function name strings instead of handles, but this 
        % is less efficient and apparently a security risk. 
        critFuncs = struct(...
            'Fpeak',        @TWD_Common.ClickDiscriminator.assess_Fpeak,...
            'F0',           @TWD_Common.ClickDiscriminator.assess_F0,...
            'bw3db',        @TWD_Common.ClickDiscriminator.assess_bw3db,...
            'bw3dbLower',   @TWD_Common.ClickDiscriminator.assess_bw3dbLower,...
            'bw3dbUpper',   @TWD_Common.ClickDiscriminator.assess_bw3dbUpper,...
            'bw10db',       @TWD_Common.ClickDiscriminator.assess_bw10db,...
            'bw10dbLower',  @TWD_Common.ClickDiscriminator.assess_bw10dbLower,...
            'bw10dbUpper',  @TWD_Common.ClickDiscriminator.assess_bw10dbUpper,...
            'ZCR',          @TWD_Common.ClickDiscriminator.assess_ZCR,...
            'dur',          @TWD_Common.ClickDiscriminator.assess_dur,...
            'slope',        @TWD_Common.ClickDiscriminator.assess_slope,...
            'slopeDur',     @TWD_Common.ClickDiscriminator.assess_slopeDur,...
            'durE50',       @TWD_Common.ClickDiscriminator.assess_durE50,...
            'deltaE',       @TWD_Common.ClickDiscriminator.assess_deltaE,...
            'nSamples',     @TWD_Common.ClickDiscriminator.assess_nSamples);
    end
    
    %% METHODS - PUBLIC ===================================================
    methods (Access = public)
        %% Constructor ----------------------------------------------------
        function obj = ClickDiscriminator(varargin)
        % Creates a "ClickDiscriminator" object.
        % Empty creation is allowed.
        % Input: "critParamsFilePath"
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            narginchk(0,1)
            
            % initialize properties
            critNames = fieldnames(obj.critFuncs);
            nCrit = numel(critNames);
            critMat = repmat([-Inf,Inf,0],nCrit,1);
            obj.critTable = array2table(critMat,'RowNames',critNames,'VariableNames',{'Threshold1','Threshold2','UseCategory'});
            obj.missTol = double.empty(0,1);
            
            % load criteria if any were specified
            if nargin > 0
                critParamsFilePath = varargin{:};
                obj.loadCritParams(critParamsFilePath)
            end
        end
        
        %% loadCritParams -------------------------------------------------
        function loadCritParams(obj,critParamsFilePath)
        % Sets criteria properties from an Excel file.
        % Sheet 1 is list of criteria, their thresholds and categories.
        % Sheet 2 is miss tolerances for each category.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % 1) read criteria table
            inCritTableRaw = readtable(critParamsFilePath,'Sheet',1,'ReadRowNames',true);
            
            % make sure the criteria names are consistent, remove the ones
            % that are unknown.
            critNames = fieldnames(TWD_Common.ClickDiscriminator.critFuncs);
            inCritNamesRaw = inCritTableRaw.Properties.RowNames;
            badCrit = ~ismember(inCritNamesRaw,critNames);
            badCritNames = inCritNamesRaw(badCrit);
            for ii = 1:sum(badCrit)
                warning('Ignoring unknown click discrimination criterion "%s"',badCritNames{ii})
            end
            inCritTable = inCritTableRaw(~badCrit,:);
            inCritNames = inCritTable.Properties.RowNames;
            
            % convert threshold cellstrs to numeric
            for ii = 1:2
                fieldii = ['Threshold',num2str(ii)];
                if iscell(inCritTable.(fieldii))
                    inCritTable.(fieldii) = str2double(inCritTable.(fieldii));
                end
            end
            
            % set final table
            obj.critTable(inCritNames,:) = inCritTable;
            
            
            % 2) read miss tolerance table
            missTolTable = readtable(critParamsFilePath,'Sheet',2);
            
            % convert to array
            missTolVec = zeros(max(missTolTable.Category),1);
            missTolVec(missTolTable.Category) = missTolTable.N;
            obj.missTol = missTolVec;
        end
        
        %% setCritProperties ----------------------------------------------
        function setCritProperties(obj,critName,varargin)
        % Call this method to update any property of a criterion (i.e.
        % Threshold1, Threshold2, or UseCategory). Specify properties as
        % Name-Value pairs.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % parse input
            assert(ismember(critName,fieldnames(obj.critFuncs)),sprintf('"%s" is not a valid criterion',critName))
            oldVals = obj.critTable{critName,:};
            p = inputParser;
            p.addParameter('Threshold1',oldVals(1))
            p.addParameter('Threshold2',oldVals(2))
            p.addParameter('UseCategory',oldVals(3))
            p.parse(varargin{:})
            
            % set new values
            newVals = [...
                p.Results.Threshold1,...
                p.Results.Threshold2,...
                p.Results.UseCategory];
            obj.critTable{critName,:} = newVals;
        end
        
        %% setMissTol -----------------------------------------------------
        function setMissTol(obj,category,n)
        % Call this method to update the miss tolerance for a category.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            validateattributes(n,{'numeric'},{'nonnegative','integer','scalar'})
            obj.missTol(category) = n;
        end
        
        %% evalTargetClicks -----------------------------------------------
        function [isTarget,passTable,clickParams] = evalTargetClicks(obj,data)
        % Major function - assesses if clicks meet the criteria for being
        % from the target species or not.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % determine which criteria to test
            use = obj.critTable.UseCategory > 0;
            critNames = obj.critTable.Properties.RowNames(use);
            critCats = obj.critTable.UseCategory(use);
            nCrit = sum(use);
            cats = unique(critCats);
            
            % initialize output
            nClicks = numel(data.dur);
            passTable = array2table(false(nClicks,nCrit),'VariableNames',critNames);
            valMat = NaN(nClicks,nCrit);
            
            % loop through each criterion and evaluate it
            for ii = 1:nCrit
                critii = critNames{ii};
                fcnii = obj.critFuncs.(critii);
                thii = [obj.critTable.Threshold1(critii),obj.critTable.Threshold2(critii)];
                [passii,valii] = feval(fcnii,data,thii);
                
                passTable.(critii) = passii;
                valMat(:,ii) = valii;
            end
            
            % get summary
            %%% count criteria that have passed in each category
            isTarget = true(nClicks,1);
            for ii = cats'
                iCritsii = critCats == ii;
                nCritsii = sum(iCritsii);
                passMatii = passTable{:,iCritsii};
                nPassedii = sum(passMatii,2);
                passii = nPassedii >= (nCritsii - obj.missTol(ii));
                
                % update isTarget
                isTarget = isTarget & passii;
            end
            
            % extra data: compile table of feature values
            clickParams = array2table(valMat,'VariableNames',critNames);
        end
    end
    
    %% METHODS - PUBLIC STATIC ============================================
    % These are the criteria functions. 
    % If you want to add a new criterion, implement it here. Make sure the
    % function names match the criteria names, otherwise the public
    % evaluation function won't recognize them.
    methods (Static)

        % assess_Fpeak ----------------------------------------------------
        function [pass,val] = assess_Fpeak(data,th)
            val = data.peakFr;
            pass = val >= th(1) & val <= th(2);
        end

        % assess_F0 -------------------------------------------------------
        function [pass,val] = assess_F0(data,th)
            val = data.F0;
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_bw3db ----------------------------------------------------
        function [pass,val] = assess_bw3db(data,th)
            val = data.bw3db(:,3);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_bw3dbLower -----------------------------------------------
        function [pass,val] = assess_bw3dbLower(data,th)
            val = data.bw3db(:,1);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_bw3dbUpper -----------------------------------------------
        function [pass,val] = assess_bw3dbUpper(data,th)
            val = data.bw3db(:,2);
            pass = val >= th(1) & val <= th(2);
        end

        % assess_bw10db ---------------------------------------------------
        function [pass,val] = assess_bw10db(data,th)
            val = data.bw10db(:,3);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_bw10dbLower ----------------------------------------------
        function [pass,val] = assess_bw10dbLower(data,th)
            val = data.bw10db(:,1);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_bw10dbUpper ----------------------------------------------
        function [pass,val] = assess_bw10dbUpper(data,th)
            val = data.bw10db(:,2);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_ZCR ------------------------------------------------------
        function [pass,val] = assess_ZCR(data,th)
        % measures zero crossing rate (based on TK duration)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % define variables
            fs = data.fs;
            offset = data.offset;
            yFilt = data.yFilt;
            dur = data.dur;
            [nClicks,nClickSamplesMax] = size(yFilt);
            
            % get click start and end samples
            nClickSamples = round(dur.*(fs/1000));
            clickStarts = repmat(offset+1,nClicks,1);
            clickEnds = nClickSamples + offset;
            clickEnds(clickEnds > nClickSamplesMax) = nClickSamplesMax;
            
            % loop through each click and get ZCR
            val = NaN(nClicks,1);
            for ii = 1:nClicks
                
                % isolate click based on TK duration
                durii = dur(ii);
                clickStartii = clickStarts(ii);
                clickEndii = clickEnds(ii);
                xClickii = yFilt(ii,clickStartii:clickEndii);
                
                % find zero crossing events
                ZC1ii = false; % 1st element of course cannot be a zero crossing
                ZCRemainii = xClickii(1:end-1).*xClickii(2:end) < 0; % Definition of zero crossing events
                ZCii = [ZC1ii,ZCRemainii];
                
                % compute rate
                nZCii = sum(ZCii);
                val(ii) = nZCii/durii;
            end
            
            % evaluate
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_dur ------------------------------------------------------
        function [pass,val] = assess_dur(data,th)
            val = data.dur;
            pass = val >= th(1) & val <= th(2);
        end

        % assess_slope ----------------------------------------------------
        function [pass,val] = assess_slope(data,th)
            val = data.slope(:,1);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_slopeDur -------------------------------------------------
        function [pass,val] = assess_slopeDur(data,th)
        % slopeDur is the period (in ms) where the spectrogram has 
        % significant energy (based on -8 dB threshold at time of writing).
        % This is a time-based version of what SBP called "nSamples".
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Fs = data.fs;
            slopeDur_samples = data.nSamples;
            
            val = slopeDur_samples/(Fs/1000);
            pass = val >= th(1) & val <= th(2);
        end
        
        % assess_durE50 ---------------------------------------------------
        function [pass,val] = assess_durE50(data,th)
        % Assesses duration based on 50% envelope energy.
        % This is an adaptation of the 2nd half of "bw_elim_false_AMAR250".
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % set thresholds
            energyThr = 0.5;
            shortDurThr_ms = th(1); %SBP ORIGINAL: 0.095; % "smaller = dolphin click"
            longDurThr_ms = th(2); %SBP ORIGINAL: 0.35;   % longer = strange dolphin click
            
            % define variables
            fs = data.fs;
            offset = data.offset;
            yFilt = data.yFilt;
            dur = data.dur;
            [nClicks,nClickSamplesMax] = size(yFilt);
            
            % loop through each click
            val = NaN(nClicks,1);
            for ii = 1:nClicks
                % compute envelope
                nClickSamplesii = round(dur(ii)*(fs/1000));
                xClickii = yFilt(ii,:);
                posStartii = offset + 1;
                posEndii = offset + nClickSamplesii;
                posEndii = min([nClickSamplesMax,posEndii]);
                clickExtractii = xClickii(posStartii:posEndii);
                xEnvii = abs(hilbert(clickExtractii));
                
                % calculate when envelope energy reaches 50% first and when it drops below last
                %%% actually relative energy - 50% within the envelope only
                xEnvNormii = xEnvii - min(xEnvii);
                xEnvNormii = xEnvNormii/max(xEnvNormii);
                aboveThrii = xEnvNormii >= energyThr;
                
                % find the first value above 50% energy with positive slope
                % and find the last above 50% energy with negative slope.
                %%% Original calculation was quite convoluted, so I'm
                %%% taking a different, simpler approach here. Results 
                %%% could potentially differ a bit in some cases, but 
                %%% probably for the best.
                %%% One thing to note, both here and in the original
                %%% calculation: gaps in energy are ignored.
                
                % get sign of slope for whole envelope
                envSlopeii = diff(xEnvNormii);
                envSlopeDirection1ii = [0,sign(envSlopeii)]; % instantaneous slope (lookback version; 1st element has a slope of 0)
                envSlopeDirection2ii = [sign(envSlopeii),0]; % instantaneous slope (lookahead version; last element has a slope of 0)
                %%% added the lookback/lookahead versions to patch a rare
                %%% bug that occured when the last sample above threshold
                %%% had instantaneous rising slope (but next point was
                %%% decreasing)
                
                % find 1st element where envelope is above threshold and
                % slope increases
                iStartAboveThrii = find(aboveThrii & (envSlopeDirection1ii > 0 | envSlopeDirection2ii > 0), 1,'first');
                if isempty(iStartAboveThrii)
                    % if there are no solutions (rare), just use 1st 
                    % sample above threshold
                    iStartAboveThrii = find(aboveThrii(1));
                end
                
                % find last element where envelope is above threshold and
                % slope decreases
                iEndAboveThrii = find(aboveThrii & (envSlopeDirection1ii < 0 | envSlopeDirection2ii < 0), 1,'last');
                if isempty(iEndAboveThrii)
                    % if there are no solutions (rare), just use last
                    % sample above threshold
                    iEndAboveThrii = find(aboveThrii(end));
                end
                
                % compute 50% energy duration (in samples)
                nDurii = (iEndAboveThrii - iStartAboveThrii) + 1;
                
                % convert to milliseconds
                val(ii) = nDurii/(fs/1000);
                
                % DEBUG
                %{
                if ismember(ii,[5933])
                    n = numel(xEnvii);
                    energyThrAbs = max(xEnvii)/2;
                    envSlopeDirection3ii = envSlopeDirection1ii;
                    envSlopeDirection3ii(envSlopeDirection1ii ~= envSlopeDirection2ii) = 0;
                    envSlopeDirection3ii(1) = envSlopeDirection2ii(1);
                    envSlopeDirection3ii(end) = envSlopeDirection1ii(end);
                    slopeRise = envSlopeDirection3ii > 0;
                    slopeFall = envSlopeDirection3ii < 0;
                    
                    clf
                    plot(clickExtractii);
                    hold on
                    plot(xEnvii)
                    plot([1,n],energyThrAbs*[1,1]);
                    if ~isempty(iStartAboveThrii)
                        plot(iStartAboveThrii*[1,1],max(xEnvii)*[-1,1],'-g')
                    end
                    if ~isempty(iEndAboveThrii)
                        plot(iEndAboveThrii*[1,1],max(xEnvii)*[-1,1],'-r')
                    end
                    stem(find(aboveThrii),energyThrAbs*ones(1,sum(aboveThrii)),'b');
                    plot(find(slopeRise),xEnvii(slopeRise),'^c');
                    plot(find(slopeFall),xEnvii(slopeFall),'vm');
                    grid
                    
                    keyboard
                end
                %}
                % END DEBUG
            end
            
            % evaluate
            pass = (val > shortDurThr_ms) & (val <= longDurThr_ms);
        end
        
        
        % LEGACY CRITERIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % assess_deltaE --------------------------------------------
        function [pass,val] = assess_deltaE(data,~)
        % Assesses change in energy of envelope over specific intervals.
        % This is an adaptation of the 1st half of "bw_elim_false_AMAR250".
        % At time of writing this is no longer a standard criterion, but 
        % the code is retained to allow replication of SBP's results.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTES: The original behaviour for this parameter is to disregard 
        % end times of clicks and assume that they are all at least 0.35ms 
        % long. This behaviour is preserved here for compatibility, but the
        % variable "useDur" can be edited to take click ends into account,
        % in which case all clicks that are too short to assess will get an
        % automatic pass.
        
            useDur = false; % original behaviour is "false"
        
            % set thresholds
            int1_ms = [0,0.1];
            int2_ms = [0.125,0.35];
        
            % define variables
            fs = data.fs;
            offset = data.offset;
            pos = data.pos;
            yFilt = data.yFilt;
            dur = data.dur;
            nClicks = size(pos,1);
            
            % get intervals in samples
            int1 = round(int1_ms*(fs/1000));
            int1(int1 == 0) = 1;
            int2 = round(int2_ms*(fs/1000));
            nClickSamplesMin = int2(end);
            
            % loop through each click
            deltaE = NaN(nClicks,1);
            for ii = 1:nClicks
                % check click duration if specified
                if useDur
                    nClickSamplesii = round(dur(ii)*(fs/1000));
                    if nClickSamplesii < nClickSamplesMin
                        deltaE(ii) = Inf;
                        continue % click too short
                    end
                end
                
                % compute envelope
                xClickii = yFilt(ii,:);
                posStartii = offset + 1;
                posEndii = offset + nClickSamplesMin;
                clickExtractii = xClickii(posStartii:posEndii);
                xEnvii = abs(hilbert(clickExtractii));

                % calculate whether energy increases or decreases after first 0.1 ms (decrease = dolphin, increase = beaked whale)
                % compare maximum in envelope of first 0.1 ms with maximum in envelope between 0.125 - 0.35 ms
                envMax1ii = max(xEnvii(int1(1):int1(2))); % max of click envelope in first 0.1 ms
                envMax2ii = max(xEnvii(int2(1):int2(2))); % max of click envelope in 0.125-0.35 ms interval
                deltaE(ii) = envMax2ii-envMax1ii;
            end
            
            % evaluate
            pass = deltaE > 0; % increase = beaked whale
            val = deltaE;
        end
        
        % assess_nSamples ------------------------------------------
        function [pass,val] = assess_nSamples(data,th)
        % The original "nSamples" parameter, now abandoned in favour of 
        % "slopeDur"
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            val = data.nSamples;
            pass = val >= th(1) & val <= th(2);
        end 
    end
end