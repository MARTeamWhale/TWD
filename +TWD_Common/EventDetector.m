%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Class "EventDetector"
%   Written by Wilfried Beslin
%   Last updated Apr. 12, 2023, using MATLAB R2018b
%
%   Description:
%   Class for controlling event detection based on click discrimination 
%   results. Major uses include:
%   - defines the event detection criteria functions (e.g. MinNumTarget)
%   - manages properties of event criteria (i.e. loading from file,
%   storing, updating, retrieving)
%   - evaluates a dataset against the detection functions. In other
%   words, use this class to determine if a recording has target whales
%   present or not.
%
%   NOTE:
%   This version adds a new criterion, which requires click start times.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef EventDetector < handle
    %% PROPERTIES =========================================================
    properties (SetAccess = private)
        critTable
        missTol
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
            'MinNumTarget',     @TWD_Common.EventDetector.assess_MinNumTarget,...
            'MinPercentTarget', @TWD_Common.EventDetector.assess_MinPercentTarget,...
            'MaxMedianICI', @TWD_Common.EventDetector.assess_MaxMedianICI)
    end
    
    %% METHODS - PUBLIC ===================================================
    methods (Access = public)
        %% Constructor ----------------------------------------------------
        function obj = EventDetector(varargin)
        % Creates an "EventDetector" object.
        % Empty creation is allowed.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            narginchk(0,1)
            
            % initialize properties
            critNames = fieldnames(obj.critFuncs);
            nCrit = numel(critNames);
            critMat = zeros(nCrit,2);
            obj.critTable = array2table(critMat,'RowNames',critNames,'VariableNames',{'Threshold','UseCategory'});
            obj.missTol = double.empty(0,1);
            
            % load criteria if any were specified
            if nargin > 0
                critParamsFilePath = varargin{:};
                obj.loadCritParams(critParamsFilePath)
            end
        end
        
        %% loadCritParams -------------------------------------------------
        function loadCritParams(obj,critParamsFilePath)
        % Sets criteria properties from an Excel file 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % 1) read criteria table
            inCritTableRaw = readtable(critParamsFilePath,'ReadRowNames',true);
            
            % make sure the criteria names are consistent, remove the ones
            % that are unknown.
            critNames = fieldnames(TWD_Common.EventDetector.critFuncs);
            inCritNamesRaw = inCritTableRaw.Properties.RowNames;
            badCrit = ~ismember(inCritNamesRaw,critNames);
            badCritNames = inCritNamesRaw(badCrit);
            for ii = 1:sum(badCrit)
                warning('Ignoring unknown event detection criterion "%s"',badCritNames{ii})
            end
            inCritTable = inCritTableRaw(~badCrit,:);
            inCritNames = inCritTable.Properties.RowNames;
            
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
        % Threshold or UseCategory). Specify properties as Name-Value
        % pairs.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % parse input
            assert(ismember(critName,fieldnames(obj.critFuncs)),sprintf('"%s" is not a valid criterion',critName))
            oldVals = obj.critTable{critName,:};
            p = inputParser;
            p.addParameter('Threshold',oldVals(1))
            p.addParameter('UseCategory',oldVals(2))
            p.parse(varargin{:})
            
            % set new values
            newVals = [...
                p.Results.Threshold,...
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
        
        %% evalTargetEvent ------------------------------------------------
        function [hasTarget,nTargetClicks,eventPassTable,clickPassTable] = evalTargetEvent(obj,data,clickDiscrimObj)
        % Major function - assesses if a recording has target whales or not
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NOTE: This is somewhat of a legacy version. It takes a
        % ClickDiscriminator object as input, and thus does not support
        % multiple discriminators.
        
            % determine which criteria to test
            use = obj.critTable.UseCategory > 0;
            critNames = obj.critTable.Properties.RowNames(use);
            critCats = obj.critTable.UseCategory(use);
            nCrit = sum(use);
            cats = unique(critCats);
            
            % do click discrimination
            [targetClick,clickPassTable] = clickDiscrimObj.evalTargetClicks(data);  
            nTargetClicks = sum(targetClick);
            
            % initialize output
            eventPassTable = array2table(false(1,nCrit),'VariableNames',critNames);
            
            % loop through each criterion and evaluate it
            for ii = 1:nCrit
                critii = critNames{ii};
                fcnii = obj.critFuncs.(critii);
                thii = obj.critTable.Threshold(critii);
                passii = feval(fcnii,targetClick,data.pos(:,1),thii);
                
                eventPassTable.(critii) = passii;
            end
            
            % get summary
            %%% count criteria that have passed in each category
            hasTarget = true;
            for ii = cats
                iCritsii = critCats == ii;
                nCritsii = sum(iCritsii);
                passMatii = eventPassTable{:,iCritsii};
                nPassedii = sum(passMatii);
                passii = nPassedii >= (nCritsii - obj.missTol(ii));
                
                % update hasTarget
                hasTarget = hasTarget & passii;
            end
        end
        
        %% evalTargetEvent_logical ----------------------------------------
        function [hasTarget,nTargetClicks,eventPassTable] = evalTargetEvent_logical(obj,targetClick,clickStart)
        % Alternate event detection function - takes logical vector as 
        % input instead of data and a ClickDiscriminator. Thus, it does not
        % execute click discrimination - this must be done beforehand. The
        % advantage of this is that multiple discriminators could be used
        % as needed.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % determine which criteria to test
            use = obj.critTable.UseCategory > 0;
            critNames = obj.critTable.Properties.RowNames(use);
            critCats = obj.critTable.UseCategory(use);
            nCrit = sum(use);
            cats = unique(critCats);
            
            % get number of target clicks
            nTargetClicks = sum(targetClick);
            
            % initialize output
            eventPassTable = array2table(false(1,nCrit),'VariableNames',critNames);
            
            % loop through each criterion and evaluate it
            for ii = 1:nCrit
                critii = critNames{ii};
                fcnii = obj.critFuncs.(critii);
                thii = obj.critTable.Threshold(critii);
                passii = feval(fcnii,targetClick,clickStart,thii);
                
                eventPassTable.(critii) = passii;
            end
            
            % get summary
            %%% count criteria that have passed in each category
            hasTarget = true;
            for ii = cats
                iCritsii = critCats == ii;
                nCritsii = sum(iCritsii);
                passMatii = eventPassTable{:,iCritsii};
                nPassedii = sum(passMatii);
                passii = nPassedii >= (nCritsii - obj.missTol(ii));
                
                % update hasTarget
                hasTarget = hasTarget & passii;
            end
        end
    end
    
    %% METHODS - PUBLIC STATIC ============================================
    % These are the criteria functions. 
    % If you want to add a new criterion, implement it here. Make sure the
    % function names match the criteria names, otherwise the public
    % evaluation function won't recognize them.
    methods (Static)
        
        % assess_MinNumTarget ---------------------------------------------
        function pass = assess_MinNumTarget(clickPassed,~,th)
            nTarget = sum(clickPassed);
            pass = nTarget >= th;
        end

        % assess_MinPercentTarget -----------------------------------------
        function pass = assess_MinPercentTarget(clickPassed,~,th)
            nTarget = sum(clickPassed);
            nTotal = numel(clickPassed);
            percentPassed = (nTarget/nTotal)*100;
            pass = percentPassed >= th;
        end
        
        % assess_MaxMedianICI ---------------------------------------------
        function pass = assess_MaxMedianICI(clickPassed,clickStart,th)
            ICI = diff(clickStart(clickPassed));
            pass = median(ICI) <= seconds(th);
        end
        
    end
end