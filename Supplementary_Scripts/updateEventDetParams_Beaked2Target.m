%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "updateEventDetParams_Beaked2Target" 
%   Written by Wilfried Beslin
%   Last Updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   Replaces the old beaked-whale-centric names of the criteria within 
%   EventDetParams files to more generic ones. Specifically, this function 
%   converts the "MinNumBeaked" and "MinPercentBeaked" criteria to
%   "MinNumTarget" and MinPercentTarget", respectively.
%
%   Use this function to update the event detection criteria spreadsheets
%   for compatibility with TWD version 1.3.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% INPUT - CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full path to root folder containing EventDetParams spreadsheets
%%% This will typically be a DetectionCriteria folder
dirPath = 'D:\TWD-1.3.0\_BWD\DetectionCriteria';

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


do_updateEventDetParams_Beaked2Target(dirPath);


%% ------------------------------------------------------------------------
function do_updateEventDetParams_Beaked2Target(dirPath)    
    % get all EventDetParams Excel files
    paramFilePaths = TWD_Common.Utilities.listFiles(dirPath, 'xlsx', 'Recursive',true, 'NameMustContain','EventDetParams');
    numFiles = numel(paramFilePaths);
    uniquePathStart = numel(dirPath) + 2; 
    
    % loop through each file and process it
    for ii = 1:numFiles
        filePath_ii = paramFilePaths{ii};
        fprintf('Processing file "%s"\n', filePath_ii(uniquePathStart:end))
        
        try
            % read file
            critNames_ii = readtable(filePath_ii, 'Sheet',1, 'Range','A:A');
            
            % update criteria names
            critNames_ii{:,1} = strrep(critNames_ii{:,1}, 'Beaked', 'Target');
            
            % save changes to the spreadsheet file
            writetable(critNames_ii, filePath_ii, 'Sheet',1, 'Range','A:A')
            
        catch ME
            warning('Failed to process file: %s', ME.message)
        end
    end
    
    disp('Done')
end