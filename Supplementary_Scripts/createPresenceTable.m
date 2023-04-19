%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "createPresenceTable"
%   Written by Wilfried Beslin
%   Last updated Oct. 13, 2022, using MATLAB R2018b
%
%   Description:
%   Reads in a "Validated" spreadsheet produced by the "identifySpecies"
%   script, interprets the Species ID codes, and produces a "Presence"
%   spreadsheet containing presence/absence/unsure scores for each species
%   where 1 = present, 0 = absent, and -1 = uncertain.
%
%   - Works using SpeciesCodes.xlsx spreadsheet in BWD/SWD code directory
%   - Only species whose code digit was used in the dataset will have a 
%   column
%   - 'Other' indicates that there were no beaked/sperm whales
%   - 'Uncertain' code can be used alone or can succeed code digits for 
%   individual species
%       -- If alone, then all species will be marked uncertain
%       -- If succeeding a species code, then that species will be marked
%       uncertain
%   - Unrecognized codes will result in undefined (NaN) scores for all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
% 21/10/2019
% - in future perhaps the detector module used (i.e. BWD or SWD) could be 
%   logged in the results, this way there would be no need to for the user 
%   to specify the folder here

function createPresenceTable()

    % INPUT - CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% full path to input "Validated" spreadsheet
    inFilePath = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\Testing_BWD_1-3\results01\TEST_Beaked_Validated.xlsx';
    %%% name of BWD or SWD folder, whichever is appropriate
    detectorDir = 'BWD_v1-3a';
    % END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % define other variables
    %%% path to code directory
    dirPath_root = mfilename('fullpath');
    [dirPath_root,~,~] = fileparts(dirPath_root);
    %%% species ID file path
    IDFilePath = fullfile(dirPath_root,detectorDir,'SpeciesCodes.xlsx');
    
    % make sure input file is a "Validated" spreadsheet 
    % (filename must contain "Validated")
    [dirPath,inFileName,fileExt] = fileparts(inFilePath);
    assert(contains(inFileName,'Validated'),'Input file must be a "Validated" spreadsheet ("Validated" must be in the file name)')
    
    % setup corresponding output file path 
    outFileName = strrep(inFileName,'Validated','Presence');
    outFilePath = fullfile(dirPath,[outFileName,fileExt]);
    
    % load validation spreadsheet
    inTable = readtable(inFilePath);
    
    % load species ID codes
    IDTable = readtable(IDFilePath);
    
    % get species ID codes
    IDNums = inTable.Species;
    
    % convert codes to presence/absence/uncertain
    presenceTable = readSpeciesIDs(IDNums,IDTable);
    
    % create output table
    outTable = [inTable(:,1:2),presenceTable];
    
    % save output table (overwrite if it already exists)
    if logical(exist(outFilePath,'file'))
        delete(outFilePath)
    end
    writetable(outTable,outFilePath)
end

%% readSpeciesIDs ---------------------------------------------------------
function presence = readSpeciesIDs(IDNums,IDTable)
% Converts a vector of species ID codes into a table of 
% presence/absence/unsure scores for each species.
% Columns are added only for species that have been mentioned (i.e. their
% code digits occur at least once somewhere)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % get logical indices of species codes, as well as 'Other' and
    % 'Uncertain' codes. Make sure the latter two actually exist.
    isOther = strcmpi(IDTable.Species,'Other');
    isUncertain = strcmpi(IDTable.Species,'Uncertain');
    isSp = ~isOther & ~isUncertain;
    nSp = sum(isSp);
    assert(sum(isOther)==1 && sum(isUncertain)==1, '"SpeciesCodes" spreadsheet must include codes for "Other" and "Uncertain"');
    
    % get list of species and numeric codes (converted to char)
    spList = IDTable.Species(isSp);
    spCodeStrList = cellstr(num2str(IDTable.Code(isSp)));
    codeNumStr_Sp = num2str(IDTable.Code(isSp))';
    codeNumStr_Other = num2str(IDTable.Code(isOther));
    codeNumStr_Uncertain = num2str(IDTable.Code(isUncertain));
    
    % define regular expression for finding all digits of interest that
    % may or may not be succeeded by the Uncertain digit
    expr = ['[',codeNumStr_Sp,'](',codeNumStr_Uncertain,')?'];
    
    % get number of events and initialize output
    n = numel(IDNums);
    spFieldsAll = strcat('Presence_',spList);
    presence = array2table(zeros(n,nSp),'VariableNames',spFieldsAll);
    allUncertain = false(n,1);
    
    % loop through each recording and get presence/absence
    for ii = 1:n
        
        % convert numeric code to string and read it
        IDStrii = num2str(IDNums(ii));
        
        % check code:
        % if "Other", all absent
        % if "Uncertain" only (or Other + Uncertain), all unknown
        % otherwise, determine presence/absence by species
        switch IDStrii
            case codeNumStr_Other
                % do nothing
                
            case {codeNumStr_Uncertain,[codeNumStr_Other,codeNumStr_Uncertain]}
                % mark all uncertain
                presence{ii,:} = -1;
                allUncertain(ii) = true;
                
            otherwise
                % interpret code
                regexMatchii = regexp(IDStrii,expr,'match');
                
                % make sure the number of characters in regex match is
                % equivalent to number of characters in the code. If not,
                % then this is an unrecognized code.
                nCharCodeii = numel(IDStrii);
                nCharRGMii = numel(strjoin(regexMatchii,''));
                goodCodeii = nCharCodeii == nCharRGMii;
                
                if goodCodeii
                    % translate code to presence/absence/unknown for each
                    % species
                    nSppCodesii = numel(regexMatchii);
                    for jj = 1:nSppCodesii
                        % get code part corresponding to one species
                        codejj = regexMatchii{jj};
                        
                        % get table fieldname for species
                        spCodejj = codejj(1);
                        fieldNamejj = spFieldsAll{strcmp(spCodejj,spCodeStrList)};
                        
                        % determine if species is present or uncertain
                        % based on number of digits in code part
                        if isscalar(codejj)
                            % only one digit -> present
                            scorejj = 1;
                        else
                            % two digits -> uncertain
                            scorejj = -1;
                        end
                        
                        % assign score to table
                        presence.(fieldNamejj)(ii) = scorejj;
                    end
                else
                    % unrecognized code -> set all to NaN
                    presence{ii,:} = NaN;
                end
        end
    end
    
    % find species that were not mentioned and remove them from table
    mentionMat = presence{:,:} ~= 0;
    mentionMat(allUncertain,:) = false;
    spMentioned = any(mentionMat,1);
    presence = presence(:,spMentioned);
end