%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "validate"
%   Written by WB
%   Last updated Oct. 13, 2022, using MATLAB R2018b
%
%   Description:
%   Prepares and launches the interactive event validation/species ID
%   process for BWD. Includes creation of ClickDiscriminator objects and
%   loading of reference spectra.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function validate(...
    dirPath_root,...
    dirPath_analysis,...
    dirName_detResults,...
    targetName,...
    eventMergeOpt,...
    eventMergeVal,...
    nClicksMax,...
    fRange,...
    monID)

    % set subfolder paths
    dirPath_mat = fullfile(dirPath_analysis,'mat');
    dirPath_out = fullfile(dirPath_analysis,dirName_detResults);
    dirPath_refSpec = fullfile(dirPath_root,'ReferenceSpectra');

    % get path to 1st MAT file, and extract deployment name, Fs, and nfft
    matNames = TWD_Common.Utilities.getFileNames(dirPath_mat,'mat');
    matName1 = matNames{1};
    %%% depName
    depExpr = '.*(?=_\d{8}_\d{6})'; % regular expression for extracting 1st part of MAT file name, i.e. deployment
    depName = regexp(matName1,depExpr,'match');
    assert(isscalar(depName),'MAT file names not recognized')
    depName = depName{1};
    %%% Fs and nfft
    data1 = load(fullfile(dirPath_mat,matName1),'fs','specClick');
    Fs = data1.fs;
    nfft = 2*size(data1.specClick,2);
    clear data1
    
    % determine detection protocol by reading text file
    detProtocol = TWD_Common.readProtocolTextFile(dirPath_out);
    dirPath_detCrit = fullfile(dirPath_root,'DetectionCriteria',detProtocol);
    
    % create click discriminator objects
    clickDiscrimData = createDiscriminators(dirPath_detCrit,targetName);
    
    % load reference spectra
    refSpecData = loadRefSpectra(dirPath_refSpec,Fs,nfft);

    % run validation app
    TWD_Common.runValidation(...
        'BWD',...
        dirPath_root,...
        dirPath_analysis,...
        dirName_detResults,...
        depName,...
        targetName,...
        eventMergeOpt,...
        eventMergeVal,...
        nClicksMax,...
        fRange,...
        monID,...
        clickDiscrimData,...
        refSpecData,...
        Fs)
end

%% createDiscriminators ---------------------------------------------------
function clickDiscrimData = createDiscriminators(dirPath_detCrit,targetName)
% Creates one or more ClickDiscriminator objects for event validation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % define file naming conventions
    clickCritFileFlag = 'ClickDiscrimParams_Validation';
    defaultsFileFlag = 'ValidationDefaults';
    
    % define regex to extract suffix from filename if there is one
    % - corresponds to discriminator name
    discrimNameExpr = ['(?<=',clickCritFileFlag,'_)\w*(?=\.xlsx)'];
    
    % read all spreadsheet names
    dirPath_target = fullfile(dirPath_detCrit,targetName);
    fileNames = TWD_Common.Utilities.getFileNames(dirPath_target,'xlsx');
    
    % get validation discriminator files
    iClickDiscrim = find(contains(fileNames,clickCritFileFlag));
    assert(isvector(iClickDiscrim),'In "%s":\nExpected at least one Excel file containing the name "%s"',dirPath_target,clickCritFileFlag)
    nClickDiscrim = numel(iClickDiscrim);
    discrimFileNames = fileNames(iClickDiscrim);
    
    % get default settings file if there is one
    iDefaults = find(contains(fileNames,defaultsFileFlag));
    if isscalar(iDefaults)
        try
            % load defaults file
            defaultsFilePath = fullfile(dirPath_target,fileNames{iDefaults});
            defaultsData = readtable(defaultsFilePath);
            
            % initialize
            discrimDefaults = false(nClickDiscrim,3);
            defaultOptions = {'Positive','Negative','Disabled'};
            
            % match discriminator filenames with table entries, and set
            % default settings
            for ii = 1:nClickDiscrim
                [~,discrimFileNameii,~] = fileparts(discrimFileNames{ii});
                
                rowii = strcmp(discrimFileNameii,defaultsData.File);
                colii = strcmpi(defaultsData.DefaultSetting{rowii},defaultOptions);
                
                discrimDefaults(rowii,colii) = true;
            end
            
            % make sure the defaults are valid
            assert(all(sum(discrimDefaults,2) == 1))
            validDefaults = true;
        catch
            warning('Invalid discrimination defaults, setting all to "Positive"')
            validDefaults = false;
        end
    else
        validDefaults = false;
    end
    if ~validDefaults
        discrimDefaults = repmat([true,false,false],nClickDiscrim,1);
    end
    discrimDefaults = num2cell(discrimDefaults);
    
    % initialize output
    discrimObjCell = cell(nClickDiscrim,1);
    discrimNames = cell(nClickDiscrim,1);
    
    % loop through each file and create discriminators
    for ii = 1:nClickDiscrim
        clickDiscrimFileii = discrimFileNames{ii};
        clickDiscrimFilePathii = fullfile(dirPath_target,clickDiscrimFileii);
        clickDiscrimii = TWD_Common.ClickDiscriminator(clickDiscrimFilePathii);
        discrimObjCell{ii} = clickDiscrimii;
        
        % get name
        regexOutii = regexp(clickDiscrimFileii,discrimNameExpr,'match');
        if isscalar(regexOutii)
            discrimNameii = regexOutii{:};
        else
            discrimNameii = targetName;
        end
        discrimNames{ii} = discrimNameii;
    end
    discrimObj = vertcat(discrimObjCell{:});
    
    % save output to struct
    clickDiscrimData = struct();
    clickDiscrimData.Obj = discrimObj;
    clickDiscrimData.Names = discrimNames;
    clickDiscrimData.Defaults = discrimDefaults;
end

%% loadRefSpectra ---------------------------------------------------------
function refSpecData = loadRefSpectra(dirPath_refSpec,Fs,nfft)
% loads reference spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEV NOTES:
%   - I think this desperately needs a better system

    % load MATfiles
    load(fullfile(dirPath_refSpec,'meanSpecMe'))
    load(fullfile(dirPath_refSpec,'meanSpecZc'))
    load(fullfile(dirPath_refSpec,'meanSpecMd'))
    load(fullfile(dirPath_refSpec,'meanSpecBW31'))
    load(fullfile(dirPath_refSpec,'meanSpecBW38'))
    load(fullfile(dirPath_refSpec,'meanSpecMb_Clarke2019'))
    load(fullfile(dirPath_refSpec,'meanSpecHa_Gully'))
    load(fullfile(dirPath_refSpec,'Mm_template_2016-2017_NEFSC_towed_array')) % True's template from NEFSC towed array data

    % set sampling rates for equipment; used for plotting mean spectra.
    %%% Maybe it would be best if that information was contained in the spectra 
    %%% MAT file(s)? For the reference spectra anyway.
    fAMAR=0:(Fs/2000)/(nfft/2-1):Fs/2000; % for current AMAR dataset

    fsHARP = 200000; % sampling rate for HARP reference spectra
    fsAMARref = 187500; % AMAR sampling rate for Ha reference spectrum
    fsARRAY = 192000; % array sampling rate for Mm reference spectrum

    N2=512; %%% watch out for this

    fHARP=0:(fsHARP/2000)/(N2/2-1):fsHARP/2000;
    fAMARref=0:(fsAMARref/2000)/(N2/2-1):fsAMARref/2000;
    fARRAY=0:(fsARRAY/2000)/(N2/2-1):fsARRAY/2000;

    %normalize min 10<f<35 (pos 26:91); max 20<f<80 (pos 52:205) for FFT 512
    %%% some work to be done here? Could this also be incorporated in mean
    %%% spectra MATs?
    meanSpecMe = meanSpecMe - min(meanSpecMe(26:91));
    meanSpecMe = meanSpecMe/max(meanSpecMe(52:205));
    meanSpecZc = meanSpecZc - min(meanSpecZc(26:91));
    meanSpecZc = meanSpecZc/max(meanSpecZc(52:205));
    meanSpecMd = meanSpecMd - min(meanSpecMd(26:91));
    meanSpecMd = meanSpecMd/max(meanSpecMd(52:205));
    meanSpecBW31 = meanSpecBW31 - min(meanSpecBW31(26:91));
    meanSpecBW31 = meanSpecBW31/max(meanSpecBW31(52:205));
    meanSpecBW38 = meanSpecBW38 - min(meanSpecBW38(26:91));
    meanSpecBW38 = meanSpecBW38/max(meanSpecBW38(52:205));
    meanSpecHa = meanSpecHa - min(meanSpecHa(26:91));
    meanSpecHa = meanSpecHa/max(meanSpecHa(52:205));
    meanSpecMm = Mm_normalized_meanSpec; %already normalized
    meanSpecMb = meanSpecMb - min(meanSpecMb(26:91));
    meanSpecMb = meanSpecMb/max(meanSpecMb(52:205));
    
    % structify!
    f = struct(...
        'fAMAR',fAMAR,...
        'fAMARref',fAMARref,...
        'fHARP',fHARP,...
        'fARRAY',fARRAY);
    refSpec = struct(...
        'meanSpecMe',meanSpecMe,...
        'meanSpecZc',meanSpecZc,...
        'meanSpecMd',meanSpecMd,...
        'meanSpecBW31',meanSpecBW31,...
        'meanSpecBW38',meanSpecBW38,...
        'meanSpecHa',meanSpecHa,...
        'meanSpecMm',meanSpecMm,...
        'meanSpecMb',meanSpecMb);
    
    % create final output
    refSpecData = struct();
    refSpecData.f = f;
    refSpecData.spec = refSpec;
end
