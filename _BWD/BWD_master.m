%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% "BWD_master"
%   Written by Wilfried Beslin
%   For TWD version 1.3
%   Revised 2023/04/19
%
%   Despcription:
%   Interface script for running WB's version of the beaked whale detector.
%   Edit the variables in this file as needed, then run to launch automatic 
%   beaked whale detection.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deployment name
depName = 'DEP';

% path to analysis folder
dirPath_analysis = 'D:\BWD_Results'; 

% path to WAV file folder
%%% For step 1: click compilation. Set empty to ignore this step.
dirPath_audio = 'D:\wav_files';

% name of output folder for event detection results
%%% For step 2: beaked whale detection. Set empty to ignore this step.
dirName_detResults = 'results01';

% name of detection protocol
%%% Used in Step 2 only
detProtocol = 'GeneralWithDiscriminators';

% data segmentation 
%%% best for continuous data. Set to duration(Inf,Inf,Inf) for duty-cycled.
%%% Must be a MATLAB duration object. Syntax is "duration(h,m,s)".
segDur = duration(Inf,Inf,Inf);

% FFT samples
nfft = 512;  

% Max clicks per MAT file
nMATClicksMax = 10000; 

% END CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% determine path to code directory
dirPath_root = mfilename('fullpath');
[dirPath_root,~,~] = fileparts(dirPath_root);

% run detector
BWD.run(...
    dirPath_root,...
    depName,...
    dirPath_analysis,...
    dirPath_audio,...
    dirName_detResults,...
    detProtocol,...
    segDur,...
    nfft,...
    nMATClicksMax);