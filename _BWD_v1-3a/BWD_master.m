%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% "BWD_master"
%   Written by Wilfried Beslin
%   For BWD version 1.3a
%   Revised 2022/10/13
%
%   Despcription:
%   Interface script for running WB's version of the beaked whale detector.
%   Edit the variables in this file as needed, then run to launch automatic 
%   beaked whale detection.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deployment name
depName = 'TEST';

% path to analysis folder
dirPath_analysis = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\Testing_BWD_1month'; 

% path to WAV file folder
%%% For step 1: click compilation. Set empty to ignore this step.
dirPath_audio = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\audio_1month';

% name of output folder for event detection results
%%% For step 2: beaked whale detection. Set empty to ignore this step.
dirName_detResults = 'results02';

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