%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% "SWD_master"
%   Written by Wilfried Beslin
%   For TWD version 1.3
%   Revised 2023/04/19
%
%   Despcription:
%   Interface script for running experimental sperm whale detector.
%   Edit the variables in this file as needed, then run to launch automatic 
%   sperm whale detection.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% CHANGE AS NEEDED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deployment name
depName = 'TEST';

% path to analysis folder
dirPath_analysis = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\Testing_SWD_1-3';

% path to WAV file folder
%%% For step 1: click compilation. Set empty to ignore this step.
dirPath_audio = 'C:\Users\BeslinW\Desktop\Beakies\TWD_Testing\audio';

% name of output folder for event detection results
%%% For step 2: sperm whale detection. Set empty to ignore this step.
dirName_detResults = 'results03';

% detection parameter directory path
%%% Used in both steps 1 and 2
detProtocol = 'Sperm';

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
SWD.run(...
    dirPath_root,...
    depName,...
    dirPath_analysis,...
    dirPath_audio,...
    dirName_detResults,...
    detProtocol,...
    segDur,...
    nfft,...
    nMATClicksMax);