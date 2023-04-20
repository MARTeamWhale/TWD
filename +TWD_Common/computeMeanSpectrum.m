%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function "computeMeanSpectrum"
%   Written by WB, based on code written by SBP and JES.
%   Last updated Dec. 28, 2018, using MATLAB R2016b
%
%   Computes the average power spectrum given a matrix of log-transformed 
%   spectra. Rows = clicks, columns = magnitudes (log-transformed).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meanSpec = computeMeanSpectrum(specs)

    % transform spectrum to linear scale
    specsTrans = 10.^(specs/20);
    
    % take mean
    meanSpecTrans = mean(specsTrans,1);
    
    % revert to log scale
    meanSpec = 20*log10(meanSpecTrans);
end