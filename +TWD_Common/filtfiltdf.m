%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function "filtfiltdf"
%   Written by WB
%   Last updated Sep. 19, 2019, using MATLAB R2018b
%
%   Description:
%   This utility function is intended to replace the native "filtfilt" 
%   method for digitalFilter objects.
%   - I've discovered that using filtfilt with digitalFilter objects as input
%   may result in numerical instabilities. A better approach is to first
%   extract coefficients from the digital filter, and use those as input.
%   - The best coefficients to use are the zero-pole-gain representation 
%   [z,p,k], which must then be converted to second-order-sections [sos,g] 
%   in order to use filtfilt. 
%   - It's also possible to extract the common transfer function 
%   coefficients [b,a] and use those with filtfilt. However this approach
%   is known to be less stable and Mathworks even discourages it: In the
%   documentaion for the "butter" function, it is written:
%       "In general, use the [z,p,k] syntax to design IIR filters. To 
%       analyze or implement your filter, you can then use the [z,p,k] 
%       output with zp2sos. If you design the filter using the [b,a] 
%       syntax, you might encounter numerical problems. These problems are 
%       due to round-off errors and can occur for n as low as 4."
%   - Interestingly, despite the known instabilities of [b,a], it seems 
%   that using these coefficients is even more reliable than digitalFilter
%   directly. So I don't know what filtfilt(d,x) is doing to be so
%   terrible...
%
%   WARNING: still haven't fully solved filter instabilities.
%   There can be issues with SOS as well. In some cases e.g. when trying 
%   highpass FIR and some IIR filters, using SOS (either by conversion 
%   from transfer function or zero-pole-gain coefficients) results in 
%   wonky waveform with crazy amplitude. Observations:
%   - there are two different SOS matrices, based on output arguments:
%       [sosg,g] = zpk2sos(z,p,k)
%           and
%       sos = zpk2sos(z,p,k);
%   These two sos matrices are not the same. g seems to be some sort of
%   weight, and if g is not called, this weight is integrated into the sos
%   matrix. (try fvtool(sos) and fvtool(sosg))
%   - filtfilt only accepts (sosg,g). Does not take sos alone. However, it
%   still results in wonky scales sometimes.
%   - sosfilt is the unidirectional option for filtering with sos, but 
%   input is opposite to filtfilt: this one does NOT take g at all. To get 
%   correct scale, do:
%       xFilt = g*sosfilt(sosg,x) OR xFilt = sosfilt(sos,x)
%   - generally it seems like filtfilt does not work well with sos. If
%   doing the following:
%       xFilt = sosfilt(sos,x);
%       xFilt = fliplr(sosfilt(sos,fliplr(xFilt)));
%   This results in a more correct-looking waveform with no phase delay, 
%   (more faithful to TF version), but the last few points will be off.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xFilt = filtfiltdf(d,x)

    % filter with second-order sections 
    % (Mathworks recommended method, but doesn't always work)
    %[z,p,k] = zpk(d); % zero-pole-gain
    %[sos,g] = zp2sos(z,p,k);
    %xFilt = filtfilt(sos,g,x);
    
    % filter with transfer function coefficients
    % (supposedly less stable than SOS, but sometimes works better)
    [b,a] = tf(d); 
    xFilt = filtfilt(b,a,x);
end

% steps/notes for recreating filtfilt with SOS
%   - get filter order with:
%       ord = filtord(sos);
%   - it would seem that all SOS matrices have dimensions k-by-6
%       - each row corresponds to tf coefficients for different sections...
%       - 1st 3 columns, i.e. sos(:,1:3) are 'b' coefficients
%       - last 3 columns, i.e. sos(:,4:6) are 'a' coefficients
%   - in filtfilt lines 140-148, it seems SOS matrix is transformed
%   according to g weights. It seems the calculation, generally, is:
%       for ii = 1:numel(g)
%           sos(ii,1:3) = g(ii)*sos(ii,1:3)
%       end
%   This implies that g is normally a vector where length equals number of
%   rows in SOS
%   - lines 160-168 is where "initial conditions" are obtained from SOS - 
%   important lines
%       - for each section (SOS row), compute the following:
%           rhs  = b(2:3) - b(1)*a(2:3))
%           lhs = eye(2) - [-a(2:3),[1;0]]
%           zi = lhs / rhs
%       Note that b and a in this case are each 3-element column vecs.
%   - initial conditions zi are specified as the last argument in the
%   'filter' function. Unfortunately, 'sosfilt' does not take zi.
%
%   - I give up. filtfilt is rather complex and sosfilt uses
%   MEX files, so can't easily tell how it works.