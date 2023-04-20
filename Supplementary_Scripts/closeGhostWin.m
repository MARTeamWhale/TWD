%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "closeGhostWin"
%   Written by Wilfried Beslin
%   Last Updated Apr. 20, 2023, using MATLAB R2018b
%
%   Description:
%   Removes invisible windows that cannot be deleted with a standard call 
%   to the close() function. Such windows may remain after the BWD/SWD
%   validation tool encounters a hard crash; run this script to delete
%   them.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function closeGhostWin()
    w = findobj('Type','Figure','Visible','off');
    nw = numel(w);
    for ii = 1:nw
        delete(w(ii))
    end
    fprintf('Closed %d hidden window(s)\n',nw)
end