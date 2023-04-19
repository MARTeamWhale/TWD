% script for closing invisible unclosable windows that may remain after 
% BWD/SWD encounters a hard crash.
% 2019/09/23

function closeGhostWin()
    w = findobj('Type','Figure','Visible','off');
    nw = numel(w);
    for ii = 1:nw
        delete(w(ii))
    end
    fprintf('Closed %d hidden window(s)\n',nw)
end