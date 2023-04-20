%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function "blue"
%   Last updated presumably Sep. 23, 2006 by Mellisa Soldevilla
%   (Added to TWD by Wilfried Beslin on Dec. 3, 2019)
%
%   Note that this is a Triton file. Apart from the rebranded header and
%   removal of what I assume is an automatic timestamp and author tracker
%   (relating to "CVS"), this file is unmodified.
%
%   ORIGINAL HEADER:
%BLUE    Blue-tone color map.
%   BLUE(M) returns an M-by-3 matrix containing a "blue" colormap.
%   BLUE, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(blue)
%
%   See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG, 
%   COLORMAP, RGBPLOT.
%
%   9-21-06 mss Ripped from:
%   HOT.m
%   C. Moler, 8-17-88, 5-11-91, 8-19-92.
%   Copyright 1984-2002 The MathWorks, Inc. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bl = blue(m)

if nargin < 1, m = size(get(gcf,'colormap'),1); end
n = fix(1/3*m);

r = [zeros(2*n,1); (1:m-2*n)'/(m-2*n)];
g = [zeros(n,1); (1:n)'/n; ones(m-2*n,1)];
b = [(1:n)'/n; ones(m-n,1)];


bl = [r g b];