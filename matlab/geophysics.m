function [m] = geophysics(m)
%   GEOPHYSICS(M) returns an M-by-3 matrix containing the a geophysics 
%   appropriate colormap, a variant of JET(M). The colors begin with 
%   dark blue, range through shades of blue, cyan, green, yellow and red, 
%   and end with pale pink. GEOPHYSICS, by itself, is the same length as 
%   the current figure's colormap. If no figure exists, MATLAB uses the 
%   length of the default colormap.
%
%   See also JET, PARULA, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.

%   Copyright 1984-2015 The MathWorks, Inc.
%   Modified by James Wookey

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end


j = jet(64);

p = [1.0 96/255 244/255] ;

j(57:72,1) = 1 ;
j(57:72,2) = linspace(0,p(2),16) ;
j(57:72,3) = linspace(0,p(3),16) ;

j(73:88,1) = 1 ;
j(73:88,2) = linspace(p(2),1,16) ;
j(73:88,3) = linspace(p(3),1,16) ;

j = j(1:84,:) ;

m = [ interp1(linspace(1,m,84),j(:,1),1:m)', ...
      interp1(linspace(1,m,84),j(:,2),1:m)', ...
      interp1(linspace(1,m,84),j(:,3),1:m)'] ;

end

