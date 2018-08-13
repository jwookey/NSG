function [varargout]=plot_EM31(fn,varargin)
%
%  Function to load and plot a (converted) EM31 dataset (xyz file created by)
%  TrackMaker32. 
%
%  USAGE: 
%
%  plot_EM31(fn)
%     
%     Plot the raw Q and I components in a specified xyz file ('fn'). 
%     Alternatively fn can be a structure as output by this code (see 
%     later). This can be useful for (for example) performing corrections
%     on the data before gridding, interpolating or ARC file output. 
%
%  plot_EM31(...,'grid',R)
%     
%     Additionally produce a cap averaged surface plot of the EM31 data. R
%     is the cap radius (useful numbers are between 2-10m). Where no data
%     is within the averaging radius, a NaN is put in the surface grid.
%
%  plot_EM31(...,'arc')
%
%     Produce a ARC raster files of the gridded EM31 data (conductivity and 
%     in-phase). This does nothing if gridding is not specified. 
%
%  plot_EM31(...,'line',[x0 x1 y0 y1],ds,R)
%     
%     Extract a line of data from x0,y0 to x1,y1 at a spacing ds. Each
%     value is a gaussian-weighted average of FWHM R, of all recorded 
%     values of conductivity within radius R*5, or NaN otherwise. If NaN is 
%     specified for R, the closest point (at whatever distance) is used. 
%
%  plot_EM31(...,'coord',id_str)
%
%     Transform the default (UTM) coordinates into another coordinate scheme.
%     id_str specifies the transform to use (see source code for available
%     schemes).
%
%  plot_EM31(...,'arc')
%
%     Produce a ARC raster files of the gridded EM31 data (conductivity and 
%     in-phase). This does nothing if gridding is not specified.  
%
%  plot_EM31(...,'noplot')
%
%     Suppress plotting of the scatter plots. Line or grid plots (if requested)
%     are still plotted. 
%
%  [S] = plot_EM31(...)
%
%     Return whatever has been calculated (e.g., line and/or gridded data) in
%     the structure S. This has the following fields:
%     by default:
%        S.fn    - data file name
%        S.x     - data x-coordinates
%        S.y     - data y-coordinates
%        S.q     - quadrature (mS/m)
%        S.in    - in-phase (ppt)
%
%     in addition, if gridding is specifed:
%        S.GX    - interpolated grid x-matrix
%        S.GY    - interpolated grid y-matrix
%        S.GQ    - interpolated grid quadrature (mS/m)
%        Q.GI    - interpolated grid in-phase (ppt)
%         
%     in addition, if line interpolation is specifed:
%        S.ldist - interpolated line distance vector
%        S.lazi  - azimuth of the line
%        S.lx    - line x-vector
%        S.ly    - line y-vector
%        S.lq    - interpolated quadrature (mS/m)
%        S.li    - interpolated in-phase (ppt)
%

% Copyright (c) 2016-2018 James Wookey
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
%    * Redistributions of source code must retain the above copyright notice, 
%      this list of conditions and the following disclaimer.
%    * Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the 
%      documentation and/or other materials provided with the distribution.
%    * Neither the name of the University of Bristol nor the names of its 
%      contributors may be used to endorse or promote products derived from this 
%      software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


%% Grid the conductivity data using a cap average approach
sr = NaN ; % averaging radius

new_coord_id = 'utm' ;
iarc = 0 ;
iline = 0 ;

% suppress plots
iplot = 1 ;

% process the optional input arguments
iarg = 1 ;
while iarg <= (length(varargin))
   switch lower(varargin{iarg})
   case 'grid'
      sr = varargin{iarg+1} ;
      iarg = iarg + 2 ;
   case 'coord'
      new_coord_id = varargin{iarg+1} ;
      iarg = iarg + 2 ;
   case 'arc'
      iarc = 1 ;
      iarg = iarg + 1 ;
   case 'noplot'
      iplot = 0 ;
      iarg = iarg + 1 ;
   case 'line'
      iline = 1 ;
      pv = varargin{iarg+1} ;
      lspacing = varargin{iarg+2} ;
      laveR = varargin{iarg+3} ;
      iarg = iarg + 4 ;

   otherwise 
      error(['Unknown option: ' varargin{iarg}]) ;   
   end   
end 

% process the optional output arguments
if nargout>1
   error('Too many output arguments') ;
elseif nargout==1
   S.filename = fn ;
end


% load the data file or structure

if isstruct(fn) 
   Easting = fn.x ;
   Northing = fn.y ;
   q = fn.q ;
   in = fn.in ;
else % assume its a file
   [Easting,Northing,Elev,q,in,timeStr] = load_EM31_ascii(fn) ;
end

switch lower(new_coord_id)
case 'osgb_dh'
% convert UTM to OSGB (for Dartington Hall locality)
   x = Easting-170894 ;
   y = Northing-5526388 ;

case 'osgb_downs'
% convert UTM to OSGB (the Downs locality)
%        UTM        OSGB       Diff
% Eing   526202     356845     169357 
% Ning   5702516    175100     5527416     

   x = Easting-169357 ;
   y = Northing-5527416 ;

case 'osgb_penwylt'
% convert UTM to OSGB (the Downs locality)
%        UTM        OSGB       Diff
% Eing   454424     285604     168820 
% Ning   5741875    215441     5526434     
%454424 5741875
%285604 215441

   x = Easting-169357 ;
   y = Northing-5527416 ;

case 'utm'
% keep in UTM
   x = Easting ;
   y = Northing ;
otherwise
   warning('Specified coordinate scheme not supported, output in UTM') ;
   x = Easting ;
   y = Northing ;
end

if nargout
   S.x = x ;
   S.y = y ;
   S.q = q ;
   S.in = in;
end

if iplot

   %% plot the conductivity
   figure
   scatter(x,y,50,q,'filled')
   axis equal
   hold on
   title('Conductivity (mS/m)')
   xlabel('East (m)')
   ylabel('North (m)')
   colormap(jet)
   colorbar
   daspect([1 1 1])

   if iline
      hold on
      plot([pv(1) pv(3)],[pv(2) pv(4)],'w-','LineWidth',3) ;
      plot([pv(1) pv(3)],[pv(2) pv(4)],'k-','LineWidth',2) ;
   end
   
   %% Plot in-phase component
   figure
   scatter(x,y,50,in,'filled')
   axis equal
   hold on
   title('In-phase (ppt)')
   xlabel('East (m)')
   ylabel('North (m)')
   colormap(jet)
   colorbar
   daspect([1 1 1])

   if iline
      hold on
      plot([pv(1) pv(3)],[pv(2) pv(4)],'w-','LineWidth',3) ;
      plot([pv(1) pv(3)],[pv(2) pv(4)],'k-','LineWidth',2) ;
   end

end

% create line
if iline
   % distance vector
   ldist = [0:lspacing:sqrt( (pv(3)-pv(1)).^2 +(pv(4)-pv(2)).^2 )] ;
   lazi = atan2d(pv(3)-pv(1),pv(4)-pv(2)) ;
   % lx, ly
   lx = pv(1)+ldist.*sind(lazi) ;
   ly = pv(2)+ldist.*cosd(lazi) ;
   
   % interpolate data to line
   lq = lx.*NaN ;
   li = lx.*NaN ;
   

   % generate weights vector
   xwgt = 0:1/laveR:laveR*5 ;
   ywgt = gaussian(xwgt,laveR) ;

   for i=1:length(ldist) ;
      dist=sqrt( (lx(i)-x).^2 + (ly(i)-y).^2 ) ;
      % interpolate weights
      
      % list of close points
      if isnan(laveR) 
         % use closest
         [~,imin] = min(dist) ;
         lq(i)=mean(q(imin)) ;
         li(i)=mean(in(imin)) ; 
      else   
         index = find(dist<laveR*5) ;
         if ~isempty(index)
            wgt = interp1(xwgt,ywgt,dist(index),'pchip',0) ;
            %lq(i)=mean(q(index)) ;
            %li(i)=mean(in(index)) ;   
            lq(i)=sum(q(index).*wgt)./sum(wgt) ;
            li(i)=sum(in(index).*wgt)./sum(wgt) ;
         end
      end
   end

   %% Plot line data
   figure
   yyaxis left
   plot(ldist,lq,'b-') ;
   ylabel('Conductivity (mS/m)') ;
   xlabel('Distance (m)')
   yyaxis right ;
   plot(ldist,li,'r-');
   title('Extracted line data')
   ylabel('In-phase (ppt)')
   axis tight
   if nargout
      S.ldist = ldist ;
      S.lazi = lazi ;
      S.lx = lx ;
      S.ly = ly ;
      S.lq = lq ;
      S.li = li ;
   end

end


if isnan(sr)
   if nargout
      varargout = {S} ;
   end
   return % no need to grid.
end


% create grid
gx = min(x)-sr*2:sr/2:max(x)+sr*2; % x/e grid definition vector
gy = min(y)-sr*2:sr/2:max(y)+sr*2; % y/n grid definition vector
[GX,GY] = meshgrid(gx,gy);
GQ = GX.*NaN ;
GI = GX.*NaN ;

ngx=length(gx) ;
ngy=length(gy) ;

% interpolate raw data to grid
for ix=1:ngx
   for iy=1:ngy
      dist=sqrt( (gx(ix)-x).^2 + (gy(iy)-y).^2 ) ;
      % list of close points
      index = find(dist<sr) ;
      if ~isempty(index)
         GQ(iy,ix)=mean(q(index)) ;
         GI(iy,ix)=mean(in(index)) ;
      end

   end
end

%% Plot the gridded data

figure
pcolor(GX,GY,GQ)
axis equal
hold on
colormap(jet);
%plot(x,y,'.k')
shading flat
colorbar
title('Gridded conductivity (mS/m)')
xlabel('East (m)')
ylabel('North (m)')
daspect([1 1 1])

figure
pcolor(GX,GY,GI)
axis equal
hold on
colormap(jet);
%plot(x,y,'.k')
shading flat
colorbar
title('Gridded InPhase (ppt)')
xlabel('East (m)')
ylabel('North (m)')
daspect([1 1 1])

if iarc
%% Export the gridded data to an ASCII raster file in ARC format
   arcgridwrite('EM_Cond.asc',GX,GY,GQ,'nobar') ;
   arcgridwrite('EM_InPh.asc',GX,GY,GI,'nobar') ;
end

if nargout
   S.GX = GX ;
   S.GY = GY ;
   S.GQ = GQ ;
   Q.GI = GI ;
end

end

% LOAD_EM31_ASCII - load an EM31 datafile as produced by trackmaker31
% This is for GPS-enabled data. It assumes that the elevation from the 
% GPS has been included.

function [Eing,Ning,Elev,Quad,InPh,timeStr]=load_EM31_ascii(filename)

%% Initialize variables.
startRow = 3;
formatSpec = '%15f%15f%11f%11f%10f%4f%4f%7f%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '',...
                     'WhiteSpace', '', 'HeaderLines' ,startRow-1,...
                     'ReturnOnError', false);

%% Remove white space around all cell columns.
dataArray{9} = strtrim(dataArray{9});

%% Close the text file.
fclose(fileID);

%% Allocate imported array to column variable names
Eing = dataArray{:, 1};
Ning = dataArray{:, 2};
Quad = dataArray{:, 3};
InPh = dataArray{:, 4};
Elev = dataArray{:, 5};
qual = dataArray{:, 6};
nsat = dataArray{:, 7};
hdop = dataArray{:, 8};
timeStr = dataArray{:, 9};

end

function arcgridwrite(fileName,X,Y,Z,varargin)
%ARCGRIDWRITE- Write gridded data set in Arc ASCII Grid Format
%
%   arcgridwrite(fileName,X,Y,Z)- converts data in a matlab
%   grid (as produced by eg. meshgrid and griddata) into a text file
%   in Arc ASCII Grid Format. The file can also be automatically
%   converted to raster on PC's with ARCINFO installed.
%
%   INPUTS
%       fileName:  output filename including extension
%       X:         X coordinates (m x n) 
%       Y:         Y coordinates (m x n) 
%       Z:         gridded data (m x n)
%
%   SYNTAX AND OPTIONS
%       arcgridwrite('D:\tools\bathyGrid.asc',X,Y,Z)
%
%       arcgridwrite('D:\tools\bathyGrid.asc',X,Y,Z,...
%           'precision',5) - changes default floating point output
%           from 3 to 5 decimal places.
%
%       arcgridwrite('D:\tools\bathyGrid.asc',X,Y,Z,'nobar') 
%           suppress the progress bar. 
%
%   EXAMPLE - create a raster grid of the peaks function
%       [X,Y,Z]=peaks(100);
%       arcgridwrite('peaksArc.asc',X,Y,Z,'precision',3,...
%           'convert')
%
%   NOTES 
%   1) Because the Arc ASCII format has only one parameter for cell size,
%   both X and Y must have the same, non-varying grid spacing.  
%
% A.Stevens @ USGS 7/18/2007
% astevens@usgs.gov
% 
% Adapted by James Wookey, University of Bristol, 31/10/2016
% glxjw@bris.ac.uk

%check number of inputs
if nargin < 4
    error('Not enough input arguments');
end

[mz,nz,zz]=size(Z);
[mx]=size(X);
[my]=size(Y);

minX=min(X(:));
minY=min(Y(:));
dx=abs(diff(X(1,:)));
dy=abs(diff(Y(:,1)));

maxDiff=0.01; %threshold for varying dx and dy.  increase or
%decrease this parameter if necessary.  

%check input dimensions and make sure x and y 
%spacing are the same. 
if any(diff(dx)>maxDiff) || any(diff(dy)>maxDiff)
     error('X- and Y- grid spacing should be non-varying');
else
    dx=dx(1); 
    dy=dy(1);
end

if ischar(fileName)~=1
    error('First input must be a string')
elseif zz>1
    error('Z must be 2-dimensional');
elseif mx~=mz
     error('X, Y and Z should be same size');
elseif my~=mz
    error('X, Y and Z should be same size');
elseif abs(dx-dy)>maxDiff;
    error('X- and Y- grid spacing should be equal');   
end


convert=0;
nobar=0 ;
dc=3; %default number of decimals to output
if isempty(varargin)~=1
    [m1,n1]=size(varargin);
    opts={'precision';'nobar'};

    for i=1:n1;
        indi=strcmpi(varargin{i},opts);
        ind=find(indi==1);
        if isempty(ind)~=1
            switch ind
                case 1
                    dc=varargin{i+1};
                case 2
                    nobar=1;
            end
        else
        end
    end
end


%replace NaNs with NODATA value 
Z(isnan(Z))=-9999;

%Z=flipud(Z); JW -> this results in an inverted grid in ArcMap

%define precision of output file
if isnumeric(dc)
    dc=['%.',sprintf('%d',dc),'f'];
elseif isnumeric(dc) && dc==0
    dc=['%.',sprintf('%d',dc),'d'];
end

fid=fopen(fileName,'wt');

%write header
fprintf(fid,'%s\t','ncols');
fprintf(fid,'%d\n',nz);
fprintf(fid,'%s\t','nrows');
fprintf(fid,'%d\n',mz);
fprintf(fid,'%s\t','xllcorner');
fprintf(fid,[dc,'\n'],minX);
fprintf(fid,'%s\t','yllcorner');
fprintf(fid,[dc,'\n'],minY);
fprintf(fid,'%s\t','cellsize');
fprintf(fid,[dc,'\n'],dx);
fprintf(fid,'%s\t','NODATA_value');
fprintf(fid,[dc,'\n'],-9999);

%configure filename and path
[pname,fname,ext] = fileparts(fileName);
if isempty(pname)==1;
    pname=pwd;
end

if strcmpi(pname(end),filesep)~=1
    pname=[pname,filesep];
end

%wait bar
if ~nobar
h = waitbar(0,['Writing file: ',fname,ext...
    ', Please wait...']);
end
%write data
test=repmat([dc,'\t'],1,nz); 
test(end-1:end)='\n'; 
%write data 
for i=mz:-1:1 
   fprintf(fid,test,Z(i,:)); 
   if ~nobar
       %update waitbar 
       waitbar(i/mz,h,['Writing file: ',[fname,ext],... 
           sprintf(' %d%% complete...',round(i/mz*100))]) 
   end
end 


if ~nobar,close(h),end

fclose(fid);

end

function y = gaussian(x,FWHM) ;

   sig = FWHM ./ sqrt(2.*log(2)) ;
   mu = 0.0 ; % no shift
   
   x1 = -(x-mu).^2./(2.*sig.^2) ;

   y = 1./(sig.*sqrt(2.*pi)) .* exp(x1) ;

% calculate the maximum value (at zero) ;
   x0 = -(0-mu).^2./(2.*sig.^2) ;
   ymax = 1./(sig.*sqrt(2.*pi)) .* exp(x0) ;

% normalised
   y=y./ymax ;

end
