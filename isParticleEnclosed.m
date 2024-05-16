function[isEnclosed]=isParticleEnclosed(x,y,xc,yc)
%isParticleEnclosed  True for particles enclosed within a given contour.
%
%   ISENCLOSED = isParticleEnclosed(XA,YA,XC,YC), where XA and YA give
%   article positions and XC and YC determine a set of contours, returns 
%   a boolean variable, ISENCLOSED, that is true for enclosed particles. 
% 
%   XA and YA are 2D arrays of particle positions, with time along
%   the first or row dimension, and columns being different particles.  
%
%   XC and YC are cell arrays of closed contours, with length SIZE(XA,1).  
%   Thus XC and YC specify a different closed contour at each time. 
% 
%   ISENCLOSED is then a boolean array of size SIZE(XA) that is true for 
%   particles that are instantaneous enclosed within the closed contour 
%   corresponding to that time, and false otherwise.  
%   
%   To find the particles enclosed in the largest closed streamline, use 
%   [XCP,YPC]=LARGESTCLOSEDSTREAMLINE(...) followed by XC=CELLADD(XPC,XO)
%   and YC=CELLADD(YPC,YO) where XO and YO are the eddy center locations. 
%
%   Usage: isEnclosed=isParticleEnclosed(xa,ya,xc,yc)

arguments (Input)
    x (:,:) {mustBeNumeric,mustBeReal,mustBeFinite}
    y (:,:) {mustBeNumeric,mustBeReal,mustBeFinite,mustHaveSameSize(x,y)}
    xc (:,1) cell {mustHaveSameRowLength(xc,x)}
    yc (:,1) cell {mustHaveSameRowLength(yc,y)}
end

arguments (Output)
    isEnclosed {mustBeNumericOrLogical,mustBeReal,mustBeFinite}
end

isEnclosed = false(size(x));
[xo,yo]=curvemoments(xc,yc);

for k = 1:length(xc)
    %distance from the contour to its own center
    Rk = sqrt((xc{k}-xo(k)).^2+(yc{k}-yo(k)).^2);
    Rmin = min(Rk);
    Rmax = max(Rk);

    %distance from each point to the center of the contour 
    d = sqrt((x(k,:)-xo(k)).^2+(y(k,:)-yo(k)).^2);
    isEnclosed(k,d>Rmax) = false;%these cannot be enclosed
    isEnclosed(k,d<=Rmin) = true;%these must be enclosed

    %only need to use inpolygon for those in between max and min radius 
    index=find(d>Rmin&d<=Rmax);
    isEnclosed(k,index) = (inpolygon(x(k,index),y(k,index),xc{k},yc{k})==1);
end
end
