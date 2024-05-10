function[xo,yo,xc,yc] = eddyCenter(x,y,z,level,options)
%eddyCenter  Find the center of a coherent eddy from 3D (x,y,t) fields.
%
%   [XO,YO] = eddyCenter(X,Y,Z,LEVEL) finds the center of an eddy as the
%   centroid of closed contours where the field Z takes on the value LEVEL. 
%
%   X and Y are uniformly spaced 1D arrays, Z is a 3D array with LENGTH(Y)
%   rows and LENGTH(X) columns, and XO and YO are 1D arrays with SIZE(Z,3) 
%   rows. The third dimension of Z is interpreted as time. 
%
%   In the case that more than one matching closed contour of Z(:,:,K) is 
%   found for any K, the largest area closed contour at that K is chosen. 
%   Note that this assumes equal grid spacings X(2)-X(1) and Y(2)-Y(1).
%
%   Options
%
%   [XO,YO] = eddyCenter(X,Y,Z,LEVEL,method = "timesMax") instead used the 
%   closed contours at which Z equals LEVEL times its maximum at each time. 
%
%   [XO,YO] = eddyCenter(X,Y,Z,LEVEL,method = "timesMin") instead used the 
%   closed contours at which Z equals LEVEL times its minimum at each time.   
%
%   [XO,YO] = eddyCenter(..., contains = {X1,Y1}) additionally specifies a
%   point that the chosen contour must contain.  X1 and Y1 are arrays of
%   length SIZE(Z,3).  This is applied before the largest area criterion.
%
%   [XO,YO,XC,YC] = eddyCenter(...) optionally returns cell arrays of the 
%   same length as XO and YO containing the closed contours at each time. 
%   These can be plotted with cellplot(XC,YC).
%
%   Examples
%
%   [XO,YO] = eddyCenter(X,Y,ZETA,0) where ZETA is a 3D array of relative 
%   vorticity finds the eddy center as the centroid of the zero vorticity 
%   contour at each time.
%
%   [XO,YO] = eddyCenter(X,Y,ZETA,0.9,method = "timesMin") finds the eddy 
%   center as the centroid of the curve at which the vorticity is equal to
%   0.9 times its minimum at each time.
% 
%   [XO,YO] = eddyCenter(X,Y,ZETA,1.0,method = "timesMin") finds the eddy 
%   center as the curve tracing out the minimum vorticity at each time. 
%
%   Usage: [xo,yo] = eddyCenter(x,y,z,level);
%          [xo,yo,xc,yc] = eddyCenter(x,y,z,level);
%          [xo,yo,xc,yc] = eddyCenter(x,y,z,level,method = "timesMax");

arguments (Input)
    x {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(x)}
    y {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(y),mustHaveSameSpacing(x,y)}
    z {mustBeNumeric,mustBeReal,mustBeCompatible(x,y,z)}
    level (1,1) double {mustBeReal,mustBeFinite}
    options.Method (1,:) string ...
        {mustBeMember(options.Method,["timesMax","timesMin","value"])} = "value"
    options.Contains (1,2) cell = {[],[]}
end

arguments (Output)
    xo {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    yo {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    xc (:,1) cell 
    yc (:,1) cell 
end

if strcmpi(options.Method,"timesMax") || strcmpi(options.Method,"timesMin")
    [extrema,xo,yo,xc,yc] = findExtrema(x,y,z,options.Method);
    if level ~= 1 
        level = extrema * level;
        [xo,yo,xc,yc] = centerFromContours(x,y,z,level,options.Contains);
    end
else
    [xo,yo,xc,yc] = centerFromContours(x,y,z,level,options.Contains);
end

end
%--------------------------------------------------------------------------
function[extrema,xo,yo,xc,yc] = findExtrema(x,y,z,method)

[extrema,row,col] = deal(zeros(size(z,3),1));

for k = 1:size(z,3)
    if strcmpi(method,"timesMax")
        [extrema(k),index] = max(z(:,:,k),[],"all");
    elseif strcmpi(method,"timesMin")
        [extrema(k),index] = min(z(:,:,k),[],"all");
    end
    [row(k),col(k)] = ind2sub(size(z(:,:,1)),index);
end

xo=interp1(1:length(x),x,col);
yo=interp1(1:length(y),y,row);

[xc,yc] = deal(cell(size(z,3),1));
%set cells to single point for case of level = 1
for k = 1:size(z,3)
    xc{k}=xo;
    yc{k}=yo;
end

end
%--------------------------------------------------------------------------
function[xo,yo,xc,yc] = centerFromContours(x,y,z,level,contains)

if length(level)==1
    %In this case, replicate LEVEL across all times
    level = level + zeros(size(z,3),1);
end

[xo,yo] = deal(zeros(size(z,3),1));
[xc,yc] = deal(cell(size(z,3),1));

for k = 1:size(z,3)
    [col,row] = closedcurves(z(:,:,k),level(k));

    %keep only those contours containing a specified point
    if ~isempty(contains{1})
        %interpolate (x,y) for contained point onto (i,j) indices
        row1=interp1(y,1:length(y),contains{2}(k));
        col1=interp1(x,1:length(x),contains{1}(k));
        
        bool=false(length(col),1);
        for l = 1:length(col)
            bool(l)=(inpolygon(col1,row1,col{l},row{l})==1);%convert to boolean
        end
        col=col(bool);
        row=row(bool);
    end

    %for multiple closed curves, choose the one with the largest area
    if length(col) == 1
        col = col{1};
        row = row{1};
    else
        [~,~,~,R] = curvemoments(col,row);
        [~,ii] = max(R,[],"all");
        col = col{ii};
        row = row{ii};
    end

    xc{k}=interp1(1:length(x),x,col);
    yc{k}=interp1(1:length(y),y,row);
    [xo(k),yo(k)] = curvemoments(xc{k},yc{k});
end

end
