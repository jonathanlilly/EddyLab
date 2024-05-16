function[xo,yo,xc,yc] = eddyCenterFromField(x,y,z,level,options)
%eddyCenterFromField  Find the center of a coherent eddy from 3D fields.
%
%   [XO,YO] = eddyCenterFromField(X,Y,Z,LEVEL) finds the center of an eddy 
%   as the centroid of closed contours where the field Z takes on the value 
%   LEVEL. The "O" here can be thought of as standing for "origin", since 
%   it will become the origin of an eddy-centered coordinate system.
%
%   X and Y are uniformly spaced 1D arrays, Z is a 3D array with LENGTH(X)
%   rows and LENGTH(Y) columns, and XO and YO are 1D arrays with SIZE(Z,3) 
%   rows. The third dimension of Z is interpreted as time. 
%
%   In the case that more than one matching closed contour of Z(:,:,K) is 
%   found for any K, the largest area closed contour at that K is chosen. 
%   Note that this assumes equal grid spacings X(2)-X(1) and Y(2)-Y(1).
%
%   Options
%
%   [XO,YO] = eddyCenterFromField(X,Y,Z,LEVEL,method = "timesMax") uses the 
%   closed contours at which Z equals LEVEL times its maximum at each time. 
%
%   [XO,YO] = eddyCenterFromField(X,Y,Z,LEVEL,method = "timesMin") uses the
%   closed contours at which Z equals LEVEL times its minimum at each time.   
%
%   [XO,YO] = eddyCenterFromField(..., contains = {X1,Y1}) specifies a
%   point that the chosen contour must contain.  X1 and Y1 are arrays of
%   length SIZE(Z,3).  This is applied before the largest area criterion.
%
%   [XO,YO,XC,YC] = eddyCenterFromField(...) optionally returns cell arrays 
%   of the same length as XO and YO containing the closed contours at each
%   time. These can be plotted with cellplot(XC,YC).
%
%   Examples
%
%   [XO,YO] = eddyCenterFromField(X,Y,ZETA,0) where ZETA is a 3D array of 
%   relative vorticity finds the eddy center as the centroid of the zero 
%   vorticity contour at each time.
%
%   [XO,YO] = eddyCenterFromField(X,Y,ZETA,0.9,method = "timesMin") finds
%   the eddy center as the centroid of the curve at which the vorticity is 
%   equal to 0.9 times its minimum at each time.
% 
%   [XO,YO] = eddyCenterFromField(X,Y,ZETA,1.0,method = "timesMin") finds 
%   the eddy center as the curve tracing out the vorticity minimum at each
%   time. 
%
%   Usage: [xo,yo] = eddyCenterFromField(x,y,z,level);
%          [xo,yo,xc,yc] = eddyCenterFromField(x,y,z,level);
%          [xo,yo,xc,yc] = ...
%             eddyCenterFromField(x,y,z,level,method = "timesMax");

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

xo=interp1(1:length(x),x,row);
yo=interp1(1:length(y),y,col);

[xc,yc] = deal(cell(size(z,3),1));
%set cells to single point for case of level = 1
for k = 1:size(z,3)
    xc{k}=xo;
    yc{k}=yo;
end

end
%--------------------------------------------------------------------------
function[xo,yo,xc,yc] = centerFromContours(x,y,z,level,contains)

if isscalar(level)
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
        row1=interp1(x,1:length(x),contains{1}(k));
        col1=interp1(y,1:length(y),contains{2}(k));
        
        bool=false(length(col),1);
        for l = 1:length(col)
            bool(l)=(inpolygon(row1,col1,row{l},col{l})==1);%convert to boolean
        end
        col=col(bool);
        row=row(bool);
    end

    %for multiple closed curves, choose the one with the largest area
    if isscalar(col)
        col = col{1};
        row = row{1};
    else
        [~,~,~,R] = curvemoments(row,col);
        [~,ii] = max(R,[],"all");
        col = col{ii};
        row = row{ii};
    end

    xc{k}=interp1(1:length(x),x,row);
    yc{k}=interp1(1:length(y),y,col);
    [xo(k),yo(k)] = curvemoments(xc{k},yc{k});
end

end
