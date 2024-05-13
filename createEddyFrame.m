function[x1,y1,fieldOut]=createEddyFrame(x,y,xo,yo,L,fieldIn,options)
%createEddyFrame  Create a frame centered on an eddy through interpolation.
%
%   [XP,YP,ZP] = createEddyFrame(X,Y,XO,YO,L,Z) interpolates the field Z 
%   defined at spatial locations (X,Y) into a frame moving with the eddy
%   center point (XO,YO), in a square domain with side lengths set by L.
% 
%   The "P" in XP, YP, and ZP could be read as "prime", distinguishing 
%   variables centered on the eddy from their original versions.
%
%   [XP,YP,ZP1,ZP2,...,ZPN] = createEddyFrame(X,Y,XO,YO,L,Z1,Z2,...,ZN) 
%   also works.  That is, Z and ZP are repeating arguments. 
% 
%   X and Y are uniformly spaced 1D arrays, Z is a 3D array with LENGTH(Y)
%   rows and LENGTH(X) columns, and XO and YO are 1D arrays with SIZE(Z,3) 
%   rows, see EDDYCENTER. The third dimension of Z is interpreted as time. 
%
%   XP and YP are uniformly spaced 1D arrays centered on zero, and with
%   XP(END)-XP(1) and YP(END)-YP(1) being as close as possible to L without
%   exceeding it given their sampling interval, set as discussed shortly.
%   The units of L are interpreted as being the same as those of X and Y.
%
%   Options
%
%   [XP,YP,ZP] = createEddyFrame(...,dx = DX) sets the sampling interval of
%   the output coordinates, XP(2)-XP(1) and YP(2)-YP(1). By default, the
%   the sampling interval of the input coordinates (X,Y) will be used.  
%
%   [XP,YP,ZP] = createEddyFrame(...,method = METHOD) sets the 
%   interpolation method, with valid choices specified in INTERP2. The
%   default choice is method = "linear". 
% 
%   Usage: [xp,yp,zp] = createEddyFrame(x,y,xo,yo,L,z);
%          [xp,yp,zp] = createEddyFrame(x,y,xo,yo,L,z,dx = 2);
%          [xp,yp,zp] = createEddyFrame(x,y,xo,yo,L,z,method = "spline");


arguments (Input)
    x {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(x)}
    y {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(y),mustHaveSameSpacing(x,y)}
    xo {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    yo {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector} 
    L (1,1) double {mustBeReal,mustBeFinite}
end
arguments (Input, Repeating)
    fieldIn {mustBeNumeric,mustBeReal,mustBeCompatible(x,y,fieldIn)} 
end
arguments (Input)
    options.dx(1,1) double {mustBeReal,mustBeFinite} 
    options.Method (1,:) string ...
        {mustBeMember(options.Method,["nearest","linear","spline","cubic","makima"])} = "linear"
end

arguments (Output)
    x1 {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(x1)}
    y1 {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(y1)}
end
arguments (Output,Repeating)
    fieldOut {mustBeNumeric,mustBeReal} 
end

if isfield(options,'dx')
    dx = options.dx;
else 
    dx = (x(2)-x(1));
end

%create a new, centered grid 
x1 = (0:dx:L/2)';
x1 = [-flipud(x1(2:end));x1];
y1 = x1;

[xg,yg]=meshgrid(x1,y1);

for l=1:length(fieldIn)
    for k = 1:size(fieldIn{l},3)
        fieldOut{l}(:,:,k) = interp2(x,y,fieldIn{l}(:,:,k),xg+xo(k),yg+yo(k),options.Method);
    end
end

end


