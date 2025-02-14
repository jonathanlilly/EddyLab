function eddy_model = analyticalEddyModel(eddyPath,params,eddyShapeString)
% analyticalEddyModel Creates an analytical model of an eddy
%
% Inputs:
%   eddyPath - Struct containing path functions:
%     .xe - Function handle @(t) for x position
%     .ye - Function handle @(t) for y position
%   params - Struct containing eddy parameters:
%     .A - Amplitude (meters)
%     .L - Length scale (meters)
%   eddy_function - (optional) Function handle of form 
%                   @(x,y,t,A,L,xe,ye) for custom eddy shape
%
% Output:
%   eddy_model - Function handle @(x,y,t) representing the eddy

% arguments
%     eddyPath struct
%     params struct
%     % eddyShapeString {mustBeString}
% end
%take variables out of the structure
use params
% set default eddy function as Gaussian
if isempty(eddyShapeString)
    eddyShape = @(x,y,t,A,L,xe,ye) A.*exp(-((x-xe(t)).^2 + (y-ye(t)).^2)/L^2);
end

if eddyShapeString == 'Gaussian'
    eddyShape = @(x,y,t,A,L,xe,ye) A.*exp(-((x-xe(t)).^2 + (y-ye(t)).^2)/L^2);

else eddyShapeString == 'Ellipse'
    if isempty(params.thetaDot)
        thetaDot=0;
    end
    eddyShape = @(x,y,t,A,La,Lb,thetaDot,xe,ye) A.*exp(-(...
        (((x-xe(t)).*cos(t*thetaDot) + (y-ye(t)).*sin(t*thetaDot))/La).^2 + ...
        ((-(x-xe(t)).*sin(t*thetaDot) + (y-ye(t)).*cos(t*thetaDot))/Lb).^2));
end

% make an eddy function with a chosen set of parameters (x,y,t)
eddy_model = @(x,y,t) eddyShape(x,y,t,A,L,eddyPath.xe,eddyPath.ye);
