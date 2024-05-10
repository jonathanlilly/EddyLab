function[xc,yc,levels,firstClosedLevel] = largestClosedStreamline(xp,yp,psip,ci)
%largestClosedStreamline  Largest closed streamline in the eddy frame.
%
%   [XPC,YPC] = largestClosedStreamline(XP,YP,PSIP,CI) finds the largest 
%   closed streamline of PSIP, the streamfunction of the flow associated 
%   with an eddy relative to a frame of reference moving with the eddy. 
% 
%   XP and YP are coordinate arrays in the co-moving frame, with the origin
%   (0,0) corresponding to the eddy center, created by CREATEEDDYFRAME.
%
%   PSIP is a 3D array of streamfunction values for the velocity in the 
%   co-moving frame, created by EDDYFRAMESTREAMFUNCTION. 
%
%   CI is the desired contour interval, a scalar. 
%
%   XPC and YPC are cell arrays of length SIZE(PSIP,3) containing the 
%   largest closed streamline at each time, expressed in the eddy-centered
%   coordinate system.  These can be plotted with CELLPLOT(XPC,YPC).
%
%   [XPC,YCP,LEVELS] = largestClosedStreamline(...) also returns the 
%   streamline contour levels that are searched, with spacing CI. 
%
%   [XPC,YPC,LEVELS,LEVELINDEX] = largestClosedStreamline(...) also returns 
%   LEVELINDEX, a length SIZE(PSIP,3) index into LEVELS.  The streamline
%   values for each contour are given by LEVELS(LEVELINDEX);
%
%   In the case that no closed streamline is found at a particular time,
%   the corresponding cells of XPC and YPC will be empty, and the
%   corresponding value of INDEX will be NaN.
%
%   Usage: [xpc,ypc,levels,index] = largestClosedStreamline(xp,yp,psip,0.5)

%

arguments (Input)
    xp {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(xp)}
    yp {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(yp),mustHaveSameSpacing(xp,yp)}
    psip {mustBeNumeric,mustBeReal,mustBeCompatible(xp,yp,psip)} 
    ci (1,1) double {mustBeReal,mustBeFinite}
end

arguments (Output)
    xc (:,1) cell 
    yc (:,1) cell 
    levels {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    firstClosedLevel {mustBeNumeric,mustBeReal,mustBeVector}
end

[xg,yg] = meshgrid(xp,yp);
[xc,yc] = deal(cell(size(psip,3),1));

for k = 1:size(psip,3)
    xc{k} = [];
    yc{k} = [];
end

%ensure levels at standard intervals
a = floor(min(psip,[],"all")./ci)*ci;
b = ceil(max(psip,[],"all")./ci)*ci;
levels = a:ci:b;

%we want to proceed from outward in, so if the central value of the
%streamfunction is near the *mininum* streamfunction, flip the levels
psipo = mean(psip((end+1)/2,(end+1)/2,:),3);
if abs(levels(1)-psipo) < abs(levels(end)-psipo)
    levels = fliplr(levels);
end

firstClosedLevel = nan*ones(size(psip,3),1);

for k = 1:size(psip,3)
    %disp(['Computing time step ' int2str(k) ' of ' int2str(size(psip,3)) '.'])
    %just need to find the first or outermost level with a closed contour
    %containing the origin, since they are nested

    if k == 1
        l = 2;
    else 
        l = firstClosedLevel(k-1);
    end

    done = false;

    xck = [];
    xckf = [];
    
    while ~done && (l>1 && l<= length(levels))
        
        %previous, or more distant, contour level
        %(last one that was not closed from the previous time step)
        if isempty(xckf)
            %slight hack, adding epsilon, because I need to update closedcurves... 
            %an error is arising when levels = 1.0000 
            [xckf,yckf] = closedcurves(xg,yg,psip(:,:,k),levels(l-1)+1e-10);  
            [xckf,yckf] = subsetEnclosingOrigin(xckf,yckf);
        end
      
        %current contour level
        if isempty(xck)
            [xck,yck] = closedcurves(xg,yg,psip(:,:,k),levels(l)+1e-10);
            [xck,yck] = subsetEnclosingOrigin(xck,yck);
        end
       % [k,length(xck),length(xckf)]

        if isempty(xckf) & ~isempty(xck) 
             %last one open and current one closed; we're done
             xc{k} = xck{1};
             yc{k} = yck{1};
             firstClosedLevel(k) = l;
             done = true;
        elseif isempty(xckf) & isempty(xck)
             %neither closed; move in
             l = l+1;
             xckf = xck;
             yckf = yck;
             xck = [];%compute new current level
        elseif ~isempty(xckf) & ~isempty(xck)
             %neither open; move out
             l = l-1;
             xck = xckf;
             yck = yckf;
             xckf = [];%compute new previous level
        end
        if ~(l>1 && l<= length(levels))
            xc{k} = [];
            yc{k} = [];
            firstClosedLevel(k) = nan;
        end
    end
end

end
%--------------------------------------------------------------------------
function [xc,yc] = subsetEnclosingOrigin(xc,yc)
%return only those contours that enclose the origin 
if ~isempty(xc)
    bool = false(length(xc),1);
    for m = 1:length(xc)
        bool(m) = (inpolygon(0,0,xc{m},yc{m})==1); %This converts to a boolean
    end
    xc = xc(bool);
    yc = yc(bool);
end

end
