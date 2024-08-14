function[mz,rmid,numz,stdz1,stdz2] = radialStatisticsFromScatter(x,y,t,z,rbin,tbin,options)
%radialStatisticsFromScatter  Radial statistics from scattered data.
%
%   [MZ,RMID,NUMZ,STDZ1,STDZ2] = ...
%        radialStatisticsFromScatter(X,Y,T,Z,RBIN,TBIN) computes radial 
%   profiles of azimuthal and temporal statistics from scattered data.
% 
%   Here X, Y, T, and Z are arrays all having the same size, while RBIN and
%   TBIN are 1D arrays giving bin edges for radius R=SQRT(X.^2+Y.^2) and 
%   time T, respectively.  
% 
%   All output fields are 1D arrays of the same size.  MZ is the mean in
%   radial bins with midpoints RMID and STDZ1 and STDZ2 are two different
%   standard deviations, computed as described subsequently.  NUMZ is the 
%   number of data points found in each bin.
%
%   The output fields have two different interpretations depending on
%   whether the function is called with the firstAverage = "azimuthal" or
%   firstAverage = "temporal" option.  See the documentation at 
%   radialStatisticsFromField for further details. 
%
%   The functionality of radialStatisticsFromScatter is the same as that
%   of radialStatisticsFromField, but the former is intended to be used
%   with irregularly sampled data, such as alongtrack altimeter data, and
%   the latter with gridded model or observational data.
%
%   See also radialStatisticsFromField.
%
%   Usage:  [mz,rmid,numz,stdz1,stdz2] = ...
%                 radialStatisticsFromScatter(x,y,t,z,rbin,tbin)


arguments (Input)
    x {mustBeNumeric,mustBeReal}
    y {mustBeNumeric,mustBeReal,mustHaveSameSize(x,y)}
    t  {mustBeNumeric,mustBeReal,mustHaveSameSize(x,t)}
    z  {mustBeNumeric,mustBeReal,mustHaveSameSize(x,z)}
    rbin  {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    tbin  {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(tbin)}
    options.firstAverage (1,:) string ...
        {mustBeMember(options.firstAverage,["temporal","azimuthal"])} = "azimuthal"
end

arguments (Output)
    mz    {mustBeNumeric,mustBeReal,mustBeVector}
    rmid  {mustBeNumeric,mustBeReal,mustBeVector}
    numz  {mustBeNumeric,mustBeReal,mustBeVector}
    stdz1 {mustBeNumeric,mustBeReal,mustBeVector}
    stdz2 {mustBeNumeric,mustBeReal,mustBeVector}
end

if strcmpi(options.firstAverage,"temporal")
    %mean and standard deviation in xy bins
    dx = rbin(2)-rbin(1);
    xbin = -rbin(end):dx:rbin(end);
    [mzxy,xmid,ymid,~,szxy] = twodstats(x,y,z,xbin,xbin);
    [xg,yg] = meshgrid(xmid,ymid);
    rg = sqrt(xg.^2+yg.^2);
    %these are functions of xy but not time

    %azimuthal average
    [mz,~,rmid,numz,stdz1] = twodstats(rg,rg,mzxy,[-inf inf],rbin); %don't use the first index
    stdz2 = sqrt(twodstats(rg,rg,szxy.^2,[-inf inf],rbin)); %azimuthal average of temporal variance
elseif strcmpi(options.firstAverage,"azimuthal")
    %statistics in r and t
    [mzrt,rmid,~,numz,szrt] = twodstats(sqrt(x.^2+y.^2),t,z,rbin,tbin);
    mz = mean(mzrt,1,"omitnan")';%average over azimuth then time
    numz =sum(numz,1,"omitnan")';%total number of points in each bin
    stdz1 = std(mzrt,1,1,"omitnan")';%temporal standard deviation of azimuthal average
    stdz2 = sqrt(mean(szrt.^2,1,"omitnan"))';%sqrt of temporal average of azimuthal variance
end

end