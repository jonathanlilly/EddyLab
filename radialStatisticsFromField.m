function[mz,rmid,numz,stdz1,stdz2] = radialStatisticsFromField(x,y,z,options)
%radialStatisticsFromField  Radial statistics from 3D (x,y,t) fields.
%
%   [MZ,RMID,TMID,NUMZ,STDZ1,STDZ2] = radialStatisticsFromField(X,Y,Z)
%   computes radial profiles of azimuthal and temporal statistics from 
%   gridded data.
% 
%   Here X and Y are uniformly spaced 1D arrays, and Z is a 3D array 
%   with LENGTH(Y) rows and LENGTH(X) columns. The third dimension of Z is 
%   interpreted as time. 
% 
%   Radial bins are chosen automatically to have edges running from zero to
%   the maximum radial distance from the origin implied by X and Y. The 
%   bin spacing is chosen to be the same as that of X and Y, X(2)-X(1).
%
%   All output fields are 1D arrays of the same size.  MZ is the mean in
%   radial bins with midpoints RMID and STDZ1 and STDZ2 are two different
%   standard deviations, computed as described subsequently.  NUMZ is the 
%   number of data points found in each bin.
%
%   Averaging order
%
%   This function returns output fields with different interpretations
%   depending on whether azimuthal or temporal averages are computed first.
%
%   [...] = radialStatisticsFromField(...,firstAverage = "azimuthal") 
%   performs an azimuthal average first.  The output statistics are 
%
%      MZ    — Azimuthal average of Z followed by temporal average
%      STDZ1 — Square root of temporal variance of azimuthal average
%      STDZ2 — Square root of temporal average of azimuthal variance
%
%   [...] = radialStatisticsFromField(...,firstAverage = "temporal") 
%   performs a temporal average first.  The output statistics are 
%
%      MZ    — Temporal average of Z followed by azimuthal average
%      STDZ1 — Square root of azimuthal variance of temporal average
%      STDZ2 — Square root of azimuthal average of temporal variance
%
%   The firstAverage argument is optional, and defaults to "azimuthal". 
%
%   In both cases the total variance about the time-mean, azimuthal-mean 
%   profile of Z is given by STDZ1.^2 + STDZ2.^2.
%
%   See also radialStatisticsFromScatter.
%
%   Usage:  [mz,rmid,numz,stdz1,stdz2] = radialStatisticsFromField(x,y,z)

arguments (Input)
    x {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(x)}
    y {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector,mustBeUniform(y),mustHaveSameSpacing(x,y)}
    z {mustBeNumeric,mustBeReal,mustBeCompatible(x,y,z)} 
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

[xg,yg] = meshgrid(x,y);
rg = sqrt(xg.^2+yg.^2);

dx = x(2)-x(1);
rbin = [-dx/2:dx:max(rg(:))]';

if strcmpi(options.firstAverage,"temporal")
    %mean and standard deviation in xy bins
    mzxy = mean(z,3,"omitnan");
    szxy = std(z,1,3,"omitnan");
    %these are functions of xy but not time

    %azimuthal average
    [mz,~,rmid,numz,stdz1] = twodstats(rg,rg,mzxy,[-inf inf],rbin); %don't use the first index
    stdz2 = sqrt(twodstats(rg,rg,szxy.^2,[-inf inf],rbin)); %azimuthal average of temporal variance
elseif strcmpi(options.firstAverage,"azimuthal")
    %statistics in r and t
    tbin = [0.5:1:size(z,3)+0.5]';
    rg=vrep(rg,size(z,3),3);
    tg = vrep(permute(1:size(z,3),[3 1 2]),[size(z,1) size(z,2)],[1 2]);
    [mzrt,rmid,~,numz,szrt] = twodstats(rg,tg,z,rbin,tbin);
    mz = mean(mzrt,1,"omitnan")';%average over azimuth then time
    numz =sum(numz,1,"omitnan")';%total number of points in each bin
    stdz1 = std(mzrt,1,1,"omitnan")';%temporal standard deviation of azimuthal average
    stdz2 = sqrt(mean(szrt.^2,1,"omitnan"))';%sqrt of temporal average of azimuthal variance
end
 
end