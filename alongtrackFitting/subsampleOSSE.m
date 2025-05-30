function alongtrack = subsampleOSSE(alongtrack,eddy_field)
arguments
    alongtrack struct
    eddy_field struct
end
elapsed_time = alongtrack.t - min(alongtrack.t);
% subsampling alongtrack from QG model(i.e. eddy_field)
[xMat, yMat, tMat] = ndgrid(eddy_field.x, eddy_field.y, eddy_field.t-1);
alongtrack.ssh = interpn(xMat, yMat, tMat, eddy_field.ssh, alongtrack.x, alongtrack.y, elapsed_time, 'linear', 0);