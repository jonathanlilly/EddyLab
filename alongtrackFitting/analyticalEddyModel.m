function eddy_model = analyticalEddyModel(eddyPath,params,eddy_function)
% Input:
% xe(t), ye(t) are function handles that represents the eddy paths as a function of t
% params are parameters for the analytical eddy
arguments (Input)
    eddyPath.xe (1,1) function_handle %function of t
    eddyPath.ye (1,1) function_handle %function of t
    params struct %A L
    eddy_function function_handle
end
% set default eddy function as Gaussian
if isempty(eddy_function)
    eddy_function = @(x,y,t,A,L,xe,ye) A.*exp(-((x-xe(t)).^2 + (y-ye(t)).^2)/L^2);
end

% make an eddy function with a chosen set of parameters (x,y,t)
eddy_model = @(x,y,t) test_eddy(x,y,t,0.15,80e3,x_lin,y_lin);
