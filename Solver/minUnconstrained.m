function [x,fval,exitflag,output] = minUnconstrained(x0,MaxIter_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimset;
%% Modify options setting
options = optimset(options,'Display', 'iter');
options = optimset(options,'MaxIter', MaxIter_Data);
options = optimset(options,'TolFun', 1e-10);
options = optimset(options,'TolX', 1e-10);
options = optimset(options,'PlotFcns', {  @optimplotx @optimplotfval });
[x,fval,exitflag,output] = fminsearch(@computeCostNew,x0,options);
end