function [x,fval,exitflag,output,lambda,grad,hessian] = constrainedMinimize(x0,lb,ub,MaxIterations_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions(options,'Display', 'iter-detailed');
options = optimoptions(options,'MaxIterations', MaxIterations_Data);
options = optimoptions(options,'PlotFcn', {  @optimplotx @optimplotfval });
options = optimoptions(options,'Diagnostics', 'off');
options = optimoptions(options,'StepTolerance', 1e-12);
options = optimoptions(options,'ConstraintTolerance', 1e-12);
options = optimoptions(options,'FiniteDifferenceType', 'forward');
[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(@computeCostNew,x0,[],[],[],[],lb,ub,[],options);
