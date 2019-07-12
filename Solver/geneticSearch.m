function [x,fval,exitflag,output,population,score] = geneticSearch(nvars,lb,ub,PopulationSize_Data)
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimoptions('ga');
%% Modify options setting
options = optimoptions(options,'PopulationSize', PopulationSize_Data);
options = optimoptions(options,'Display', 'iter');
options = optimoptions(options,'PlotFcn', {  @gaplotbestf @gaplotbestindiv @gaplotscores });
[x,fval,exitflag,output,population,score] = ...
ga(@computeCost,nvars,[],[],[],[],lb,ub,[],[],options); %computeCost fitnessFunc
