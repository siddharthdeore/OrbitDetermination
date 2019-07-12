clc; clear; close all;
addpath("../common","vallado");
global tle1 tle2 opsmode whichconst terms typerun typeinput min_in_day timezone rho_obs stationLat stationLon stationAlt xp yp dut1 lod ddpsi ddeps dx dy dat jd_obs

%% recorded on 2019 05 08
% USA 99 (MILSTAR-1 1)
%1 22988U 94009A   19128.29665602 -.00000262  00000-0  00000+0 0  9998
%2 22988  13.4618  55.8107 0002338   5.4549 232.5550  1.00270087  5601

%% load data
data=load('data.txt'); % data file format |jd lat lon az el rtasc dec rho drho|
jd_obs=data(:,1); % julian date
RA_obs=data(:,2);
DEC_obs=data(:,3);
rho_obs=getRhoVector(RA_obs,DEC_obs);
%% Ground Station
stationLat=41.89*pi/180; % rad
stationLon=12.49*pi/180; % rad
stationAlt=0; % m
%% SGP4 configration
% input arguments
opsmode= 'a'; whichconst = 84; terms = 2; typerun='c'; typeinput = 'e';
min_in_day=1440;
timezone=0;

tle1='1 22988U 94009A   19128.29665602 -.00000262  00000-0  00000+0 0  9998';
tle2='2 22988  13.4618  56.8107 0003338   5.5549  232.555  1.03270087  5601';
%                inc     RA      ecc       argp     mo      no
%tle2='2 22988  13.4618  55.8107 0002338   5.4549 232.5550  1.00270087 5601'


% fetch EOP data
[xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(datetime([2019,05,08]),'yyyy mm dd'));
[~, ~, ~, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
% Apriori estimate
Xest=[ satrec.ecco; satrec.inclo; satrec.nodeo; satrec.argpo; satrec.mo;str2double(tle2(52:63));]';
% Apriori cost
fprintf('initial cost = %f',computeCost(Xest));
%% Optimization function
%Xest = gamultiobj(@computeCost,6);
lb=[0.0001 0.0001 0.0001 0.0001 1, 0.5];
ub=[0.5 1.14 1.14 1.14 1 2.5];
options = optimoptions('gamultiobj','PlotFcn',@gaplotparetodistance);
[Xest,fval,exitflag,output,population,scores]  = gamultiobj(@computeCost,6,[],[],[],[],lb,ub,options);
fprintf(2,'%s\n',tle2);
[tle1,tle2]=coe2tle(tle1,tle2,Xest(2),Xest(3),Xest(1),Xest(4),Xest(5),Xest(6));
fprintf('2 22988  13.4618  55.8107 0002338   5.4549 232.5550  1.00270087  5601')

%% some own functions
function [rho] = getRhoVector(rtasc,dec)
% retruns a vector from RA and DEC
cd=cos(dec);
rho=[cos(rtasc).*cd, cd.*sin(rtasc), sin(dec)];
end