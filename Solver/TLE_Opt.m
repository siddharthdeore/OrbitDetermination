%%
% first run TLE_simulate.m to generate observables
% then run this file to see if solution converges
% for constraint minimization check lower and upper bounds since my
% selection of lower and upper bound might not be perfect

clc; clear; close all;
addpath("../common","vallado");
global tle1 tle2 opsmode whichconst terms typerun typeinput min_in_day timezone rho_obs stationLat stationLon stationAlt xp yp dut1 lod ddpsi ddeps dx dy dat jd_obs  RA_comp DEC_comp

%% load data
data=load('data.txt'); % data file format |jd lat lon az el rtasc dec rho drho|
jd_obs=data(1:20,1); % julian date
RA_obs=data(1:20,2);
DEC_obs=data(1:20,3);
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

tle1='1 25544U 98067A   19185.38548443  .00000071  00000-0  90797-5 0  9998';
tle2='2 25544  51.7428 275.8460 0007157 105.3654 340.4454 15.50956609177925';
% see difference between tle2 and tle_true
tle_true='2 25544  51.6428 275.9460 0007157 105.5654 340.6454 15.50956609177925';
%                   inc     RA      ecc       argp     mo      no

% fetch EOP data
[xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(datetime([2019,08,06]),'yyyy mm dd'));
[~, ~, ~, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
% Apriori estimate
Xest=[ satrec.ecco; satrec.inclo; satrec.nodeo; satrec.argpo; satrec.mo;str2double(tle2(52:63));satrec.bstar]';

% Apriori cost
fprintf('initial cost = %f',computeCost(Xest));
%% Optimization function
lb=[0.0000001,   0,    0,    0,    0,  0.01,       0]; % lower bound
ub=[0.0010000,  pi, 2*pi, 2*pi, 2*pi,   30,    0.01]; % upper bound
galb=[satrec.ecco-0.01, satrec.inclo-1*pi/180, satrec.nodeo-1*pi/180, satrec.argpo-1*pi/180, satrec.mo-1,str2double(tle2(52:63))-1,satrec.bstar-1e-4]; % lower bound
gaub=[satrec.ecco+0.01, satrec.inclo+1*pi/180, satrec.nodeo+1*pi/180, satrec.argpo+1*pi/180, satrec.mo+1,str2double(tle2(52:63))+1,satrec.bstar+1e-4]; % upper bound
MaxIterations=200;
typeAlgo=input('enter algo type \n 1. fminsearch \n 2. fmincon \n 3. GA \n');
switch typeAlgo
    case 1
        Xest= minUnconstrained(Xest,MaxIterations);
    case 2
        Xest=constrainedMinimize(Xest,lb,ub,MaxIterations);
    case 3
        PopulationSize=5;
        Xest=geneticSearch(7,galb,gaub,PopulationSize);
end
fprintf(1,'%s\n',tle2);
[tle1,tle2]=coe2tle(tle1,tle2,Xest(1),Xest(2),Xest(3),Xest(4),Xest(5),Xest(6));
fprintf(2,'%s\n',tle2);
fprintf(tle_true)
figure('Name','RA Residuals')
plot(RA_comp-RA_obs,'r*')
xlabel('observations')
ylabel('rtasc [rad]')
figure('Name','DEC Residuals')
plot(DEC_comp-DEC_obs,'r*')
xlabel('observations')
ylabel('decl [rad]')

fprintf('\necc  : %8.6f \n',Xest(1));
fprintf('incl : %8.6f  rad \n',Xest(2));
fprintf('node : %8.6f  rad \n',Xest(3));
fprintf('argp : %8.6f  rad \n',Xest(4));
fprintf('M0   : %8.6f \n',Xest(5));
fprintf('n0   : %8.6f  rev/day \n',Xest(6));
fprintf('bstar: %8.6f   \n',Xest(7));

%% some own functions
function [rho] = getRhoVector(rtasc,dec)
% retruns a vector from RA and DEC
cd=cos(dec);
rho=[cos(rtasc).*cd, cd.*sin(rtasc), sin(dec)];
end
