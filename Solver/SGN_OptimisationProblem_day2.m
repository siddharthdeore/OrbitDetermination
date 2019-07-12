%% solver for observation taken on 22-03-2019
clc; clear; close all;
addpath("../common","vallado");
global tle1 tle2 opsmode whichconst terms typerun typeinput min_in_day timezone rho_obs stationLat stationLon stationAlt xp yp dut1 lod ddpsi ddeps dx dy dat jd_obs DEC_comp RA_comp

%% load data
data=load('observationsAll.txt'); % data file format |jd lat lon az el rtasc dec rho drho|
jd_obs=data(13:25,1); % julian date
RA_obs=data(13:25,2);
DEC_obs=data(13:25,3);
%%
% 1:14 -21
% 14:26 -22
% 27:end -23
rho_obs=getRhoVector(RA_obs,DEC_obs);
%% Ground Station
stationLat=41.95778*pi/180; % rad
stationLon=12.50556*pi/180; % rad
stationAlt=76; % m
%% SGP4 configration
% input arguments
opsmode= 'a'; whichconst = 84; terms = 2; typerun='c'; typeinput = 'e';
min_in_day=1440;
timezone=0;

tle1='1 27422U 02021B   19080.74912824 -.00000022  00000-0  75424-5 0  9999';
tle2='2 27422  98.4176  57.1153 0008279 149.2393 167.4261 14.09156672878468';
tle_true='2 27422  98.2911  57.2753 0013167 129.0062 231.2296 14.29268859878468';
%                   inc     RA      ecc       argp     mo      no
% fetch EOP data
[xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(datetime([2019,03,22]),'yyyy mm dd'));
[~, ~, ~, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
% Apriori estimate
Xest=[ satrec.ecco; satrec.inclo; satrec.nodeo; satrec.argpo; satrec.mo;str2double(tle2(52:63));satrec.bstar]';

% Apriori cost
fprintf('initial cost = %f',computeCostNew(Xest));
%% Optimization function
tic
lb=[0.0000001, 0,    0,    0,    0,  0.0001,       0]; % lower bound
ub=[0.9999999,  2*pi-0.001, 2*pi, 2*pi, 2*pi,   30,    0.01]; % upper bound
MaxIterations=200;
typeAlgo=input('enter algo type \n 1. fminsearch \n 2. fmincon \n 3. GA \n');
switch typeAlgo
    case 1
        Xest= minUnconstrained(Xest,MaxIterations);
    case 2
        Xest=constrainedMinimize(Xest,lb,ub,MaxIterations);
    case 3
        PopulationSize=20;
        Xest=geneticSearch(7,lb,ub,PopulationSize);
end
toc
fprintf(1,'%s\n',tle2);
[tle1,tle2]=coe2tle(tle1,tle2,Xest(1),Xest(2),Xest(3),Xest(4),Xest(5),Xest(6));
fprintf(2,'%s\n',tle2);
fprintf(tle_true)

fprintf('\necc  : %8.6f \n',Xest(1));
fprintf('incl : %8.6f  rad \n',Xest(2));
fprintf('node : %8.6f  rad \n',Xest(3));
fprintf('argp : %8.6f  rad \n',Xest(4));
fprintf('M0   : %8.6f \n',Xest(5));
fprintf('n0   : %8.6f  rev/day \n',Xest(6));
fprintf('bstar: %8.6f   \n',Xest(7));

figure('Name','RA Residuals')
plot(RA_comp-RA_obs,'r*')
xlabel('observations')
ylabel('rtasc [rad]')
figure('Name','DEC Residuals')
plot(DEC_comp-DEC_obs,'r*')
xlabel('observations')
ylabel('decl [rad]')
%% some own functions
function [rho] = getRhoVector(rtasc,dec)
% retruns a vector from RA and DEC
cd=cos(dec);
rho=[cos(rtasc).*cd, cd.*sin(rtasc), sin(dec)];
end
