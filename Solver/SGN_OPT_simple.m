clc; clear; close all;
addpath("../common","vallado");
global tle1 tle2 opsmode whichconst terms typerun typeinput min_in_day timezone rho_obs stationLat stationLon stationAlt  jd_obs
global  xp  yp  dut1  lod  ddpsi  ddeps  dx  dy dat

%% load data
data=load('observationsAll.txt'); % data file format |jd lat lon az el rtasc dec rho drho|
%jd_obs=data(1:12,1); % julian date
%RA_obs=data(1:12,2);
%DEC_obs=data(1:12,3);
%%
% 1:12 -21
% 13:26 -22
% 27:end -23
%% Ground Station
stationLat=41.95778*pi/180; % rad
stationLon=12.50556*pi/180; % rad
stationAlt=76; % m
%% SGP4 configration
% input arguments
opsmode= 'a'; whichconst = 84; terms = 2; typerun='c'; typeinput = 'e';
min_in_day=1440;
timezone=0;
dayCoice = input('\n Enter day numer 1,2,3 \n');
switch dayCoice
    case 1
        jd_obs=data(1:12,1); % julian date
        RA_obs=data(1:12,2);
        DEC_obs=data(1:12,3);

        tle1='1 27422U 02021B   19080.74912824 -.00000022  00000-0  75424-5 0  9999';
        tle2='2 27422  98.2911  57.2753 0013167 129.0062 231.2296 14.29268859878468';
        tle_true='2 27422  98.2911  57.2753 0013167 129.0062 231.2296 14.29268859878468';
        [xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(datetime([2019,03,21]),'yyyy mm dd'));
    case 2
        jd_obs=data(13:26,1); % julian date
        RA_obs=data(13:26,2);
        DEC_obs=data(13:26,3);
        tle1='1 27422U 02021B   19079.90905414 -.00000024 +00000-0 +67363-5 0  9997';
        tle2='2 27422 098.2909 056.4759 0013128 131.2716 228.9599 14.29268806880239';
        tle_true='2 27422 098.2909 056.4759 0013128 131.2716 228.9599 14.29268806880239';
         [xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(datetime([2019,03,22]),'yyyy mm dd'));
    case 3
        jd_obs=data(27:end,1); % julian date
        RA_obs=data(27:end,2);
        DEC_obs=data(27:end,3);
        tle1='1 27422U 02021B   19081.86922708 -.00000029  00000-0  45726-5 0  9994';
        tle2='2 27422  98.2913  58.3413 0013185 125.9815 234.2592 14.29268846878622';
        tle_true='2 27422  98.2913  58.3413 0013185 125.9815 234.2592 14.29268846878622';

        [xp,yp,dut1,lod,ddpsi,ddeps,dx,dy,dat] = fetchSpaceData(datestr(datetime([2019,03,23]),'yyyy mm dd'));
end
rho_obs=getRhoVector(RA_obs,DEC_obs);

%1 27422U 02021B   19080.74912824 -.00000022  00000-0  75424-5 0  9999
%2 27422  98.2911  57.2753 0013167 129.0062 231.2296 14.29268859878468
%                   inc     RA      ecc       argp     mo      no
[~, ~, ~, satrec] = twoline2rv(tle1, tle2, typerun, typeinput, opsmode, whichconst);
% Apriori estimate
Xest=[ satrec.ecco; satrec.inclo; satrec.nodeo; satrec.argpo; satrec.mo;str2double(tle2(52:63))]';

% Apriori cost
fprintf('initial cost = %f',computeCostNew(Xest));
%% Optimization function
tic
lb=[0.00001, -Inf, -Inf, -Inf, -Inf,  0.001]; % lower bound
ub=[0.99999,  Inf,  Inf,  Inf,  Inf,   125]; % upper bound
MaxIterations=1000;
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
Xest = constraintCOE(Xest);
[tle1,tle2]=coe2tle(tle1,tle2,Xest(1),Xest(2),Xest(3),Xest(4),Xest(5),Xest(6));
fprintf(2,'%s\n',tle2);
fprintf(tle_true)

[RA_comp,DEC_comp] = getComputedRADEC(Xest);

figure
subplot(2,1,1)
plot(180/pi*(RA_comp-RA_obs),'r*');
subplot(2,1,2)
plot(180/pi*(DEC_comp-DEC_obs),'r*');

%% some own functions
function [rho] = getRhoVector(rtasc,dec)
% retruns a vector from RA and DEC
cd=cos(dec);
rho=[cos(rtasc).*cd, cd.*sin(rtasc), sin(dec)];
end

% constraint the orbital elements
function Xest = constraintCOE(Xest)
if Xest(1) >=1
    Xest(1) = 0.99999; %ecc
elseif Xest(1) <=0
    Xest(1)=0.00001;
end
Xest(2)=rem(Xest(2),2*pi);  % incl
Xest(3)=rem(Xest(3),2*pi);  % argp
Xest(4)=rem(Xest(4),2*pi);  % omega
Xest(5)=rem(Xest(5),2*pi);  % omega
if Xest(6) >=25
    Xest(6) = 25;
elseif Xest(6) <=0
    Xest(6)=0.0000001;
end
end
