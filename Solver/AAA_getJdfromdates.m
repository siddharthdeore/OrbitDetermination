clc; clear; close all;
addpath("../common","vallado");

dts=load('dates.txt');
day = dts(:,1);
mon = dts(:,2);
year = dts(:,3);
hr = dts(:,4);
min  = dts(:,5);
sec = dts(:,6) + 0.5;
format longG
for i = 1:length(dts)
[jd,jdfrac] = jday(year(i), mon(i), day(i), hr(i), min(i), sec(i));
%[ayear,amon,aday,ahr,amin,asec] = invjday ( jd, jdfrac );
JD(i)=jd;
JDF(i)=jdfrac;
JDN(i,1)=jd+jdfrac;
end
