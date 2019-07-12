close; clear; clc;
addpath("celestrack");
ob=load('obs_utc.txt');
yr=ob(:,1);
mon=ob(:,2);
day=ob(:,3);
hr=ob(:,4);
min=ob(:,5);
sec=ob(:,6);
for i =1:length(ob)
[jd, jdfrac] = jday(yr(i), mon(i), day(i), hr(i), min(i), sec(i));
JDN(i,1)=jd;
JDF(i,1)=jdfrac;
end