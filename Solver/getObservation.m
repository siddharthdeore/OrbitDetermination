%%
function [y_obs] = getObservation(state)
global terms min_in_day timezone  stationLat stationLon stationAlt xp yp dut1 lod ddpsi ddeps dat jd_obs t_obs delayTerm
        a=[0, 0, 0];
    [~, ~, ~, satrec] = twoline2rvDouble2(state); % without bstar
    
    epoch_jd=satrec.jdsatepoch + satrec.jdsatepochf;
    min_since_epoch = (jd_obs(1)-epoch_jd)*min_in_day+delayTerm;
    
    t_mins=(jd_obs-jd_obs(1))*min_in_day;
    %% Propogator
    NN=length(t_mins);
    %rho_comp=zeros(size(rho_obs));
    RA_comp = zeros(size(t_obs));
    DEC_comp = RA_comp;
    %theta = RA_comp;
    for i=1:NN
        [satrec, r, v] = sgp4(satrec, min_since_epoch+t_mins(i));	% r v is in TEME
        [year,mon,day,hr,min,sec] = invjday ( epoch_jd,(min_since_epoch+t_mins(i))/min_in_day);
        %% convert time
        [~, ~, jdut1, jdut1frac, ~, ~, ~, ttt, ~,~, ~, ~, ~,~, ~, ~,~, ~, ~,~ ] = convtime ( year, mon, day, hr, min, sec, timezone, dut1, dat );
        jd=jdut1+ jdut1frac;
        %% ECI
        [reci, veci, ~] = teme2eci  ( r', v', a', ttt, ddpsi, ddeps);
        %% Geodectic Lat Lon  -- Corrected
        [lst,~] = lstime ( stationLon, jd);  % local and Greanwitch sideral time
        %% Aximuth Elevation -- Corrected
        [~,az,el,~,~,~] = rv2razel ( reci,veci, stationLat,stationLon,stationAlt,ttt,jd,lod,xp,yp,terms,ddpsi,ddeps );
        az=mod(az,2*pi);
        %% RA DEC -- Corected
        [rtasc,dec] = azel2radec(az,el,stationLat,stationLon,lst);
        
        RA_comp(i,:)=rtasc;
        DEC_comp(i,:)=dec;
                        
        %rho_comp(i,:)=getRhoVector(RA_comp(i,:),DEC_comp(i,:));
    end
    y_obs=[RA_comp;DEC_comp];
end
%% some own functions
function [rho] = getRhoVector(rtasc,dec)
% retruns a vector from RA and DEC
rho=[cos(rtasc).*cos(dec), cos(dec).*sin(rtasc), sin(dec)];
end
