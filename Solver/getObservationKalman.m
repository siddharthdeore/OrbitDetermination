%%
function [y_obs] = getObservationKalman(state)
global terms min_in_day timezone  stationLat stationLon stationAlt xp yp dut1 lod ddpsi ddeps dat  current_time
        a=[0, 0, 0];
    [~, ~, ~, satrec] = twoline2rvDouble1(state);
        epoch_jd=satrec.jdsatepoch + satrec.jdsatepochf;

    %theta = RA_comp;
        [satrec, r, v] = sgp4(satrec,current_time);	% r v is in TEME
        [year,mon,day,hr,min,sec] = invjday ( epoch_jd,(current_time)/min_in_day);
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
                        
        %rho_comp(i,:)=getRhoVector(RA_comp(i,:),DEC_comp(i,:));
    y_obs=[rtasc;dec];
end
