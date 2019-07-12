%% search and fetch two line element string from text file
% INPUT
% search_date -  string in format 'yyyy mm dd' e.g. '2019 02 23'
% OUTPUT
% EOP String
% yyyy dd mm MJD      x         y        UT1-UTC     LOD        dPsi    dEpsilon    dX        dY     DAT
% example:
% [name,tle1,tle2] = fetchTLE('ISS','active.txt');
% deore.in

function [xp,yp,dut1,lod,dPsi,dEpsilon,dx,dy,dat] = fetchSpaceData(search_date)
if(nargin<2)
    % if filename not specified, search in active.txt
    filename='EOP-Last5Years.txt';
end
 file = dir(strcat('SpaceData\',filename));
% if(hours(datetime('now')-file.date)>24)
%     %cprintf('*red','WARNING : file is 24+ hrs old, please update\n')
%     fprintf(2,'WARNING : file is 24+ hrs old, please update\n')
% end

 search_date=upper(search_date);
 fprintf(2,'Searching "%s" in %s \n',search_date,filename);
% cprintf('blue','Searching "%s" in %s \n',search_date,filename)
 % read text file delimited with '\n' as cell array
 cell_array=textread(strcat('SpaceData\',filename),'%s','delimiter','\n');
 % find the line number of first occurence in search string
 line_index=find(~cellfun(@isempty,strfind(cell_array,search_date)));
 % if search string exists in file then return tle1 and tle2 else return
 % first tle in file
 if(length(line_index)>0)
    EOP_String=char(cell_array(line_index(1)));
%fprintf('# ----------------------------------------------------------------------------------------------------\n');
%fprintf('#   Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT\n');
%fprintf('# (0h UTC)           "         "          s          s          "        "          "         "     s \n');
%fprintf('# ----------------------------------------------------------------------------------------------------\n');
%fprintf('# y4 mm dd nnnnn +n.nnnnnn +n.nnnnnn +n.nnnnnnn +n.nnnnnnn +n.nnnnnn +n.nnnnnn +n.nnnnnn +n.nnnnnn nnn\n');
%fprintf('# ----------------------------------------------------------------------------------------------------\n');
%fprintf('%s \n',strtrim(EOP_String));
   xp = str2num(EOP_String(18:26))* pi / (3600.0*180.0); % arc sec to rad
   yp = str2num(EOP_String(28:36))* pi / (3600.0*180.0); % arc sec to rad
   dut1 = str2num(EOP_String(38:47));
   lod = str2num(EOP_String(49:58));
   dPsi = str2num(EOP_String(60:68))* pi / (3600.0*180.0);
   dEpsilon = str2num(EOP_String(70:78))* pi / (3600.0*180.0);
   dx = str2num(EOP_String(80:88));%* pi / (3600.0*180.0);
   dy = str2num(EOP_String(90:98));%* pi / (3600.0*180.0);
   dat = str2num(EOP_String(100:102));

    %fprintf([0,0.5,0],' fetched sucessfully!\n');
else
%    EOP_Stirng=char(cell_array(1));
%    tle1=char(cell_array(1+1));
%    tle2=char(cell_array(1+2));
    %cprintf('*red','No results found for searched string !!\n',EOP_Stirng);
%    cprintf('*blue','Loading default TLE of Satellite "%s" \n',EOP_Stirng);
 end
 
end