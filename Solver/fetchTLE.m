%% search and fetch two line element string from text file
% INPUT
% search_query - first few characters of satellite name
% filename     - TLE filename e.g.
% OUTPUT
% name - satellite name
% tle1 - two line element string 1
% tle2 - two line element string 2
% example:
% [name,tle1,tle2] = fetchTLE('ISS','active.txt');
% siddharth@deore.in

function [name,tle1,tle2] = fetchTLE(search_query,filename)
if(nargin<2)
    % if filename not specified, search in active.txt
    filename='active.txt';
end
 file = dir(strcat('TLE\',filename));
 if(hours(datetime('now')-file.date)>24)
     %cprintf('*red','WARNING : file is 24+ hrs old, please update\n')
     fprintf(2,'WARNING : file is 24+ hrs old, please update\n')
 end

 search_query=upper(search_query);
 fprintf(2,'Searching "%s" in %s \n',search_query,filename);
% cprintf('blue','Searching "%s" in %s \n',search_query,filename)
 % read text file delimited with '\n' as cell array
 cell_array=textread(strcat('TLE\',filename),'%s','delimiter','\n');
 % find the line number of first occurence in search string
 line_index=find(~cellfun(@isempty,strfind(cell_array,search_query)));
 % if search string exists in file then return tle1 and tle2 else return
 % first tle in file
 if(length(line_index)>0)
    name=char(cell_array(line_index(1)));
    tle1=char(cell_array(line_index(1)+1));
    tle2=char(cell_array(line_index(1)+2));
    % cprintf([0,0.5,0],'TLE of Object '); cprintf(-[0,0.5,0],'%s ',strtrim(name)); cprintf([0,0.5,0],' fetched sucessfully!\n');
    
    fprintf('TLE of Object ');
    fprintf(2,'%s ',strtrim(name));
    fprintf(' fetched sucessfully!\n');
else
    name=char(cell_array(1));
    tle1=char(cell_array(1+1));
    tle2=char(cell_array(1+2));
    cprintf('*red','No results found for searched string !!\n',name);
    cprintf('*blue','Loading default TLE of Satellite "%s" \n',name);
end
end