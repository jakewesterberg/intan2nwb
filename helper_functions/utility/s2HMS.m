%SECS2HMS - converts a time in seconds to a string giving the time in hours, minutes and second
%Usage TIMESTRING = SECS2HMS(TIME)]);
%Example 1: >> secs2hms(7261)
%>> ans = 2 hours, 1 min, 1.0 sec
%Example 2: >> tic; pause(61); disp(['program took ' secs2hms(toc)]);
%>> program took 1 min, 1.0 secs

function time_string=s2HMS(time_in_secs)
time_string='';
nhours = 0;
nmins = 0;
nhours = floor(time_in_secs/3600);
hour_string = ':';
time_string = [num2str(nhours) hour_string];
nmins = floor((time_in_secs - 3600*nhours)/60);
minute_string = ':';
time_string = [time_string num2str(nmins) minute_string];
nsecs = time_in_secs - 3600*nhours - 60*nmins;
time_string = [time_string sprintf('%2.1f', nsecs) ''];
end