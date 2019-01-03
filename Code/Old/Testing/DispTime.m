function time_string = DispTime(t)
%Function DispTime(T) to display
if t >= 1800
    time_string = datestr(datenum(0,0,0,0,0,t),'MM:SS');
elseif t >= 3600
    time_string = datestr(datenum(0,0,0,0,0,t),'HH:MM:SS');
else
    time_string = datestr(datenum(0,0,0,0,0,t),'SS');
end
end