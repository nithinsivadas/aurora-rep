 function mtime = unixtime2matlab(x)
%% unixtime2matlab.m converts unix time to matlab time
mtime = x/(24*3600) + datenum('jan-01-1970');