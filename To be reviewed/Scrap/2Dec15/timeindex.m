function [ n_themis,n_pfisr ] = timeindex( time, thdtime, pfisrtime )
%TIMEINDEX Estimates the index of the time array corresponding to the
%input time 'time'
%   The outputs are integer numbers., 
%   n_themis - corresponds to the index of the thdtime array for value 'time'
%   n_pfisr  - corresponds to the index of the pfisrtime array for value 'time'

% Calculating the index of the THEMIS data time
[Data_X_Time, index]=unique(thdtime);
n_themis=floor(interp1(thdtime,(index),time));
% Calculaing the index of the PFISR data time
[Data_X_Time2, index2]=unique(pfisrtime);
n_pfisr=floor(interp1(pfisrtime,(index2),time));

end

