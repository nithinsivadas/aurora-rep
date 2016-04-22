function [] = plot_q(q_1,q_thd_N,nt1,nt2,altitude,time)
%PLOT_Q Simply plots the production rates vs. altitude at time instants
% represented by the time indices nt1 and nt2
%   Detailed explanation goes here

plot(log10(q_1(:,nt2)), altitude,'b'); hold on; 
plot(log10(q_thd_N(nt1,:)),altitude,'g'); 
xlim([0 7]);
title(datestr(time,'HH:MM'));

end

