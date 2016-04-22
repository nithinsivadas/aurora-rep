%Routine to Plot the Electron Density vs. Altitude
clear;

DateNumBeg=datenum('16-Feb-2008 08:00');
DateNumEnd=datenum('16-Feb-2008 13:00');
agrid=[70:1:200]';
cmin=0;             cmax=6;
minz=70;            maxz=200;
minx=DateNumBeg;    maxx=DateNumEnd;

fileNameStr ='DataFile_2008_16_Feb.h5';
% altitude_grid=65:1:220;
[total_electron, altitude_grid, t ] = pfisr_avg_e_v2(fileNameStr,agrid);

[T1,Z]=meshgrid(t,altitude_grid);	
surf(T1,Z,real(log10(total_electron)),'EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('PFISR Electron Density in log10 [cm^-^3 s^-^1] on 16th Feb 2008');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
% caxis([cmin cmax]);
