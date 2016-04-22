% A script to plot the production rates 
% of an energy spectrum (with a particular pitch angle range)
% and compare them with the production rates derived from PFISR data

% by Nithin Sivadas
% Dec 4th, 2015

% Same as Routine 3, but with pitch angle < 22.5 deg
clear all;

load('/home/nithin/Documents/git-repos/aurora-rep/To be reviewed/Scrap/2Dec15/thd.mat');
DateNumBeg=datenum('26-Mar-2008 08:00');
DateNumEnd=datenum('26-Mar-2008 13:00');

%% THEMIS D SST
emin_bin=16;

eb=thdesa.ebin(emin_bin:1:end);
eb(length(eb)+1:1:length(eb)+length(thdsst.ebin))=thdsst.ebin;

thdeflux = interp1(thdesa.time,thdesa.eflux(:,emin_bin:1:end),thdsst.time);
thdeflux(:,size(thdeflux,2)+1:1:size(thdeflux,2)+size(thdsst.eflux,2))=thdsst.eflux;

% Description of data available in the structure file
% thdsst.time  - Time in matlab units
% thdeflux     - Combined Energy Flux vales from esa and sst in [ev/ cm^2 sec sr eV]
% eb           - Combined energy values per bin from esa and sst [eV]

%% Calculating A from Rees Model and Inverting q to Flux in the Ionosphere

fileName{1} ='A thd E Normal 60to185km.h5';
data=convertA2E_v2(char(fileName{1}));

thdeflux(isnan(thdeflux)) = 10^-1 ;

% Calculating produxtion rates from the THEMIS D Spectrum (without
% considering the Auroral Acceleration Region
q_thd_N=zeros(length(thdsst.time),length(data.z));

for i=1:1:length(thdsst.time)
    q_thd_N(i,:)=(data.A(:,:))*(thdeflux(i,:)'./eb(:)); % Production rate derived from the THEMIS flux observed
end;

q_thd_N(q_thd_N==0)=10^-5;






%% Calculating low and high energy flux 
timerange=920:1:1090;
ebin_yi  = interp1(1:1:size(eb,2),eb,0.5:1:(size(eb,2)+0.5),'linear','extrap');
ebin_high= ebin_yi(2:end);
ebin_low = ebin_yi(1:end-1);
eb_binsize=ebin_high-ebin_low;
eb_ln=1:20;
eb_hn=21:length(eb);
thdloweflux=sum(thdeflux(timerange,eb_ln).*repmat(eb_binsize(eb_ln),length(timerange),1),2);
thdhigheflux=sum(thdeflux(timerange,eb_hn).*repmat(eb_binsize(eb_hn),length(timerange),1),2);

% figure;plot(thdsst.time(timerange),10*log10(thdloweflux/max(thdloweflux))); hold on; plot(thdsst.time(timerange), 10*log10(thdhigheflux/max(thdhigheflux)));datetick('x','HH:MM');
timerange1=240:1:375;
data.eflux1=data.eflux';
pfisrloweflux =sum(data.eflux1(timerange1,eb_ln).*repmat(eb_binsize(eb_ln),length(timerange1),1),2);
pfisrhigheflux=sum(data.eflux1(timerange1,eb_hn).*repmat(eb_binsize(eb_hn),length(timerange1),1),2);

figure;
subplot(2,1,1);
plot(data.time_q(timerange1),10*log10(pfisrloweflux/max(pfisrloweflux))); hold on; plot(data.time_q(timerange1), 10*log10(pfisrhigheflux/max(pfisrhigheflux)));datetick('x','HH:MM');
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Ratio of Electron Energy Flux to the Maximum Electron Energy Flux [dB]');
legend('PFISR LOW 0.36 - 81.4 keV','PFISR HIGH 81.4 - 894 keV');

subplot(2,1,2);
plot(data.time_q(timerange1),(pfisrloweflux)); hold on; plot(data.time_q(timerange1), (pfisrhigheflux));
hold on; plot(thdsst.time(timerange),(thdloweflux),'.-'); hold on; plot(thdsst.time(timerange), (thdhigheflux),'.-');datetick('x','HH:MM');
set(gca,'YScale','log');
legend('PFISR LOW 0.36 - 81.4 keV','PFISR HIGH 81.4 - 894 keV', 'THEMIS D LOW 0.36 - 81.4 keV','THEMIS D HIGH 81.4 - 894 keV');
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Integrated Electron Energy Flux [eV/cm^2 sr sec]');

%% Plots the way I like it!
cmin=0;             cmax=6;
minz=60;            maxz=185;
minx=DateNumBeg;    maxx=DateNumEnd;

%% Prodiction rates
figure;
% Production rate plot from PFISR Measuremets
[T1,Z]=meshgrid(data.time_q,data.z);	
subplot(3,1,1);
surf(T1,Z,real(log10(data.q_PFISR)),'EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rates derived from PFISR Electron Density in log10 [cm^-^3 s^-^1]');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax]);

% Production rate from THEMIS flux derived from PFISR Measurements
[T1,Z]=meshgrid(data.time_q,data.z);	
subplot(3,1,2);
surf(T1,Z,real(log10(data.q_THEMIS)),'EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rates from THEMIS flux derived from PFISR measurements in log10 [cm^-^3 s^-^1]');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax])

% Production rate plot from THEMIS Flux Measuremets
[T1,Z]=meshgrid(thdsst.time,data.z);	
subplot(3,1,3);
surf(T1,Z,real(log10(q_thd_N))','EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rates derived from the THEMIS-D electron energy flux (p.a <22.5 deg) in log10 [cm^-^3 s^-^1]');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax])

%% Electron Density and Reduced Chi-squared of Fits
cmin=6;             cmax=12;
minz=60;            maxz=185;

figure;
% Production rate plot from THEMIS Flux Measuremets
subplot(2,1,1)
[T1,Z]=meshgrid(data.time_ne,data.z);
surf(T1,Z,real(log10(data.ne)),'EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Electron density measured using PFISR in log10 [m^-^3]');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax]);

% Chi sqaured
subplot(2,1,2)
plot(data.time_q,data.chi);
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','\chi^2');
title('The \chi^2 value of Maximum Entropy Inversion - Fit');
xlim([DateNumBeg DateNumEnd]);



%% Kernel A 
cmin=-10;             cmax=4;
minz=60;            maxz=185;
minx=data.Ebin(1);  maxx=data.Ebin(end);

figure;
[Z,E]=meshgrid(data.z,data.Ebin);
surf(E,Z,real(log10(data.A))','EdgeColor','none');
set(get(gca,'XLabel'),'String','Energy [eV]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rate per Energy Flux [cm^-^1 Sr eV/ eV]');
view([0 90]);
set(gca,'YScale','log');
set(gca,'XScale','log');
xlim([minx maxx]);
ylim([minz maxz]);
caxis([cmin cmax]);

%% Create Emperical A
datetime=datenum('26-03-2008 09:45','dd-mm-yyyy HH:MM');
thdphi=interp1(data.time_q,data.eflux(:,:)',datetime)'./eb(:);
% thdphi=interp1(thdsst.time,thdeflux(:,:),datetime)'./eb(:);
dataq=interp1(data.time_q,data.q_PFISR(:,:)',datetime)';
% dataq=interp1(data.time_q,data.eflux(:,:)',datetime)'./eb(:);
A_emp = (dataq*thdphi')*(inv(thdphi*thdphi'+10^-10));

cmin=-10;             cmax=4;
minz=60;            maxz=185;
minx=data.Ebin(1);  maxx=data.Ebin(end);
figure;
[Z,E]=meshgrid(data.z,data.Ebin);
surf(E,Z,real(log10(A_emp))','EdgeColor','none');
set(get(gca,'XLabel'),'String','Energy [eV]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rate per Energy Flux [cm^-^1 Sr eV/ eV]');
view([0 90]);
set(gca,'YScale','log');
set(gca,'XScale','log');
xlim([minx maxx]);
ylim([minz maxz]);
% caxis([cmin cmax]);

%% Emperical A Kernal for phi1 and phi2
datetime=datenum('26-03-2008 09:45','dd-mm-yyyy HH:MM');
thdphi=interp1(thdsst.time,thdeflux(:,:),datetime)'./eb(:);
% dataq=interp1(data.time_q,data.q_PFISR(:,:)',datetime)';
dataq=interp1(data.time_q,data.eflux(:,:)',datetime)'./eb(:);
A_emp = (dataq*thdphi')*(inv(thdphi*thdphi'+10^-10));

minx=eb(1);            maxx=eb(end);
minz=data.Ebin(1);     maxz=data.Ebin(end);
figure;
[Z,E]=meshgrid(eb,data.Ebin);
surf(E,Z,real(log10(A_emp))','EdgeColor','none');
set(get(gca,'XLabel'),'String','Energy [eV] at the plasmasheet from THEMIS');
set(get(gca,'YLabel'),'String','Energy [eV] at Ionosphere from PFISR');
title('Ratio of Differential Number Flux between Plasmasheet and the Ionosphere');
view([0 90]);
set(gca,'YScale','log');
set(gca,'XScale','log');
xlim([minx maxx]);
ylim([minz maxz]);
% caxis([cmin cmax]);


%% Differential Electron Energy Flux
cmin=2;             cmax=9;
minz=eb(1);            maxz=eb(end);

figure;
% THEMIS Energy Flux measured
[T1,Z]=meshgrid(thdsst.time,eb);	
subplot(2,1,1);
surf(T1,Z,real(log10((thdeflux+10^-6))'),'EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Energy [eV]');
title('Electron Energy Flux from THEMIS D in log10 [eV/ cm^2 sec sr eV] for pitch angle <22.5 deg');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax])

% THEMIS Energy Flux derived from PFISR Measurements
[T1,Z]=meshgrid(data.time_q,data.Ebin);	
subplot(2,1,2);
surf(T1,Z,real(log10((data.eflux))),'EdgeColor','none');
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Energy [eV]');
title('Electron Energy Flux derived from PFISR Measurments in log10 [eV/ cm^2 sec sr eV] for pitch angle <22.5 deg');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax]);
%% Movie
datetime1=datenum('26-03-2008 08:01','dd-mm-yyyy HH:MM');
datetime2=datenum('26-03-2008 12:59','dd-mm-yyyy HH:MM');

k=1;
figure, set(gcf, 'Color','white'); axis tight;
set(gca, 'nextplot','replacechildren', 'Visible','off'); 
for datetime=linspace(datetime1,datetime2,500);
    thdeflux1=interp1(thdsst.time,thdeflux(:,:),datetime);
    dataeflux1=interp1(data.time_q,data.eflux(:,:)',datetime);
    plot(eb,thdeflux1,'.-'); 
    hold on; 
    plot(data.Ebin,dataeflux1);
    set(gca,'YScale','log');
    set(gca,'XScale','log');
    xlim([10^2 10^7]);
    ylim([10^2 10^9]);
    title([datestr(datetime),' U.T.']);
    set(get(gca,'XLabel'),'String','Energy [eV]');
    set(get(gca,'YLabel'),'String','\phi [eV cm^-^2 sec^-^1 sr^-^1 eV^-^1] ');
    grid on;
    legend('original (p.a. <1.5 deg)','inverted');
    M(k)=getframe(gcf); k=k+1;
    clf;
 end
fps=10;
movie2avi(M(2:end),'pa_1_5.avi');
