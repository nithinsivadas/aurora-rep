% A script to plot the production rates 
% of an energy spectrum (with a particular pitch angle range)
% and compare them with the production rates derived from PFISR data
% The following correction was made: eflux/ebin was multiplied by 2*pi.

% by Nithin Sivadas
% Dec 1st, 2015

clear all;

load('thd.mat');

addpath('/home/nithin/Documents/git-repos/energy-height-conversion/Rees_Ionization')
addpath('/home/nithin/Documents/aurora-rep-backup-Jan-16/aurora-rep/To be reviewed')
DateNumBeg=datenum('26-Mar-2008 08:00');
DateNumEnd=datenum('26-Mar-2008 13:00');

latitude    = 65; 								% [Deg]
longitude   = -148; 							% [Deg]

%% THEMIS D SST
% Calculating the lower and higher limits of the energy bins
ebin_yi  = interp1(1:1:size(thdsst.ebin,2),thdsst.ebin,0.5:1:(size(thdsst.ebin,2)+0.5),'linear','extrap');
ebin_high= ebin_yi(2:end);
ebin_low = ebin_yi(1:end-1);

% Description of data available in the structure file
% thdsst.time  - Time in matlab units
% thdsst.eflux - Energy Flux in [ev/ cm^2 sec sr eV]
% thdsst.ebin  - Energy values per bin [eV]

%% Calculating A from Rees Model

fileName{1} ='A thd E Normal.h5';
fileName{2} ='A thd E Low.h5';
fileName{3} ='A thd E High.h5';

for i=1:1:1
    q 		  = hdf5read(char(fileName{i}),'prod/eigenprofile');%Production rates (Nreactionz x Nalt x NEnergy)
    q         = q+min(min(q(q>0)));
    Ebins	  = hdf5read(char(fileName{i}),'Ebins');	% Energy bin values
    z         = hdf5read(char(fileName{i}),'altitude');	% Altitude
    % diffflux  = hdf5read(fileName,'diffnumflux');		% Differential number flux (cm^-2 s^-1 eV^-1)
    diffflux  = ones(size(Ebins,1),1);					% In Sergeinko Iranov model 
                                                        % diffnumber flux is 1 cm^-2 s^-1 eV^-1
    sensorloc = hdf5read(char(fileName{i}),'sensorloc');% Unit in degrees, Geographical coordinates
    Edim0     = 1:1:length(Ebins);						% Energy available for computation
    zdim0     = 1:1:length(z);							% Altitude available for computation
    [phi]     = diag(diffflux(Edim0));					% Total Flux
    phi_inv   = inv(phi);
    [Z,X]=meshgrid(z(zdim0),Ebins(Edim0)/(1000));		%For surface plots
    % Finding A
    % A = squeeze(sum(q(:,zdim0,Edim0)))*(phi_inv); 
    A(i,:,:) = (q(zdim0,Edim0))*(phi_inv); 
    %Plotting A
    figure; surf(X,Z,log10(squeeze(A(i,(zdim0-zdim0(1)+1),Edim0)))','EdgeColor','none');
    set(get(gca,'XLabel'),'String','Energy [keV]');
    set(get(gca,'YLabel'),'String','Altitude [km]');
    title(['Plot of A matrix with filename: ',char(fileName{i})]);
    view([0 90]);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    xlim([min(min(X)) max(max(X))]);
    ylim([min(min(Z)) max(max(Z))]);
end

%% Downloading the PFISR Data
% fileNameStr=downloadPFISR('2008-03-26 06:00','2008-03-26 15:00','yyyy-mm-dd HH:MM');
[ne,altitude,t]=pfisr_mag_e('/home/nithin/Documents/aurora-rep-backup-Jan-16/aurora-rep/To be reviewed/DataFile_2008_1.h5',1);
% ne is in m^-3 so converting to cm^-3
ne=ne*10^-6;

%% MSIS-90 Model

[D, T, F10_7_USED, AP_USED] = msis(DateNumBeg, latitude, longitude, altitude);

nO      = D(:,1);
nN2     = D(:,2);
nO2     = D(:,3);
density = D(:,4);
x=altitude;
for ty=1:1:size(ne,2)
    y=ne(:,ty);
    xi=x(find(~isnan(y)));
    yi=y(find(~isnan(y)));
    ne(:,ty)=real(interp1(xi,yi,x,'spline','extrap'));
end;
%% Using the AIDA Tools to get A matrix
[A1] = ionization_profile_matrix(altitude,thdsst.ebin',nO,nN2,nO2,density,1);

thdsst.eflux(isnan(thdsst.eflux)) = 10^-1 ;
for i=1:1:length(thdsst.time)
    diff_flux(1,:)=2*pi*thdsst.eflux(i,:)'./thdsst.ebin(:); % [cm-2 s-1 eV-1]
    q_thd_N(i,:)=(A1(:,:))*(diff_flux(1,:))';
end;

q_thd_N(q_thd_N==0)=10^-5;

[Data_X_Time, index]=unique(thdsst.time);
t1=floor(interp1(thdsst.time,(index),DateNumBeg));
t2=ceil(interp1(thdsst.time,(index),DateNumEnd));

%% Figure Plot q_thd
cmin=0;
cmax=6;
minz=50;
maxz=100;

[T1,Z]=meshgrid(thdsst.time,z(zdim0));	
figure; 
subplot(2,1,1);
surf(T1(:,(t1:1:t2)),Z(:,(t1:1:t2)),log10(q_thd_N((t1:1:t2),:))','EdgeColor','none');
set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rates derived from the THEMIS D Energy Spectrum with Normal Energy Bins in log10 [cm^-^3 s^-^1]');
view([0 90]);
set(gca,'YScale','log');
xlim([min(min(T1(:,(t1:1:t2)))) max(max(T1(:,(t1:1:t2))))]);
% ylim([min(min(Z(:,(t1:1:t2)))) max(max(Z(:,(t1:1:t2))))]);
ylim([minz maxz]);
caxis([cmin cmax])

%% Estimating q_1 - from PFISR measurements of electron density
alpha   = 2.5*10^-12*exp(-altitude(zdim0)./51.2); % [m^3/s] %z is in [km]
% Converting to cm^3/s
alpha   = alpha*10^6;
t1      = t(1:1:end-2);
t2      = t(3:1:end);
dt      = 86400*(t2-t1);
dn      = ne(zdim0,3:1:end)-ne(zdim0,1:1:end-2);
[DT,DNE]= meshgrid(dt,ne(zdim0,1));
dn_dt   = dn./(2*DT);

%phi_new = 10^13*ones(size(diag(phi)',2),1);   %[cm^-2 s^-1]
beta    = 20;
z_2     = z(zdim0);

% for i = 300:1:size(t1)
tdim = 1:1:size(t1);
q_1 = dn_dt(:,tdim) + diag(alpha)*(ne(zdim0,tdim).^2);

%% Figure Plot q_new
[T1,Z]=meshgrid(t2,z_2);	
subplot(2,1,2);
surf(T1,Z,real(log10(q_1)),'EdgeColor','none');
set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rates derived from PFISR Electron Density in log10 [cm^-3 s^-1]');
view([0 90]);
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
ylim([minz maxz]);
caxis([cmin cmax])