%% GenerateA.m
% Forward model to estimate the matrix A that relates
% ion production rate (q) with number flux of precipitating electrons (phi)
% 
% q = A * phi
%
%
clear;

% Loading data from the HDF5 file created by a Glow run of precipitating electrons
fileName  ='/home/nithin/Dropbox/aurora_data (1)/eigenprofiles/GLOWrates.h5';
% fileName  ='q_2008_03_26.h5';
q 		  = hdf5read(fileName,'prod/eigenprofile');	%Production rates (Nreactionz x Nalt x NEnergy)
% The reactions are as follows
% 0  O+(2P)
% 1  O+(2D)
% 2  O+(4S)
% 3  N+
% 4  N2+
% 5  O2+
% 6  NO+
% 7  N2(A)
% 8  N(2P)
% 9  N(2D)
% 10 O(1S)
% 11 O(1D)

Ebins	  = hdf5read(fileName,'Ebins');				%Energy bin values
z         = hdf5read(fileName,'altitude');			%Altitude
diffflux  = hdf5read(fileName,'diffnumflux');		%Differential number flux (cm^-2 s^-1 eV^-1)
sensorloc = hdf5read(fileName,'sensorloc'); 		%Unit in degrees, Geographical coordinates
Edim	  = 5:1:55;									%Energy of interest
Edim0     = 1:1:60;									%Energy available for computation
zdim      = 31:1:102;								%Altitude of interest
zdim0     = 1:1:length(z);							%Altitude available for computation
%[phi]     = diag(diffflux(Edim0).*(-Ebins(Edim0)+Ebins(Edim0+ones(1,length(Edim0)))));
[phi]     = diag(diffflux(Edim0));					%Total Flux
phi_inv   = inv(phi);

[Z,X]=meshgrid(z(zdim),Ebins(Edim)/(1000));			%For surface plots

% Finding A
A = squeeze(sum(q(:,zdim0,Edim0)))*(phi_inv); 

%Plotting A
figure; surf(X,Z,log10(A(zdim,Edim))','EdgeColor','none');
set(get(gca,'XLabel'),'String','Energy [keV]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('A Matrix');

%Plotting q
figure;surf(X,Z,log10(squeeze(sum(q(:,zdim,Edim))))','EdgeColor','none');
set(get(gca,'XLabel'),'String','Energy [keV]');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Production rate');

%Calculating q from A for particular energies
phi_test=zeros(length(Edim),1);
index=5:10:45;
for i=1:length(index)
	phi_temp=phi_test;
	phi_temp(index(i)-4)=diffflux(index(i));
	q_test(i,:)=A(zdim,Edim)*phi_temp;
    q_glow(i,:)=squeeze(sum(q(:,zdim,index(i))));
end;

%Plotting production rates estimated from the A evaluated above
%comparison with the production rates from glow
figure;
a1=plot(q_test,z(zdim));
legend(strcat(num2str(Ebins(index)/1000,'%10.2d\n'),' keV'));
title('Production rates estimated using A');
set(get(gca,'XLabel'),'String','Production Rate');
set(get(gca,'YLabel'),'String','Altitude (in km)');
figure;
a2=plot(q_glow,z(zdim));
legend(strcat(num2str(Ebins(index)/1000,'%10.2d\n'),' keV'));
title('Production rates from Glow');
set(get(gca,'XLabel'),'String','Production Rate');
set(get(gca,'YLabel'),'String','Altitude (in km)');

