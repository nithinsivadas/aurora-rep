%% Plotting data from filtered experiments that satisfy necessary criteria
% A lower version of madrigalPlotSubstorm.m
% Input variable: expFileArrayStore08_11.mat [Contains all the Baker Code PFISR range and popl data]
 %                ae_mod_08_11.mat  [Contains the auroral electrojet index
 %                                   from 2008-2011]
%% Output: A plot of PFISR electron density vs. Altitude with time AND AE index vs. time
% Must incorporate the sunrise time here 
% Identify how to look up the year so that the correct AL index can be used

load /home/nithin/NithinBU/PFISR/expFileArratStore08_11.mat
load /home/nithin/NithinBU/PFISR/AURORAL_ELECTROJET/ONE_MINUTE/ae_mod_08_11.mat

cgiurl='http://isr.sri.com/madrigal/cgi-bin/';
[d1, d2] = size(expFileArrayStore);


for i=1:1:d2

	%% Reading particluar files     
	%pfisrGD = GeoData(@readMadhdf5,fileNameStr,{'range','popl'});

	%theta_min = min(pfisrGD.dataloc(:,3));
	%theta_max = max(pfisrGD.dataloc(:,3));

	%phi_min = min(pfisrGD.dataloc(:,2));
	%phi_max = max(pfisrGD.dataloc(:,2));

	% Identifying the number of beams
	%temp = pfisrGD.dataloc(1,1);
	%i=1;

	%while (pfisrGD.dataloc(i+1,1)==temp)
	%    i=i+1;
	%end;    
	%nb=i;

	%range_data_l = size(pfisrGD.dataloc);
	%time_data_l= size(pfisrGD.data.popl);

	%%Summing up the electrons across latitude and longitude

	%popl_sum=zeros(range_data_l(1)/nb-1,time_data_l(2));
	%total_electron=zeros(range_data_l(1)/nb-1,time_data_l(2));
	%altitude=zeros(range_data_l(1)/nb-1,1);
	%
	%for i=1:1:(range_data_l(1)/nb - 1)
	%popl_sum(i,:)=popl_sum(i,:)+nansum(10.^pfisrGD.data.popl((i-1)*nb+1:1:i*nb,1:1:time_data_l(2)));    
	%popl_mean(i,:)=popl_sum(i,:)./(nb*ones(1,time_data_l(2))-(sum(isnan(pfisrGD.data.popl((i-1)*nb+1:1:i*nb,1:1:time_data_l(2))))));
	%r1=pfisrGD.dataloc((i-1)*nb+1,1);
	%%r2=pfisrGD.dataloc(i*nb+1,1);
	%altitude(i)=r1;
	%%dv= abs((r2^3-r1^3)*(cos(theta_max)-cos(theta_min))*(phi_max-phi_min)/3);
	%%total_electron(i,:)=total_electron(i,:)+popl_mean(i,:)*dv;
	%total_electron(i,:)=total_electron(i,:)+popl_mean(i,:);
	%end;

	%% Definig low Altitude section
	%altitude_low=1:1:50;

	%% Sum & Plot

	%% Matlab time of the event
	%x=unixtime2matlab(pfisrGD.times(:,1));

	%% Total high energy electron population with time
	%sum_low_altitude_electron=sum(total_electron(altitude_low,:));
	%
	%% Plotting the image of totale electron density vs. altitude vs
	%figure;  imagesc(x,altitude,log10(total_electron)); datetick('x', 'HH:MM');
	%colormap('jet');

	%% Listing down the parameters in the files\
	%parmArray = getParametersWeb(cgiurl, expFileArrayStore(i).name);

	% Collecting the necessary data
	% By downloading file
	fileNameStr=sprintf('DataFile_2008_%d',i);

	% Downloading the first file in the expFileArray
	result=madDownloadFile(cgiurl, expFileArrayStore(i).name,fileNameStr,'Nithin Sivadas','nithin@bu.edu','Boston University','hdf5');

	% Using GeoData class made by John to extract range and density
	pfisrGD = GeoData(@readMadhdf5,fileNameStr,{'range','popl'});

	% Finding the range of theta and phi within the data
	theta_min = min(pfisrGD.dataloc(:,3));
	theta_max = max(pfisrGD.dataloc(:,3));

	phi_min = min(pfisrGD.dataloc(:,2));
	phi_max = max(pfisrGD.dataloc(:,2));

	% Identifying the number of beams
	temp = pfisrGD.dataloc(1,1);
	i=1;

	while (pfisrGD.dataloc(i+1,1)==temp)
	    i=i+1;
	end;    
	nb=i;

	range_data_l = size(pfisrGD.dataloc);
	time_data_l= size(pfisrGD.data.popl);


	% Averaging electron desnity per range slice
		% Summing up the electrons across latitude and longitude
	popl_sum=zeros(range_data_l(1)/nb-1,time_data_l(2));
	total_electron=zeros(range_data_l(1)/nb-1,time_data_l(2));
	altitude=zeros(range_data_l(1)/nb-1,1);

	for i=1:1:(range_data_l(1)/nb - 1)
		popl_sum(i,:)=popl_sum(i,:)+nansum(10.^pfisrGD.data.popl((i-1)*nb+1:1:i*nb,1:1:time_data_l(2)));    
		popl_mean(i,:)=popl_sum(i,:)./(nb*ones(1,time_data_l(2))-(sum(isnan(pfisrGD.data.popl((i-1)*nb+1:1:i*nb,1:1:time_data_l(2))))));
		r1=pfisrGD.dataloc((i-1)*nb+1,1);
		%r2=pfisrGD.dataloc(i*nb+1,1);
		altitude(i)=r1;
		%dv= abs((r2^3-r1^3)*(cos(theta_max)-cos(theta_min))*(phi_max-phi_min)/3);
		%total_electron(i,:)=total_electron(i,:)+popl_mean(i,:)*dv;
		total_electron(i,:)=total_electron(i,:)+popl_mean(i,:);
	end;

	%Altitude
	altitude_low=1:1:50;

	% Plot the total electron density in the entire volume
	%x=datenum(1970,1,1,0,0,0)+pfisrGD.times(:,1)/(60*60*24);
	x=unixtime2matlab(pfisrGD.times(:,1));
	%figure;
	%plot(x',sum(total_electron(altitude_full,:))); datetick('x', 'HH:MM');
	%figure;
	%plot(x',sum(total_electron(altitude_low,:)),'g'); datetick('x', 'HH:MM');

	% Total high energy electron population with time
	sum_low_altitude_electron=sum(total_electron(altitude_low,:));

	% Plotting electron density along altitude as well as the AE index
	figure;  
		subplot(2,1,1);
		imagesc(x,altitude,log10(total_electron)); datetick('x', 'HH:MM');
		colormap('jet');
		subplot(2,1,2)
		plot(x,al_x,x,au_x); datetick('x', 'HH:MM');

	clearvars popl_sum popl_mean total_electron altitude;
end;