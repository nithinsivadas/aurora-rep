%% Script file: Filtermadrigal1mar2011.m
%  This script file downloads the required PFISR .hdf5 data file from madrigal
%  and processes it to plot average elctron desity varying with only range and time.
%% Output: A plot of PFISR electron density vs. Altitude

clear;

% Defining the URL of the Madrigal Database
cgiurl='http://isr.sri.com/madrigal/cgi-bin/';

% Collecting the experiment arrays that satisfy the following data and time
% criterion
expArray = getExperimentsWeb(cgiurl, 61, datenum('03/01/2011 08:00:00'),...
    datenum('03/01/2011 12:00:00'), 1);

% Collecting the file names
expFileArray = getExperimentFilesWeb(cgiurl,expArray(1).id);


% Listing down the parameters in the files
parmArray = getParametersWeb(cgiurl, expFileArray(1).name);

% Downloading the first file in the expFileArray
result=madDownloadFile(cgiurl, expFileArray(1).name,'test01','Nithin Sivadas','nithin@bu.edu','Boston University','hdf5');

% Using GeoData class made by John to extract range and density
pfisrGD = GeoData(@readMadhdf5,'test01',{'range','popl'});


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
popl_sum=zeros(range_data_l(1)/nb-1,time_data_l(2)); % Sum of density along each range slice
total_electron=zeros(range_data_l(1)/nb-1,time_data_l(2)); % Average of density per range slice
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

% Altitude
altitude_full=1:1:range_data_l(1)/nb - 1;
altitude_low=1:1:50;

% Plot the total electron density in the entire volume
%x=datenum(1970,1,1,0,0,0)+pfisrGD.times(:,1)/(60*60*24);
x=unixtime2matlab(pfisrGD.times(:,1));
figure;
plot(x',sum(total_electron(altitude_full,:))); datetick('x', 'HH:MM');
%figure;
%plot(x',sum(total_electron(altitude_low,:)),'g'); datetick('x', 'HH:MM');

% Total high energy electron population with time (i.e. low altitude)
sum_low_altitude_electron=sum(total_electron(altitude_low,:));

% Plotting the average electron density per range slice
figure;  imagesc(x,altitude,log10(total_electron)); datetick('x', 'HH:MM');

