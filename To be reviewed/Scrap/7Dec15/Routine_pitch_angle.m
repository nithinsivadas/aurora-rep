%% Routine to Plot Normalized pitch angle distribution
%  THEMIS D - SST  Pitch angle distribution using thd.mat derived from
%  SPEDAS

%  First half generates a plot

%  Second half generates a movie

%% Loading the File
load('/home/nithin/Documents/git-repos/aurora-rep/To be reviewed/Scrap/2Dec15/thd.mat');

%% Plotting the Normalized Pitch Angle Distribution

% Substorm Growth 
Date(1)=datenum('26-Mar-2008 09:00');
Date(2)=datenum('26-Mar-2008 10:00');
Date(3)=datenum('26-Mar-2008 11:00)');

% Substorm Growth - crossing the PSBL
Date(4)=datenum('26-Mar-2008 11:01');
Date(5)=datenum('26-Mar-2008 11:02');
Date(6)=datenum('26-Mar-2008 11:03');
Date(7)=datenum('26-Mar-2008 11:10');

% Substorm Onset and Expansion
Date(8)=datenum('26-Mar-2008 11:55');
Date(9)=datenum('26-Mar-2008 12:03');
Date(10)=datenum('26-Mar-2008 12:30');

paefluxsst=interp1(thdsst.time,thdsst.paeflux,Date,'linear');
flux=sum(paefluxsst,2);
paefluxsst=paefluxsst./repmat(max(paefluxsst,[],2),1,length(thdsst.pa));

figure;
subplot(1,3,1)
n=1:3;
plot(thdsst.pa,paefluxsst(n,:),'.-');
legend1=datestr(Date(n),'HH:MM');
legend(legend1);
xlim([0 180]);
title('Growth');
set(get(gca,'YLabel'),'String','Normalized Pitch Angle Distribution');

subplot(1,3,2)
n=[4:7];
plot(thdsst.pa,paefluxsst(n,:),'.-');
legend1=datestr(Date(n),'HH:MM');
legend(legend1);
xlim([0 180]);
title('Growth - PSBL & South Lobe');
set(get(gca,'XLabel'),'String','Pitch Angle [deg] ');

subplot(1,3,3)
n=8:10;
plot(thdsst.pa,paefluxsst(n,:),'.-');
legend1=datestr(Date(n),'HH:MM');
legend(legend1);
xlim([0 180]);
title('Onset & Expansion');


%% Generating Movie
datetime1=datenum('26-03-2008 08:01','dd-mm-yyyy HH:MM');
datetime2=datenum('26-03-2008 12:59','dd-mm-yyyy HH:MM');
Date=linspace(datetime1,datetime2,1000);
paefluxsst=interp1(thdsst.time,thdsst.paeflux,Date,'linear');
flux=sum(paefluxsst,2);
paefluxsst=paefluxsst./repmat(max(paefluxsst,[],2),1,length(thdsst.pa));

k=1;
figure, set(gcf, 'Color','white'); axis tight;
set(gca, 'nextplot','replacechildren', 'Visible','off'); 
for datetime=linspace(datetime1,datetime2,1000);
    plot(thdsst.pa,paefluxsst(k,:),'.-');
    xlim([0 180]);
    % Indicating reliable and unreliable fluxes
    if flux(k)<10000
        str1='Unreliable flux';
    else
        str1=['Total electron energy flux = ',num2str(round(flux(k)/1000)),' [keV/ cm^2 sec sr eV]'];
    end
    title({[datestr(Date(k),'HH:MM'),' U.T.'],str1});
    ylim([0 1]);
    set(get(gca,'YLabel'),'String','Normalized Pitch Angle Distribution');
    set(get(gca,'XLabel'),'String','Pitch Angle [deg] ');
    grid on;
    M(k)=getframe(gcf); k=k+1;
    clf;
 end
movie2avi(M(2:end),'normalized_pa.avi');



