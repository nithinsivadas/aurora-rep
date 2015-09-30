%% Plotting data of substorms from 2008 to 2011
%  Requires input data from PFISR and Auroral-Electrojet Index
%  The file plots the electron density observed using PFISR, and it's corresponding AE index
%  The code also marks sunset and sunrise times on the plots. 
%% Output: A plot of PFISR electron density vs. Altitude with time AND AE index vs. time with 
        %  sunrise and sunset markers

clearvars;

% Loading downloaded data
load '/home/nithin/NithinBU/Aug 3 2015/expFileArratStore08_14.mat'
load '/home/nithin/NithinBU/PFISR/AURORAL_ELECTROJET/ONE_MINUTE/ae_mod_08_14.mat'

cgiurl='http://isr.sri.com/madrigal/cgi-bin/';
[d1, d2] = size(expFileArrayStore);

for k=51:1:88

    % Guessing the filenames of downloaded data
    fileNameStr=sprintf('DataFile_%s_%d',expFileArrayStore(k).name(27:30),k);

    % Invoking the GeoData class to extract range and density from the HDF5 file
    pfisrGD = GeoData(@readMadhdf5,fileNameStr,{'range','popl'});

    % Identifying the range of theta and phi within the data
    theta_min = min(pfisrGD.dataloc(:,3));
    theta_max = max(pfisrGD.dataloc(:,3));

    phi_min = min(pfisrGD.dataloc(:,2));
    phi_max = max(pfisrGD.dataloc(:,2));

    % Identifying the number of beams
    temp = pfisrGD.dataloc(1,1);
    j=1;
    while (pfisrGD.dataloc(j+1,1)==temp)
        j=j+1;
    end;    
    nb=j;

    range_data_l = size(pfisrGD.dataloc);
    time_data_l= size(pfisrGD.data.popl);

    % Averaging electron desnity per range slice
    popl_sum=zeros(range_data_l(1)/nb-1,time_data_l(2));
    total_electron=zeros(range_data_l(1)/nb-1,time_data_l(2));
    altitude=zeros(range_data_l(1)/nb-1,1);

    % Plotting the Electron density with sunrise and sunset times
    %% Commented the following plotting script because it seems to be repeated below
    %figure('units','normalized','outerposition',[0 0 1 1],'Position',[0 0 1 1]);
    %fig=gcf;
    %set(fig,'PaperUnits','normalized','PaperOrientation','Landscape','PaperType','usletter','PaperPositionMode','Manual','PaperPosition',[-0.05 0 1.05 1],'Visible','off');
    %ax1=subplot(2,1,1); imagesc(x,altitude,log10(total_electron)); colormap('jet');datetick('x', 'HH:MM');
    %set(ax1,'XTick',[]);
    %colorbar('manual','position',[0.93 0.55 0.02 0.35]);
    %ax2=subplot(2,1,2); plot(x,al_x,x,au_x);datetick('x', 'HH:MM');
    %linkaxes([ax1 ax2],'x');
    %ax1.Position=[0.1 0.5 0.8 0.40];
    %set(ax2,'YAxisLocation','right','Position',[0.1 0.1 0.8 0.40],'XLim',[x(1) x(end)]);
    %% Possible error in plotting sunrise and sunset times here?
    %line([floor(x(5))+sunrise(4)/24 floor(x(5))+sunrise(4)/24],ax2.YLim,'parent',ax2,'linewidth',2,'color','black');
    %line([floor(x(5))+sunrise(4)/24 floor(x(5))+sunrise(4)/24],ax1.YLim,'parent',ax1,'linewidth',2,'color','black');
    %line([floor(x(5))+sunset(4)/24 floor(x(5))+sunset(4)/24],ax2.YLim,'parent',ax2,'linewidth',2,'color','m');
    %line([floor(x(5))+sunset(4)/24 floor(x(5))+sunset(4)/24],ax1.YLim,'parent',ax1,'linewidth',2,'color','m');
    %%text(floor(x(5))+(sunrise(4))/24, min(ax2.YLim)*1.13, datestr(floor(x(5))+(sunrise(4))/24,'HH:MM'));
    %%text(floor(x(5))+(sunrise(4))/24, min(ax2.YLim)*1.2, 'Sunrise Hour');
    %xlabel(ax2,'Time [UT]');
    %ylabel(ax1,'Altitude [Km]');
    %ylabel(ax2,'[nT]');
    %title(ax1,strcat('PFISR Electron density (popl) Date: ',datestr(floor(x(5)),'dd mmm yyyy'),' Altitude of Sun Hour - 150 km'));
    %legend1 = legend(ax2,'show','AL index','AU index', strcat('Sunrise Hour --',datestr(floor(x(5))+(sunrise(4))/24,' HH:MM')),strcat('Sunset Hour --',datestr(floor(x(5))+(sunset(4))/24,' HH:MM')));
    %set(legend1,'Location','southwest');
    %print(fig,strcat(num2str(k),'_',datestr(floor(x(5)),'dd_mmm_yyyy'),'.pdf'),'-dpsc2','-r300');
    %strcat('Exp_No_',num2str(k),'_Date_',datestr(floor(x(5)),'dd_mmm_yyyy'))

    for i=1:1:(range_data_l(1)/nb - 1)
        popl_sum(i,:)=popl_sum(i,:)+nansum(10.^pfisrGD.data.popl((i-1)*nb+1:1:i*nb,1:1:time_data_l(2)));    
        popl_mean(i,:)=popl_sum(i,:)./(nb*ones(1,time_data_l(2))-(sum(isnan(pfisrGD.data.popl((i-1)*nb+1:1:i*nb,1:1:time_data_l(2))))));
        r1=pfisrGD.dataloc((i-1)*nb+1,1);
        altitude(i)=r1;
        total_electron(i,:)=total_electron(i,:)+popl_mean(i,:);
    end;

    %Altitude
    altitude_low=1:1:50;

    % Sum & Plot
    x=unixtime2matlab(pfisrGD.times(:,1));

    % Identifying the year of aquisition of data    

    if condition
        body
    elseif condition
        body
    else
        body
    end

    if expFileArrayStore(k).name(27:30)=='2008'
        al_x=interp1(ae_mod_2008.al(:,1),ae_mod_2008.al(:,2),x);
        au_x=interp1(ae_mod_2008.au(:,1),ae_mod_2008.au(:,2),x);
    elseif expFileArrayStore(k).name(27:30)=='2009'
        al_x=interp1(ae_mod_2009.al(:,1),ae_mod_2009.al(:,2),x);
        au_x=interp1(ae_mod_2009.au(:,1),ae_mod_2009.au(:,2),x);
    elseif expFileArrayStore(k).name(27:30)=='2010'
        al_x=interp1(ae_mod_2010.al(:,1),ae_mod_2010.al(:,2),x);
        au_x=interp1(ae_mod_2010.au(:,1),ae_mod_2010.au(:,2),x);
    elseif expFileArrayStore(k).name(27:30)=='2011'
        al_x=interp1(ae_mod_2011.al(:,1),ae_mod_2011.al(:,2),x);
        au_x=interp1(ae_mod_2011.au(:,1),ae_mod_2011.au(:,2),x);
    elseif expFileArrayStore(k).name(27:30)=='2012'
        al_x=interp1(ae_mod_2012.al(:,1),ae_mod_2012.al(:,2),x);
        au_x=interp1(ae_mod_2012.au(:,1),ae_mod_2012.au(:,2),x);
    elseif expFileArrayStore(k).name(27:30)=='2013'
        al_x=interp1(ae_mod_2013.al(:,1),ae_mod_2013.al(:,2),x);
        au_x=interp1(ae_mod_2013.au(:,1),ae_mod_2013.au(:,2),x);
    elseif expFileArrayStore(k).name(27:30)=='2014'
        al_x=interp1(ae_mod_2014.al(:,1),ae_mod_2014.al(:,2),x);
        au_x=interp1(ae_mod_2014.au(:,1),ae_mod_2014.au(:,2),x);
    else
        error('Year beyond the scope of the AE indices for this version of the code. (Sorry).');
    end;

    % Total high energy electron population with time
    sum_low_altitude_electron=sum(total_electron(altitude_low,:));

    % Finding the sunrise and sunset hour
    sunrise = madCalculatorWeb(cgiurl,x(5),pfisrGD.sensorloc(1),pfisrGD.sensorloc(1),0,-360+pfisrGD.sensorloc(2),-360+pfisrGD.sensorloc(2),0,150,150,0,'SUNRISE_HOUR,SDWHT');
    if isnan(sunrise(4))==1
        if sunrise(5)<250
            sunrise(4)=0;
        else
            sunrise(4)=23.99;
        end;
    end;

    sunset = madCalculatorWeb(cgiurl,x(5),pfisrGD.sensorloc(1),pfisrGD.sensorloc(1),0,-360+pfisrGD.sensorloc(2),-360+pfisrGD.sensorloc(2),0,150,150,0,'SUNSET_HOUR,SDWHT');
    if isnan(sunset(4))==1
        if sunset(5)<250
            sunset(4)=23.99;
        else
            sunset(4)=0;
        end;
    end;

    % Plotting PFISR electron density and AE index with sunrise and sunset times
    figure('units','normalized','outerposition',[0 0 1 1],'Position',[0 0 1 1]);
    fig=gcf;
    set(fig,'PaperUnits','normalized','PaperOrientation','Landscape','PaperType','usletter','PaperPositionMode','Manual','PaperPosition',[-0.05 0 1.05 1],'Visible','off');
    ax1=subplot(2,1,1); imagesc(x,altitude,log10(total_electron)); colormap('jet');datetick('x', 'HH:MM');
    set(ax1,'XTick',[]);
    colorbar('manual','position',[0.93 0.55 0.02 0.35]);
    ax2=subplot(2,1,2); plot(x,al_x,x,au_x);datetick('x', 'HH:MM');
    linkaxes([ax1 ax2],'x');
    ax1.Position=[0.1 0.5 0.8 0.40];
    set(ax2,'YAxisLocation','right','Position',[0.1 0.1 0.8 0.40],'XLim',[x(1) x(end)]);
    line([floor(x(5))+sunrise(4)/24 floor(x(5))+sunrise(4)/24],ax2.YLim,'parent',ax2,'linewidth',2,'color','black');
    line([floor(x(5))+sunrise(4)/24 floor(x(5))+sunrise(4)/24],ax1.YLim,'parent',ax1,'linewidth',2,'color','black');
    line([floor(x(5))+sunset(4)/24 floor(x(5))+sunset(4)/24],ax2.YLim,'parent',ax2,'linewidth',2,'color','m');
    line([floor(x(5))+sunset(4)/24 floor(x(5))+sunset(4)/24],ax1.YLim,'parent',ax1,'linewidth',2,'color','m');
    %text(floor(x(5))+(sunrise(4))/24, min(ax2.YLim)*1.13, datestr(floor(x(5))+(sunrise(4))/24,'HH:MM'));
    %text(floor(x(5))+(sunrise(4))/24, min(ax2.YLim)*1.2, 'Sunrise Hour');
    xlabel(ax2,'Time [UT]');
    ylabel(ax1,'Altitude [Km]');
    ylabel(ax2,'[nT]');
    title(ax1,strcat('PFISR Electron density (popl) Date: ',datestr(floor(x(5)),'dd mmm yyyy'),' Altitude of Sun Hour - 150 km'));
    legend1 = legend(ax2,'show','AL index','AU index', strcat('Sunrise Hour --',datestr(floor(x(5))+(sunrise(4))/24,' HH:MM')),strcat('Sunset Hour --',datestr(floor(x(5))+(sunset(4))/24,' HH:MM')));
    set(legend1,'Location','southwest');
    print(fig,strcat(num2str(k),'_',datestr(floor(x(5)),'dd_mmm_yyyy'),'.pdf'),'-dpsc2','-r300');
    strcat('Exp_No_',num2str(k),'_Date_',datestr(floor(x(5)),'dd_mmm_yyyy'))

    clearvars popl_sum popl_mean total_electron altitude;

end;
