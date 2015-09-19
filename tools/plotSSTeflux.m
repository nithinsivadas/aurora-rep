function [ Data, DateNumBeg, DateNumEnd ] = plotSSTeflux( DateBeg, DateEnd, DateFormat,...
    Probes, eflux)
%% plotSSTeflux: 2D color plot of energy flux with time and energy value 
%% measured by the Solid State Telescope on-board the THEMIS spacecraft 
%   DateBeg    : String that specifies the beginning date of collecting data
%   (UT) '26 Mar 2008 8:00'
%   DateEnd    : String specifying end date '26 Mar 2008 8:00'
%   DateFormat : String specifying the format in which the date is written
%   'dd mm yyyy HH:MM'
%   Probes     : String specifying the probe within the mission {'thd','thc'}  
%   eflux      : String specifying the particular data set within a THEMIS Probe ('psif' or 'psef')  
% Example: 
% [Data, DateNumBeg, DateNumEnd]=plotSSTeflux( '26 Mar 2008 8:00', '26 Mar 2008 15:00', 'dd mm yyyy HH:MM', char('thd'), char('psef','psif'));

%   To be done
%   Comment it!

%%
Mission='themis';
DataLevel='l2';
Instrument='sst';
TargetPath=pwd;

DateNumBeg=datenum(DateBeg,DateFormat);
DateNumEnd=datenum(DateEnd,DateFormat);

ExpFileArray = downloadSPDF( DateBeg, DateEnd, DateFormat,...
    Mission, Probes, DataLevel, Instrument, TargetPath);

k=1;

for i=1:1:size(Probes)
    for j=1:1:size(eflux)
        b=1;a=1;
        for l=1:1:ceil(DateNumEnd)-floor(DateNumBeg)
        filename=[Probes(i,:),'_l2_sst_',datestr(floor(DateNumBeg+l-1),'yyyymmdd'),'_v01.cdf'];
        eflux_Y=[Probes(i,:),'_',eflux(j,:),'_en_eflux_yaxis'];
        eflux_var=[Probes(i,:),'_',eflux(j,:),'_en_eflux'];
        time=[Probes(i,:),'_',eflux(j,:),'_time'];
        [out]= spdfcdfread(filename,'Variables',{time,eflux_Y,eflux_var});
        Data(k).Probe=Probes(i,:);
        Data(k).eflux_name=eflux(j,:);
        a=a+size(out{1});
        Data(k).X_Time(b(1):1:a(1)-1)=unixtime2matlab(out{1});
        Data(k).eflux(b(1):1:a(1)-1,:)=out{3};
        Data(k).Y_Energy=out{2}(1,:);
        b=b+size(out{1});
        end;
        p=size(Data(k).Y_Energy);
        Data(k).Y_Energy=interp1((1:1:p(2)),Data(k).Y_Energy,(1:0.1:11));
        Data(k).eflux=(interp1((1:1:p(2)),Data(k).eflux',(1:0.1:11)'))';
        k=k+1;
        
    end;
end;

k=1;
n=size(eflux);
for i=1:1:size(Probes)
    figure('units','normalized','outerposition',[0 0 1 0.5],'Position',[0 0 1 1]);
    fig=gcf;
    set(fig,'PaperUnits','normalized','PaperOrientation','Landscape','PaperType','usletter','PaperPositionMode','Manual','PaperPosition',[-0.05 0 1.05 1],'Visible','on');
    for j=1:1:n(1)
        r=size(Data(k).X_Time);
        [Data_X_Time, index]=unique(Data(k).X_Time);
        t1=floor(interp1(Data_X_Time,(index),DateNumBeg));
        t2=ceil(interp1(Data_X_Time,(index),DateNumEnd));
        ax(j)=subplot(2,1,j); pcolor(Data(k).X_Time(t1:1:t2),Data(k).Y_Energy,log10(Data(k).eflux(t1:1:t2,:))'); 
        colormap('jet'); 
        xlim([DateNumBeg, DateNumEnd]);
        caxis([2 8]);
        numberofYTicks=5;
        yAxisVals=linspace(min(Data(k).Y_Energy), max(Data(k).Y_Energy), numberofYTicks);
        shading flat;
        set(ax(j),'yscale','log','YTick',yAxisVals,'YTickLabel',num2str(yAxisVals,'%10.1e\n'));
        if(j<n(1))
        set(ax(j),'XTick',[]);
        else
        numberofXTicks=ceil((DateNumEnd-DateNumBeg)*24)*2-1;
        xAxisVals = linspace (DateNumBeg, DateNumEnd, numberofXTicks);
        set(ax(j),'XTick',xAxisVals,'XTickLabel',datestr(xAxisVals,'HH:MM')); 
        end;
        %colorbar('manual','position',[0.93 0.55 0.02 0.35]);
        colorbar();
        linkaxes(ax,'x');
        ax(j).Position=[0.1 0.9-(j*0.8/n(1)) 0.8 0.75/n(1)];
        ylabel(ax(j),[eflux(j,:),' ','energy [eV]']);
        k=k+1;
    end;
        xlabel(ax(j),'Time [UT]');
        title(ax(1),['Probe: ',Probes(i,:),' Particle Energy Spectrum Date: ',datestr(DateNumBeg,'dd-mm-yyyy HH:MM'),' to ',datestr(DateNumEnd,'dd-mm-yyyy HH:MM')]);
end

