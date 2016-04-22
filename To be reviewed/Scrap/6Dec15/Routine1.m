%% Plotting the ratio of ion and electron densities

load('/home/nithin/Documents/git-repos/aurora-rep/To be reviewed/Scrap/2Dec15/thd.mat');
load('/home/nithin/Documents/git-repos/aurora-rep/To be reviewed/Scrap/6Dec15/thddensity.mat');
% data=downloadSPDF('26 Mar 2008 00:01','26 Mar 2008 23:59','dd mm yyyy HH:MM','themis',char('thd'),'l2','esa','/home/nithin/Documents/git-repos/aurora-rep/To be reviewed/Scrap/6Dec15/');

thddensityp=interp1(thddensity.ptime,thddensity.p,thddensity.etime);
figure;
scatter(thddensityp(:,:),thddensity.e(:,:),'.');
hold on; plot((0.01:0.01:100),0.94*(0.01:0.01:100).^0.86,'b');
for i=1:1:length(thdesa.time)
    thd_high_energy_e_density(i)=sum(thdesa.eflux(i,21:end).*(thdesa.ebin(21:end)-thdesa.ebin(20:end-1))./(thdesa.ebin(21:end)),2);
end;

% Marking the Plasma Sheet Boundary Layer (PSBL)
mark_PSBL=thddensity.e>((thddensityp));
mark_taillobe=thd_high_energy_e_density<1*10^-4;


%% Plotting
DateNumBeg=datenum('26-Mar-2008 08:00');
DateNumEnd=datenum('26-Mar-2008 13:00');
figure;
subplot(2,1,1)
plot(thddensity.etime,thddensity.e);
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','Density [cm^-3]');
title('Categorization of regions within the plasmasheet');
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);
hold on;
plot(thddensity.etime,thddensityp);
hold on;
plot(thddensity.etime,mark_PSBL);

subplot(2,1,2)
plot(thdscp.time,thdscp.scp);
set(gca,'XTick',(DateNumBeg:(0.5/24):DateNumEnd));
set(gca,'XTickLabel',{datestr(DateNumBeg:(0.5/24):DateNumEnd,'HH:MM')});
set(get(gca,'XLabel'),'String','Universal Time [HH:MM]');
set(get(gca,'YLabel'),'String','SC potential [V]');
title('Spacecraft Potential');
set(gca,'YScale','log');
xlim([DateNumBeg DateNumEnd]);

