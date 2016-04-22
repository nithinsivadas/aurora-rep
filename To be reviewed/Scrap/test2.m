clear;

[ne,altitude,t]=pfisr_mag_e('DataFile_2008_1.h5',1);
% ne is in m^-3 so converting to cm^-3
ne=ne*10^-6;
% [ne1,altitude1,t1]=pfisr_avg_e('DataFile_2008_1.h5');
% min_ne=min(min(ne));
% Removing NANs
% [ii jj] = find(isnan(ne));
% [ii0 jj0]=find((isnan(ne))<1);
% ne(ii,jj) = (ne(ii-1,jj)+ne(ii+1,jj))/2;

x=altitude;
for t1=1:1:size(ne,2)
    y=ne(:,t1);
    xi=x(find(~isnan(y)));
    yi=y(find(~isnan(y)));
    ne(:,t1)=interp1(xi,yi,x,'pchip','extrap');
end;
% 
% x=t;
% for a1=1:1:size(ne,1)
%     y=ne(a1,:);
%     xi=x(find(~isnan(y)));
%     yi=y(find(~isnan(y)));
%     ne(a1,:)=interp1(xi,yi,x,'spline','extrap');
% end;


% ne(ii,jj) = (ne(ii-1,jj)+ne(ii+1,jj))/2;
% for i=1:1:size(ii)
%     ne(ii(i),jj(i))=min_ne;
% end;

 

[T,ALT]=meshgrid(t,altitude);	
figure; surf(T,ALT,real(log10(ne(:,:))),'EdgeColor','none');
set(get(gca,'XLabel'),'String','Time');datetick('x', 'HH:MM');
set(get(gca,'YLabel'),'String','Altitude [km]');
title('Electron Density measured from PFISR');
view([0 90]);
set(gca,'YScale','log');
xlim([min(min(T)) max(max(T))]);
ylim([min(min(ALT)) max(max(ALT))]);