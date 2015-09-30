%% Processing AE data into data structures 
% Both AU and AL entries are sampled at 1 second intervals and 
% stored in two arrays with matlab time for each year. The output is .mat file with 
% a structure S.ae_mod_2013.au and S.ae_mod_2013.al etc.

%% Issues
    % Think of whether we can access AE data more directly, or if we can store them 
    % somewhere in this processed format with matlab time so that it can be used 
    % easily. 

clearvars;

k=1;
m=1;
n=1;
l=1;
q=1;

% The years of AE data available
date=[2012,2013,2014];

% The AE data is in the folder in the following format: ae_2012_minute.mat
% Loading the .mat file into MATLAB memory
for q=1:1:size(date')
    filename=strcat('ae_',num2str(date(q)),'_minute.mat');
    array=struct('file',load(filename));
    aeminute(q)=array;
end;

% Looping the date variable
for p=1:1:size(date')
    fields=fieldnames(aeminute(p).file);
    ae20xxminute=aeminute(p).file.(fields{1});
    k=1;
    m=1;
    n=1;
    l=1;
    
    % Looping the entries within the AE data
    for i=1:1:size(ae20xxminute(:,1))
        % Creating an array of AU entries for each second of the year vs. matlab time
        if (ae20xxminute{i,5}=='AU')
            au(k).YMD=[rem(floor(ae20xxminute{i,2}/10000),100)+2000, rem(floor(ae20xxminute{i,2}/100),100), rem(floor(ae20xxminute{i,2}),100)];
            au(k).HH=ae20xxminute{i,4};
            au(k).mm=[0:1:59];
            au(k).data=ae20xxminute(i,7:66);
                for j=1:1:60
                au_arr(l,1)=[datenum([au(k).YMD(1),au(k).YMD(2),au(k).YMD(3),au(k).HH,au(k).mm(j),0])];
                au_arr(l,2)=au(k).data{j};
                l=l+1;
                end;
         k=k+1;
        end;
    
        if (ae20xxminute{i,5}=='AL')
        % Creating an array of AU entries for each second of the year vs. matlab time
            al(m).YMD=[rem(floor(ae20xxminute{i,2}/10000),100)+2000, rem(floor(ae20xxminute{i,2}/100),100), rem(floor(ae20xxminute{i,2}),100)];
            al(m).HH=ae20xxminute{i,4};
            al(m).mm=[0:1:59];
            al(m).data=ae20xxminute(i,7:66);     
                for j=1:1:60
                al_arr(n,1)=[datenum([al(m).YMD(1),al(m).YMD(2),al(m).YMD(3),al(m).HH,al(m).mm(j),0])];
                al_arr(n,2)=al(m).data{j};
                n=n+1;
                end;
            m=m+1;
        end;
       
    end;
    varname=strcat('ae_mod_',num2str(date(p)))
    %p
    % Creating a structure containing the two arrays, and saving them as .mat file
    ae_mod_20xx.au=au_arr;
    ae_mod_20xx.al=al_arr;
    S.(varname)=ae_mod_20xx; 
    save(strcat(varname,'.mat'),'-struct','S');
end;

 
%count=1; 
%for i=1:1:size(ae20xxminute(:,1)) 
%if(ae20xxminute{i,5}=='AU') 
%count=count+1; 
%end; 
%end;
