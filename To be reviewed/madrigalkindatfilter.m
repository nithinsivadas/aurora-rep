%% Program that can Filter experiments with Baker Code : Kindat 5963
%  Which possibly has low altitude data < 100 km

clear;

% Defining the URL of the Madrigal Database
cgiurl='http://isr.sri.com/madrigal/cgi-bin/';

% Collecting the experiment arrays that satisfy the following data and time
% criterion
expArray = getExperimentsWeb(cgiurl, 61, datenum('01/01/2008 00:00:00'),...
    datenum('01/01/2012 00:00:00'), 1);

%expArray = getExperimentsWeb(cgiurl, 61, datenum('03/26/2008 10:00:00'),...
%    datenum('03/26/2008 13:00:00'), 1);

[d1,d2]=size(expArray);

% Collecting the file names
k=0;
for i=1:1:d2
expFileArray = getExperimentFilesWeb(cgiurl,expArray(i).id);
[d3,d4]=size(expFileArray);
    for j=1:1:d4
        if expFileArray(j).kindat==5963
            k=k+1;
            expFileArrayStore(k)=expFileArray(j);
        end;
    end;
end;



% Listing down the parameters in the files\
%parmArray = getParametersWeb(cgiurl, expFileArray(1).name);
