function [ ExpFileArray ] = downloadSPDF( DateBeg, DateEnd, DateFormat,...
    Mission, Probes, DataLevel, Instrument, TargetPath)

%% downloadSPDF: Downloading CDF Files for a particular mission
%   The function first produces the data file-path and file-names that 
%   contain THEMIS data within the user input data ranges. It then connects
%   to the FTP server and downloads those particular data files in CDF format.

%   Input
%   DateBeg: String that specifies the beginning date of collecting data (UT)
%   DateEnd: String specifying end date
%   DateFormat: String specifying the format in which Dates are written
%   Mission: String specifyin the name of the mission within the
%   ftp://spdf.gsfc.nasa.gov/
%   Probes: String specifying the probe within the mission
%   Date Level: String specifying the level of data
%   Instrument: String specifying the Instrument name

%   Ouput
%   ExpFileArray: Identifies and creates a file name and path array based
%   on the values of the input variables

%   v01|5th Aug 2015: Currently optimized for THEMIS satellites. Note that 
%   this code can be used to download other mission files, as long as you 
%   can estimate the file names and paths of the data files. 

%% Example of Function Call
% X=downloadSPDF('30 Mar 2012 09:00','31 Mar 2012 10:00','dd mm yyyy HH:MM','themis',char('tha','thb'),'l2','sst','/home/nithin/NithinBU/Aug 4 2015/');

% Converting the user input dates into matlab date format
DateNumBeg=datenum(DateBeg,DateFormat);
DateNumEnd=datenum(DateEnd,DateFormat);

% Error message
if(DateNumBeg>=DateNumEnd)
     error('The beginning date is later than or equal to the end date');
end;

% Calculating the file path and name within the ftp server
k=1;       
for j=1:1:size(Probes)
    for i=1:1:ceil(DateNumEnd)-floor(DateNumBeg)
    filePath(k,:)=['/pub/data/',Mission,'/',Probes(j,:),'/',DataLevel,...
        '/',Instrument,'/',datestr(floor(DateNumBeg)+i-1,'yyyy')];
    fileName(k,:)=[Probes(j,:),'_',DataLevel,'_',Instrument,'_',...
        datestr(floor(DateNumBeg)+i-1,'yyyymmdd'),'_v01.cdf'];
    k=k+1;
    end;
end;

% Connecting to the ftp server, and downloading the data files 
% into the target folder 
filePath
spdf_ftp_client=ftp('spdf.gsfc.nasa.gov');
for i=1:1:size(filePath)
    cd(spdf_ftp_client,filePath(i,:));
    mget(spdf_ftp_client,fileName(i,:),TargetPath);
    ExpFileArray(i,:)=[TargetPath,fileName(i,:)];
end;
close(spdf_ftp_client);

end

