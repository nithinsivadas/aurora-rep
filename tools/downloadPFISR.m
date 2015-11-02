function [ Data, DateNumBeg, DateNumEnd ] = downloadPFISR( DateBeg, DateEnd, DateFormat )
%downloadPFISR : Downloads PFISR hdf5 data file 
%	The function first searches out experiment data file-path and file-names that 
%   contain PFISR experiments within the user input data ranges. And then downloads the HDF5 files.

%   Input
%   DateBeg: String that specifies the beginning date of collecting data (UT)
%   DateEnd: String specifying end date
%   DateFormat: String specifying the format in which Dates are written

%	Output
%	Data: FileNameStr => A string of the downloaded HDF5 file names. 

%% Issues:
    % Currently the code selects only the Barker code data (kindat==5963)
    % Currently it only downloads PFISR data (Instrument Code = 61)
    
% Defining the URL of the Madrigal Database
cgiurl='http://isr.sri.com/madrigal/cgi-bin/';

% Converting the user input dates into matlab date format
DateNumBeg=datenum(DateBeg,DateFormat);
DateNumEnd=datenum(DateEnd,DateFormat);

% Error message
if(DateNumBeg>=DateNumEnd)
     error('The beginning date is later than or equal to the end date');
end;


% Collecting the experiment arrays that satisfy the following data and time  criterion
expArray = getExperimentsWeb(cgiurl, 61, DateNumBeg,...
    DateNumEnd, 1); % Here 61 stands for PFISR instrument

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

% Downloading the HDF5 files of all the experiments on to the active directory.
for i=1:1:k
    fileNameStr=sprintf('DataFile_%s_%d.h5',expFileArrayStore(i).name(27:30),i);
    result=madDownloadFile(cgiurl, expFileArrayStore(i).name,fileNameStr,'Nithin Sivadas','nithin@bu.edu','Boston University','hdf5');
    Data(k)=cellstr(fileNameStr);
end;

end

