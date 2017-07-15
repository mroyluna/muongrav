%this code will take the data files from the three layer runs with tilts
%and find the lowest rms
clear all

TL = load('Trend_Lower.dat');
TU = load('Trend_Upper.dat');

%%
cd Out_2130_1900_350_-1000

files=dir('AllDataSepTilt_*.dat') ;
disp(['number of files = ' num2str(length(files))])

numstations=60;
rmslist=zeros(1,length(files));

    for i=1:length(files)
        nAmes=files(i).name;
        naMes=strrep(nAmes, '.dat','');
        namEs=strrep(naMes, '-','_');
        names=strrep(namEs, '.','_');
        eval(['load ' files(i).name]);
        mats=eval([names]);
        gobs = mats(:,5);
        easting = mats(:,1);
        northing = mats(:,2);
        calcg    = mats(:,4);
        elev     = mats(:,3);
        
below_cutoff_height = elev < 2150;
UpObs   = gobs(~below_cutoff_height);
UpCalc  = calcg(~below_cutoff_height);
nUp     = numel(UpObs);
LowObs  = gobs(below_cutoff_height);
LowCalc = calcg(below_cutoff_height);
nLow    = numel(LowObs);
xU      = easting(~below_cutoff_height);
yU      = northing(~below_cutoff_height);
xL      = easting(below_cutoff_height);
yL      = northing(below_cutoff_height);

%now remove trend

UpObs = UpObs - TU(:,3) + mean(TU(:,3));
LowObs = LowObs - TL(:,3) + mean(TL(:,3));

% figure(1);clf
% plot(xU, UpObs, 'd');hold on
% plot(xU, UpCalc, '*');
% pause(0.1)
       rms=sqrt((sum((UpCalc - UpObs).^2)+ sum((LowCalc - LowObs).^2))/numstations);
      %gtrend1 comes from running bouger_trend first
       rmslist(i)=rms;
       
       
    end
cd ../

   [BestRMS, position]=min(rmslist)
   BestFit=files(position)
   
   %then i will plot obs and calc values of the best fit file

   