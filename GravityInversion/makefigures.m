%analyse.m -- 12/1/16 -- by MR
clear all
close all

load Trend_Lower.dat
load Trend_Upper.dat
TU = Trend_Upper(:,3);% - mean(Trend_Upper(:,3));
TL = Trend_Lower(:,3);% - mean(Trend_Lower(:,3));
eTU = Trend_Upper(:,4);% - mean(Trend_Upper(:,3));
eTL = Trend_Lower(:,4);% - mean(Trend_Lower(:,3));

[point_table, measured_points] = build_table();
measure_errors = point_table{measured_points, 'Errors'};
gz_error = cellfun(@norm, measure_errors);

cd results_12_1_16/
%low uniformdensity models are in 12_4; high uniform density and two-layer in 12_1
%two_layer_tilted in 12_4

files=dir('AllData_*.dat') 
%file format is 
% dataarray=[easting',northing',elevations',offset_gz_vals,gz_avg_at_stations,gz_error_at_stations];
%filename=['AllData_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.dat'];
count=1;
for i=1:length(files)
    nstr = files(i).name;
    dat = load(nstr);
    
    ind_s   = find(nstr == '_'); % find position of underscores in name
    ind_dot = find(nstr == '.'); % find position of dot in name
    inddot  = ind_dot(end); % there may be more than one dot
    n   = str2num(nstr(ind_s(1)+1:ind_s(2)-1));
    rd  = str2num(nstr(ind_s(2)+1:ind_s(3)-1));
    drd = str2num(nstr(ind_s(3)+1:ind_s(4)-1));
    le  = str2num(nstr(ind_s(4)+1:inddot-1));
    %slopeN = str2num(nstr(ind_s(5)+1:ind_s(6)-1));
    %slopeE = str2num(nstr(ind_s(6)+1:inddot-1));

    if  rd == 2300 & drd == -900 & le == 2160 
        east = dat(:,1);
        nort = dat(:,2);
        elev = dat(:,3);
        calc = dat(:,4);
        gobs = dat(:,5);
        gerr = dat(:,6);
        
        refgz = gobs(strcmp(measured_points, 'BS_TN_1'));

        below_cutoff_height = elev < 2150;
        UpObs   = gobs(~below_cutoff_height);
        %UpObsErr= gerr(~below_cutoff_height);
        UpObsErr= gz_error(~below_cutoff_height);
        nUp     = numel(UpObs);
        Upnort  = nort(~below_cutoff_height);
        Upeast  = east(~below_cutoff_height);
        
        LowObs  = gobs(below_cutoff_height);
        nLow    = numel(LowObs);
        Lonort  = nort(below_cutoff_height);
        Loeast  = east(below_cutoff_height);
        LoObsErr= gz_error(below_cutoff_height);
        
        UpCalc  = calc(~below_cutoff_height);
        LowCalc = calc(below_cutoff_height);
        rmsUp   = sqrt(sum((UpObs - UpCalc).^2)/nUp);
        rmsLow  = sqrt(sum((LowObs - LowCalc).^2)/nLow);
        
        rmsall  = sqrt(sum((gobs - calc).^2)/(nLow+nUp));
        
        errU  = sqrt(eTU.*eTU + UpObsErr.*UpObsErr);
        errL  = sqrt(eTL.*eTL + LoObsErr.*LoObsErr);
        
            
        refind = find(LowObs == refgz);
        U = UpObs-refgz-TU+TL(refind);% make everything relative to ref station
        L = LowObs-refgz-TL+TL(refind);%
        
        
        
    end
   
end

        figure(1); hold on;
            errorbar(Upnort, U, errU, 'rd','markerfacecolor','r');
            scatter(Upnort, UpCalc,'ko','markerfacecolor','k')
            %errorbar(Upnort, UpObs, UpObsErr*10, 'bo','markerfacecolor','c');
            %errorbar(Upnort, TU, UpObsErr*10, 'ro','markerfacecolor','g');
            title(['Upper n = ' num2str(n) ', top density = ' num2str(rd) ', delta density = ' num2str(drd) ', layer elevation = ' num2str(le)]);
            %legend('Observed','Regional trend','Residual','Location','northwest');
            %legend('Observed - trend','1200 kg/m^3','1300 kg/m^3','1400 kg/m^3','Location','northwest');
            legend('Observed - trend','1700 kg/m^3','1900 kg/m^3','2100 kg/m^3','Location','northwest');
            xlabel('Northing (m)'); ylabel('gz (mgal)');
            set(gca,'fontsize',[14])
            filename=['../N_Upper_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
            filename=['../N_Upper_observations.png'];
            %eval(['print ' filename ' -dpng']);
        figure(2); hold on;
            errorbar(Upeast, U, errU, 'rd','markerfacecolor','r');
            scatter(Upeast, UpCalc,'ko','markerfacecolor','k')
            %errorbar(Upeast, UpObs, UpObsErr*10, 'bo','markerfacecolor','c');
            %errorbar(Upeast, TU, UpObsErr*10, 'ro','markerfacecolor','g');
            title(['Upper n = ' num2str(n) ', top density = ' num2str(rd) ', delta density = ' num2str(drd) ', layer elevation = ' num2str(le)]);
            %legend('Calculated',  'Observed'); xlabel('Easting (m)'); ylabel('gz (mgal)');
            %legend('Observed - trend','1200 kg/m^3','1300 kg/m^3','1400 kg/m^3')
            legend('Observed - trend','1700 kg/m^3','1900 kg/m^3','2100 kg/m^3')
            %legend('Observed','Regional trend','Residual'); 
            %add tunnel location at x=4.95755e5;
            plot([4.95755e5,4.95755e5],[min(U),max(U)],'k--');
            xlabel('Easting (m)'); ylabel('gz (mgal)');
            set(gca,'fontsize',[14])
            filename=['../E_Upper_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
            filename=['../E_Upper_observations.png'];
            %eval(['print ' filename ' -dpng']);
        figure(3); hold on;
            errorbar(Lonort, L, errL, 'rd','markerfacecolor','r');
            scatter(Lonort, LowCalc,'ko','markerfacecolor','k')
            %errorbar(Lonort, LowObs, LoObsErr*2, 'bo','markerfacecolor','c');
            %errorbar(Lonort, TL, LoObsErr*2, 'ro','markerfacecolor','g');
            title(['Upper n = ' num2str(n) ', top density = ' num2str(rd) ', delta density = ' num2str(drd) ', layer elevation = ' num2str(le)]);
            %legend('Calculated',  'Observed');xlabel('Northing (m)'); ylabel('gz (mgal)');
            %legend('Observed','Regional trend','Residual'); 
            %legend('Observed - trend','Calculated');
            %legend('Observed - trend','1200 kg/m^3','1300 kg/m^3','1400 kg/m^3')
            legend('Observed - trend','1700 kg/m^3','1900 kg/m^3','2100 kg/m^3')
            xlabel('Northing (m)'); ylabel('gz (mgal)');
            set(gca,'fontsize',[14])
            filename=['../N_Lower_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
            filename=['../N_lower_observations.png'];
            %eval(['print ' filename ' -dpng']);
        figure(4); hold on;
            errorbar(Loeast, L, errL, 'rd','markerfacecolor','r');
            scatter(Loeast, LowCalc,'ko','markerfacecolor','k')
            %errorbar(Loeast, LowObs, LoObsErr*2, 'bo','markerfacecolor','g');
            %errorbar(Loeast, TL, LoObsErr*2, 'ro','markerfacecolor','g');
            title(['Upper n = ' num2str(n) ', top density = ' num2str(rd) ', delta density = ' num2str(drd) ', layer elevation = ' num2str(le)]);
            %legend('Calculated',  'Observed');xlabel('Easting (m)'); ylabel('gz (mgal)');
            %legend('Observed','Regional trend','Residual'); 
            %legend('Observed - trend','Calculated');
            %legend('Observed - trend','1200 kg/m^3','1300 kg/m^3','1400 kg/m^3')
            legend('Observed - trend','1700 kg/m^3','1900 kg/m^3','2100 kg/m^3')
            xlabel('Easting (m)'); ylabel('gz (mgal)');
            set(gca,'fontsize',[14])
            filename=['../E_lower_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
            filename=['../E_lower_observations.png'];
            %eval(['print ' filename ' -dpng']);

%         figure(1);
%             filename=['../N_Upper_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
%             eval(['print ' filename ' -dpng']);
%         figure(2); hold on
%             filename=['../E_Upper_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
%             eval(['print ' filename ' -dpng']);
%         figure(3); hold on;
%             filename=['../N_Lower_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
%             eval(['print ' filename ' -dpng']);
%         figure(4); hold on;
%             filename=['../E_lower_', num2str(n),'_',num2str(rd),'_',num2str(drd),'_',num2str(le),'.png'];
%             eval(['print ' filename ' -dpng']);
%             

cd ../
