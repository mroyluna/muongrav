%analyse.m -- 12/1/16 -- by MR
clear all
close all

load Trend_Lower.dat
load Trend_Upper.dat
TU = Trend_Upper(:,3);% - mean(Trend_Upper(:,3));
TL = Trend_Lower(:,3);% - mean(Trend_Lower(:,3));

cd results_12_2_16/

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
    le  = str2num(nstr(ind_s(4)+1:ind_s(5)-1));
    slopeN = str2num(nstr(ind_s(5)+1:ind_s(6)-1));
    slopeE = str2num(nstr(ind_s(6)+1:inddot-1));

    east = dat(:,1);
    nort = dat(:,2);
    elev = dat(:,3);
    calc = dat(:,4);
    gobs = dat(:,5);
    gerr = dat(:,6);
    
    
    below_cutoff_height = elev < 2150;
    UpObs   = gobs(~below_cutoff_height);
    nUp     = numel(UpObs);
    LowObs  = gobs(below_cutoff_height);
    nLow    = numel(LowObs);
    UpCalc  = calc(~below_cutoff_height);
    LowCalc = calc(below_cutoff_height);
    
    UpObs2  = UpObs-TU+mean(UpObs);
    LowObs2 = LowObs-TL+mean(LowObs);
    rmsUp   = sqrt(sum((UpObs2 - UpCalc).^2)/nUp);
    rmsLow  = sqrt(sum((LowObs2 - LowCalc).^2)/nLow);
    allObs  = [UpObs2;LowObs2];
    allCalc = [UpCalc;LowCalc];
    rmsall  = sqrt(sum((allObs - allCalc).^2)/(nLow+nUp)); 
    
    if rd == 2300 & drd == -900 & le == 2160
        rmsdata(count,:) = [n, rd, drd, le, rmsUp, rmsLow, rmsall, slopeN, slopeE];
        count = count + 1;
    end
    
end

cd ../

%layerelev = [2150:10:2170];
rhospaceU=[1300:50:2700];
rhospaceL=rhospaceU;
xr = rhospaceL;
yr = rhospaceU';
%rhospaceL=[min(LowDensity):50:max(LowDensity)];

slopespace=[-10:0.25:10];
ys = slopespace';
xs = [-3:0.25:5];

layerelev = 2160;

for i = 1:numel(layerelev)
    mask = find(rmsdata(:,4)==layerelev(i));
    lstr = num2str(layerelev(i));
    UpDensity    = rmsdata(mask,2);
    DeltaDensity = rmsdata(mask,3);
    rmsU         = rmsdata(mask,5);
    rmsL         = rmsdata(mask,6);
    rmsA         = rmsdata(mask,7);
    sly          = rmsdata(mask,8); % array of y-slopes
    slx          = rmsdata(mask,9); % array of x-slopes
    
    LowDensity = DeltaDensity + UpDensity;
    npts = numel(UpDensity);
    % %make array for contouring
    
    rmsUint=griddata(slx,sly,rmsU,xs,ys);
    rmsLint=griddata(slx,sly,rmsL,xs,ys);
    rmsAint=griddata(slx,sly,rmsA,xs,ys);
    
    
%     plot3(LowDensity,UpDensity,rmsU,'o');hold on
%     mesh(xr,yr,rmsUint);
    figure(); orient tall
    C3 = contourf(xs,ys,rmsAint,'w-','linewidth',[2]);hold on
    clabel(C3,'FontSize',13,'Color','w','FontWeight','bold')
    C1 = contour(xs,ys,rmsUint,'k--','linewidth',[1]);hold on
    clabel(C1,'FontSize',13,'Color','k')
    C2 = contour(xs,ys,rmsLint,'r--','linewidth',[1]);
    clabel(C2,'FontSize',13,'Color','r')
    title(['LayerElev=' lstr ])
    set(gca,'FontName','Helvetica','FontSize',[16])
    axis equal
    colormap('jet');caxis([0 2]); colorbar
    xlabel('E-slope, degrees')
    ylabel('N-slope, degrees')
    set(gca,'FontName','Helvetica','FontSize',[16])
    filename = ['LayerElev_tilt',lstr,'_RMS.png'];
    eval(['print ', filename, ' -dpng'])
end    

%legend({'Upper stations','Lower stations','All stations'},'FontSize',16,'Location','southeastoutside')
