%analyse.m -- 12/1/16 -- by MR
clear all
close all

cd results_12_2_16/

files=dir('AllData_400_threelayer*.dat') 
%file format is 
% dataarray=[easting',northing',elevations',offset_gz_vals,gz_avg_at_stations,gz_error_at_stations];
%filename=['AllData_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.dat'];

for i=1:length(files)
    nstr = files(i).name;
    dat = load(nstr);
    
    ind_s   = find(nstr == '_'); % find position of underscores in name
    ind_dot = find(nstr == '.'); % find position of dot in name

    n   = str2num(nstr(ind_s(1)+1:ind_s(2)-1));
    rd  = str2num(nstr(ind_s(2)+1:ind_s(3)-1));
    drd = str2num(nstr(ind_s(3)+1:ind_s(4)-1));
    le  = str2num(nstr(ind_s(4)+1:ind_dot-1));

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
    rmsUp   = sqrt(sum((UpObs - UpCalc).^2)/nUp);
    rmsLow  = sqrt(sum((LowObs - LowCalc).^2)/nLow);
    
    rmsall  = sqrt(sum((gobs - calc).^2)/(nLow+nUp));
    rmsdata(i,:) = [n, rd, drd, le, rmsUp, rmsLow, rmsall];
 
end

cd ../

layerelev = [2150:10:2190];
%layerelev = 2130;
figure(1);
for i = 1:numel(layerelev)
    mask = find(rmsdata(:,4)==layerelev(i));
    lstr = num2str(layerelev(i));
    UpDensity    = rmsdata(mask,2);
    DeltaDensity = rmsdata(mask,3);
    rmsU         = rmsdata(mask,5);
    rmsL         = rmsdata(mask,6);
    rmsA         = rmsdata(mask,7);
    LowDensity = DeltaDensity + UpDensity;
    rhospaceU=[min(UpDensity):50:max(UpDensity)];
    rhospaceL=[min(LowDensity):50:max(LowDensity)];
    npts = numel(UpDensity);
    % %make array for contouring
    xr = rhospaceL;
    yr = rhospaceU';
    rmsUint=griddata(LowDensity,UpDensity,rmsU,xr,yr);
    rmsLint=griddata(LowDensity,UpDensity,rmsL,xr,yr);
    rmsAint=griddata(LowDensity,UpDensity,rmsA,xr,yr);
    
    
%     plot3(LowDensity,UpDensity,rmsU,'o');hold on
%     mesh(xr,yr,rmsUint);
    subplot(2,3,i)
    C1 = contour(xr,yr,rmsUint,'k-','linewidth',[1]);hold on
    clabel(C1,'FontSize',13,'Color','k')
    C2 = contour(xr,yr,rmsLint,'r-','linewidth',[1]);
    clabel(C2,'FontSize',13,'Color','r')
    C3 = contour(xr,yr,rmsAint,'b-','linewidth',[2]);
    clabel(C3,'FontSize',13,'Color','b','FontWeight','bold')
    title(['LayerElev=' lstr ])
    set(gca,'FontName','Helvetica','FontSize',[16])
    axis equal
end    

legend({'Upper stations','Lower stations','All stations'},'FontSize',16,'Location','southeastoutside')
xlabel('Lower Layer Density')
ylabel('Upper Layer Density')
set(gca,'FontName','Helvetica','FontSize',[16])