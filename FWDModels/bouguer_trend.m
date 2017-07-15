% find trend in regional bouguer from PACES database
% Mousumi Roy -- Dec 4 2016
clear all
close all
%%
load bouguer_16_clean.xyg;
dat = bouguer_16_clean;
udat = unique(dat,'rows');

% bounds for larger region from which to get data
% -106 deg 21; -106 deg 12; 35 deg 51; 35 deg 57;

regionx = [-106-21/60 : 0.001 : -106-12.5/60];
regiony = [35+51/60   : 0.001 : 35+57/60];


inregion = udat(:,1) >= regionx(1) & udat(:,1) <= regionx(end) & udat(:,2) >= regiony(1) & udat(:,2) <= regiony(end); 
indat = udat(inregion,:);

lon = indat(:,1);
lat = indat(:,2);
gboug = indat(:,3);

[xq, yq, intboug]=griddata(lon, lat, gboug, regionx, regiony');
%study area is in this box:
%  -106.29975 35.87566
% -106.29455 35.87566
% -106.29455 35.87997
% -106.29975 35.87997
% -106.29975 35.87566
%%
box = [-106.28 35.85; -106.33 35.85; -106.33 35.9; -106.265 35.9;-106.265 35.85];

%used the website to get lat and lon of region of lidar data: 
%http://www.ngs.noaa.gov/cgi-bin/spc_getgp.prl
ll   = [35+52/60+6.76930/3600   -106-18/60-16.93468/3600];
ur   = [35+53/60+8.76740/3600   -106-17/60-21.15378/3600];

lidarbox = [ll(2) ll(1);ll(2) ur(1); ur(2) ur(1); ur(2) ll(1);ll(2) ll(1)];

%from lidar data (we know):
minX = 4.9506e+05;
maxX = 4.9646e+05;
minY = 5.3993e+05;
maxY = 5.4184e+05;

%convert interpolated bouguer inside box into northing and easting coords
xcoef = polyfit([ll(2);ur(2)],[minX;maxX],1);
ycoef = polyfit([ll(1);ur(1)],[minY;maxY],1);

east = xcoef(1)*xq + xcoef(2);
nort = ycoef(1)*yq + ycoef(2);
peast = xcoef(1)*lon + xcoef(2);
pnort = ycoef(1)*lat + ycoef(2);

figure(1)
contourf(xq, yq,intboug);hold on
colorbar
plot(lon, lat, 'o','MarkerFaceColor','w')
plot(box(:,1),box(:,2),'w-','linewidth',[2]);
plot(lidarbox(:,1),lidarbox(:,2),'r-','linewidth',[2]);
xlabel('Longitude'); ylabel('Latitude'); set(gca,'FontSize',[14])
%print regional_bouguer_trend_lat_lon_largebox.png -dpng

box = [box(:,1)*xcoef(1) + xcoef(2) box(:,2)*ycoef(1) + ycoef(2)];
lidarbox = [lidarbox(:,1)*xcoef(1) + xcoef(2) lidarbox(:,2)*ycoef(1) + ycoef(2)];
%box = lidarbox;
 
%detrend in box
intdata = [peast, pnort, gboug];
inbox = intdata(:,1) >= min(box(:,1)) & intdata(:,1) <= max(box(:,1)) & intdata(:,2) >= min(box(:,2)) & intdata(:,2) <= max(box(:,2)); 
indat = intdata(inbox,:);

[n_1,V,p_1] = affine_fit(indat);

figure(2)
contourf(east, nort,intboug);hold on
colorbar
plot(peast, pnort, 'o','MarkerFaceColor','w')
plot(box(:,1),box(:,2),'w-','linewidth',[2]);
plot(lidarbox(:,1),lidarbox(:,2),'r-','linewidth',[2]);
plot(indat(:,1),indat(:,2),'ro');
xlabel('Easting, m'); ylabel('Northing, m'); set(gca,'FontSize',[14])
print regional_bouguer_trend_largebox.png -dpng

xspace = [min(box(:,1)):50:max(box(:,1))];
yspace = [min(box(:,2)):50:max(box(:,2))];
[X, Y] = meshgrid(xspace, yspace);

plane1 = - (n_1(1)/n_1(3)*X+n_1(2)/n_1(3)*Y-dot(n_1,p_1)/n_1(3));
figure(3)
scatter3(indat(:,1),indat(:,2),indat(:,3),'bo');hold on
surf(X,Y, plane1,'facecolor','red','facealpha',0.15,'edgealpha',0.15);
plot3(p_1(1),p_1(2),p_1(3),'ro','markersize',15,'markerfacecolor','red');
xlabel('Easting, m'); ylabel('Northing, m'); 
zlabel('Bouguer Anomaly, mgals')
set(gca,'FontSize',[14])
view([16.9, 36.4])
%print fit_linear_trend_largebox.png -dpng

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elev     = eval_pts(3, :);
northing = eval_pts(2, :);
easting  = eval_pts(1,:);

measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

gz_avg_at_stations = cellfun(@mean, measured_values);
gz_error_at_stations = cellfun(@norm, measure_errors);

gobs = gz_avg_at_stations';
gerr = gz_error_at_stations';

below_cutoff_height = elev < 2150;
UpObs   = gobs(~below_cutoff_height);
nUp     = numel(UpObs);
LowObs  = gobs(below_cutoff_height);
nLow    = numel(LowObs);
xU      = easting(~below_cutoff_height);
yU      = northing(~below_cutoff_height);
xL      = easting(below_cutoff_height);
yL      = northing(below_cutoff_height);

gshift = mean(plane1(:)) - mean(UpObs);
plane2 = plane1 - gshift;
gshift = mean(plane1(:)) - mean(LowObs);
plane3 = plane1 - gshift;

F2 = scatteredInterpolant(X(:),Y(:),plane2(:));
F3 = scatteredInterpolant(X(:),Y(:),plane3(:));
gtrendU1 = F2(xU,yU);
gtrendL1 = F3(xL,yL);


%%
%repeat with lidar box sized region

lidarbox = [ll(2) ll(1);ll(2) ur(1); ur(2) ur(1); ur(2) ll(1);ll(2) ll(1)];
box = lidarbox;

box = [box(:,1)*xcoef(1) + xcoef(2) box(:,2)*ycoef(1) + ycoef(2)];
lidarbox = [lidarbox(:,1)*xcoef(1) + xcoef(2) lidarbox(:,2)*ycoef(1) + ycoef(2)];
%box = lidarbox;
 
%detrend in box
intdata = [peast, pnort, gboug];
inbox = intdata(:,1) >= min(box(:,1)) & intdata(:,1) <= max(box(:,1)) & intdata(:,2) >= min(box(:,2)) & intdata(:,2) <= max(box(:,2)); 
indat = intdata(inbox,:);

[n_1,V,p_1] = affine_fit(indat);

xspace = [min(box(:,1)):50:max(box(:,1))];
yspace = [min(box(:,2)):50:max(box(:,2))];
[X, Y] = meshgrid(xspace, yspace);

plane1 = - (n_1(1)/n_1(3)*X+n_1(2)/n_1(3)*Y-dot(n_1,p_1)/n_1(3));
figure(3)
scatter3(indat(:,1),indat(:,2),indat(:,3),'bo');hold on
surf(X,Y, plane1,'facecolor','red','facealpha',0.15,'edgealpha',0.15);
plot3(p_1(1),p_1(2),p_1(3),'ro','markersize',15,'markerfacecolor','red');
xlabel('Easting, m'); ylabel('Northing, m'); 
zlabel('Bouguer Anomaly, mgals')
set(gca,'FontSize',[14])
view([16.9, 36.4])
print fit_linear_trend_largebox.png -dpng

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elev     = eval_pts(3, :);
northing = eval_pts(2, :);
easting  = eval_pts(1,:);

measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

gz_avg_at_stations = cellfun(@mean, measured_values);
gz_error_at_stations = cellfun(@norm, measure_errors);

gobs = gz_avg_at_stations';
gerr = gz_error_at_stations';

below_cutoff_height = elev < 2150;
UpObs   = gobs(~below_cutoff_height);
nUp     = numel(UpObs);
LowObs  = gobs(below_cutoff_height);
nLow    = numel(LowObs);
xU      = easting(~below_cutoff_height);
yU      = northing(~below_cutoff_height);
xL      = easting(below_cutoff_height);
yL      = northing(below_cutoff_height);

gshift = mean(plane1(:)) - mean(UpObs);
plane2 = plane1 - gshift;
gshift = mean(plane1(:)) - mean(LowObs);
plane3 = plane1 - gshift;

F2 = scatteredInterpolant(X(:),Y(:),plane2(:));
F3 = scatteredInterpolant(X(:),Y(:),plane3(:));
gtrendU2 = F2(xU,yU);
gtrendL2 = F3(xL,yL);

%%
% take average trend and make error-bar
gtrendU = 0.5*(gtrendU1+gtrendU2-mean(gtrendU1)-mean(gtrendU2));
gtrendL = 0.5*(gtrendL1+gtrendL2-mean(gtrendL1)-mean(gtrendL2));
trenderrU = 0.5*abs(gtrendU1-mean(gtrendU1)-gtrendU2+mean(gtrendU2));
trenderrL = 0.5*abs(gtrendL1-mean(gtrendL1)-gtrendL2+mean(gtrendL2));

% add in quadrature to measurement errors
errU = sqrt(gerr(~below_cutoff_height).^2 + trenderrU.^2);
errL = sqrt(gerr(below_cutoff_height).^2 + trenderrL.^2);

%%
%now make output figures

figure(6)
errorbar(xU,UpObs-mean(UpObs),gerr(~below_cutoff_height),'bd','markerfacecolor','b','markeredgecolor','k');
hold on
errorbar(xU, gtrendU-mean(gtrendU),trenderrU,'kd')
errorbar(xU, UpObs -gtrendU-mean(UpObs)+mean(gtrendU),errU,'rd','markerfacecolor','r','markeredgecolor','k');
%plot(xU, gtrendU-mean(gtrendU),'b--')
xlabel('Easting, m'); ylabel('mgals');set(gca,'FontSize',[14])
legend('Observations','Regional trend','De-trended residual')
print Upper_detrend_largebox.png -dpng

figure(7)
errorbar(yL,LowObs-mean(LowObs),gerr(below_cutoff_height),'bd','markerfacecolor','b','markeredgecolor','k');
hold on
errorbar(yL, gtrendL-mean(gtrendL),trenderrL,'kd')
errorbar(yL, LowObs - gtrendL-mean(LowObs)+mean(gtrendL),errL,'rd','markerfacecolor','r','markeredgecolor','k');
%plot(yL, gtrendL-mean(gtrendL),'b--')
xlabel('Northing, m'); ylabel('mgals');set(gca,'FontSize',[14])
legend('Observations','Regional trend','De-trended residual')
print Lower_detrend_largebox.png -dpng

%%
trendU = [xU', yU', gtrendU', trenderrU'];
trendL = [xL', yL', gtrendL', trenderrL'];
save('Trend_Upper.dat','trendU','-ascii');
save('Trend_Lower.dat','trendL','-ascii');


