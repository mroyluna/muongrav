function forward_calculation(n,enable_plotting)
%FORWARD_CALCULATION Terrain based forward model gravity calculation
close all;
if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

load Trend_Lower.dat
load Trend_Upper.dat
TU = Trend_Upper(:,3);% - mean(Trend_Upper(:,3));
TL = Trend_Lower(:,3);% - mean(Trend_Lower(:,3));

%find best fit to both upper and lower stations

%read data from southern box 5079 rows by 4600 columns of points
[soX, soY, soElev] = BinaryTerrain.read_file('Tunnel_points_20160715.bin', 23363400, 4600, 5079);
%read data from northern box
[noX, noY, noElev] = BinaryTerrain.read_file('LIDARDATA_20161116.bin', 25091424, 6168, 4068);
 
%pick out region to work with:
%easting edges: 
edgex =[min(soX(:)), max(soX(:))];
edgey =[min(soY(:)), max(soY(:))];
% %northing edges: 
% edgey =[5.3993, 5.4147]*1e5;
[closestx,ind1] = min(abs(noX(1,:)-edgex(1)));
[closesty,ind2] = min(abs(noY(:,1)-edgey(2)));
 
%pick out the part of the N box that has the same x-width as the S box
cut1X = noX(:,ind1:ind1+4599);
cut1Y = noY(:,ind1:ind1+4599);
cut1Elev = noElev(:,ind1:ind1+4599);
%pick out the part of the N box that is N of the S box
cut2X = cut1X(1:ind2-1,:);
cut2Y = cut1Y(1:ind2-1,:);
cut2Elev = cut1Elev(1:ind2-1,:);
 
%stack the two boxes of data in order (south below north, but each matrix 
%cut2Y and soY has decreasing y values in the columns)
X = [cut2X;soX]; 
Y = [cut2Y;soY];
Elev = [cut2Elev;soElev];

% Generate submesh on which to downsample
[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));

% Interpolate linearly since we're downsampling
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');

%Constants
 rock_density       = 1900;
 delta_rock_density = 0;
% %delta rock density is used in the lower parts of the model
 LayerElev=2160;

dx = XI(1,2) - XI(1,1);
dy = YI(2,1) - YI(1,1);
assert(dx > 0 & dy > 0)

min_z = min(ElevI(:));
max_z = max(ElevI(:));

rho = repmat (rock_density, n*n,1);

rhoL = repmat(delta_rock_density, n*n, 1);

% eval_pts = [Constants.base_station, Constants.tunnel_pts];

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elevations = eval_pts(3, :);
northing = eval_pts(2, :);
easting = eval_pts(1, :);

voxel_corners = [XI(:)'; YI(:)'; repmat(min_z, 1, n*n)];

Layer=repmat(LayerElev, 1 , n*n);
comparelayer = [Layer; ElevI(:)'];
minElev=ones(1,n*n);
minElev = min(comparelayer);

eastingcut=50
%number 1 through n determines how far east the cut displayed in the
%next figure is
% clf
% figure(1)
% plot(YI(:,eastingcut),Layer(:,eastingcut))
% hold on
% plot (YI(:,eastingcut),min_z*(ones(n,1)))
% plot (YI(:,eastingcut), ElevI(:,eastingcut))
% plot (YI(:,1),minElev((eastingcut*n)-(n-1):eastingcut*n)','o')
% plot(northing, elevations,'o')

voxel_diag = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - min_z];
%voxel_diag_high = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - LayerElev];
voxel_diag_low= [repmat(dx, 1, n*n); repmat(dy, 1, n*n); minElev - min_z ];

tic;
interaction_matrixL = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag_low);
interaction_matrix = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag); 
toc;

lc = point_table{'W wall tunnel', Constants.xyz_index}';

tunnel_rooms = tunnel_spec(lc, Constants.tunnel_angle_offset_from_north, Constants.tunnel_slope);

ind = 1;
for pt = eval_pts
    for p_id = 1:4
        tunnel_effect(ind, p_id) = tunnel_rooms(p_id).eval_gz_at(pt);
    end
    ind = ind + 1;
end

rho_oriented = repmat(-(rock_density+delta_rock_density), numel(tunnel_rooms), 1);

gz_vals = interaction_matrixL * rhoL + interaction_matrix * rho + tunnel_effect * rho_oriented;
gz_vals_no_tunnel = interaction_matrixL * rhoL + interaction_matrix * rho;
%gz_vals = interaction_matrix * rho + tunnel_effect * rho_oriented;

offset_gz_vals = (gz_vals - gz_vals(strcmp(measured_points, 'BS_TN_1'))) * 1E5;
offset_gz_vals_no_tunnel = (gz_vals_no_tunnel - gz_vals_no_tunnel(strcmp(measured_points, 'BS_TN_1'))) * 1E5;
%tunnel_gz = offset_gz_vals - offset_gz_vals_no_tunnel;
tunnel_gz = tunnel_effect * rho_oriented * 1e5;

% if ~enable_plotting
%     return
% end

% Increase the fineness of the default color map
set(0, 'DefaultFigureColormap', parula(1024*16));


measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

gz_avg_at_stations = cellfun(@mean, measured_values);
gz_error_at_stations = cellfun(@norm, measure_errors);

below_cutoff_height = elevations < 2150;


   function do_plot(fig_num, fig_name, mask, N, gz_calc, gz_calc2, gz_meas, gzT, error)
        figure(fig_num); hold on;
        scatter(N(mask), gz_calc(mask),'rd')
        %scatter(N(mask), gz_calc2(mask),'bd')
        %remove trend
        gz   = gz_meas(mask) - gzT + mean(gz_meas(mask));

        %errorbar(N(mask), gz_meas(mask), error(mask), 'go','markerfacecolor','g');
        errorbar(N(mask), gz, error(mask), 'bo','markerfacecolor','b');
        
        title([fig_name ' n = ' num2str(n) ', top density = ' num2str(rock_density) ', delta density = ' num2str(delta_rock_density) ', layer elevation = ' num2str(LayerElev)]);
        legend('Calculated',  'Observed');
        %legend('Calculated', 'Calculated_no_tunnel', 'Observed');
        %xlabel('Northing (m)'); ylabel('gz (mgal)');
        % saveas(gcf, ['figures/' fig_name ' stations_' num2str(n) '_' num2str(int64(Constants.rock_density))], 'png');
   end

do_plot(10, 'Lower',  below_cutoff_height, northing, offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, TL, gz_error_at_stations);
xlabel('Northing (m)'); ylabel('gz (mgal)');
figure(10)
filename=['N_Lower_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.png'];
%eval(['print ./images/' filename ' -dpng']);

do_plot(11, 'Upper', ~below_cutoff_height, northing, offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, TU, gz_error_at_stations);
xlabel('Northing (m)'); ylabel('gz (mgal)');
figure(11)
filename=['N_Upper_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.png'];
%eval(['print ./images/' filename ' -dpng']);


do_plot(12, 'Lower',  below_cutoff_height, easting, offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, TL, gz_error_at_stations);
xlabel('Easting (m)'); ylabel('gz (mgal)');
figure(12)
filename=['E_Lower_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.png'];
%eval(['print ./images/' filename ' -dpng']);


figure(13); hold on
do_plot(13, 'Upper', ~below_cutoff_height, easting, offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, TU, gz_error_at_stations);
%input in the gravity of a buried cylinder (shape only)
tunx = 495740;
xloc = easting(~below_cutoff_height) - tunx;
tunh = 1.83; %tunnel is 3.66 m across; about 12 feet - tunnel height
tund = 100; %meters depth
grav_cyl = 1e5*(Constants.G*2*pi*(-rock_density-delta_rock_density)*tund*(tunh)^2)./(xloc.*xloc + tund*tund);
hold on;
%plot(easting(~below_cutoff_height),mean(gz_avg_at_stations(~below_cutoff_height))+grav_cyl,'mo');
%plot(easting(~below_cutoff_height),grav_cyl,'ko');
%plot(easting(~below_cutoff_height),grav_cyl,'m--','linewidth',[2]);
%plot(easting(~below_cutoff_height),mean(gz_avg_at_stations(~below_cutoff_height))+tunnel_gz(~below_cutoff_height),'ko');
plot(easting(~below_cutoff_height),tunnel_gz(~below_cutoff_height)+mean(gz_avg_at_stations(~below_cutoff_height)),'ko');%,'markerfacecolor','r','markersize',[10]);
%plot([tunx, tunx],[min(grav_cyl),max(grav_cyl)],'k--');
plot([tunx, tunx],[min(gz_avg_at_stations(~below_cutoff_height)),max(gz_avg_at_stations(~below_cutoff_height))],'k--');
xlabel('Easting (m)'); ylabel('gz (mgal)');title('Gravity Effect of Tunnel')
filename=['E_Upper_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.png'];
%filename=['TunnelGravity.png'];
eval(['print ' filename ' -dpng']);


do_plot(14, 'Lower',  below_cutoff_height, elevations, offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, TL, gz_error_at_stations);
xlabel('Elevation (m)'); ylabel('gz (mgal)');
figure(14)
filename=['El_Lower_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.png'];
%eval(['print ./images/' filename ' -dpng']);

do_plot(15, 'Upper', ~below_cutoff_height, elevations, offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, TU, gz_error_at_stations);
xlabel('Elevation (m)'); ylabel('gz (mgal)');
figure(15)
filename=['El_Upper_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.png'];
%eval(['print ./images/' filename ' -dpng']);

dataarray=[easting', northing', elevations', offset_gz_vals, offset_gz_vals_no_tunnel, gz_avg_at_stations, gz_error_at_stations];
filename=['AllData_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.dat'];
%save(filename, 'dataarray','-ascii', '-double');
%%
figure(2); hold on;
title('Elevation Data and Station Locations');
xlabel('Easting (m)'); ylabel('Northing (m)');

surf(XI, YI, ElevI, 'EdgeAlpha', 0.15);
scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:) + 2, 10, 'o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');

%for i=1:10:n*n
    %render_prism(voxel_corners(:,i),voxel_diag_low(:,i),[1;0;0],[0;1;0],[0;0;1],'red');
%end
hold on;
diag_pts = voxel_corners + voxel_diag_low;
vox_height = diag_pts(3,:);
%height = min(vox_height(:),ElevI(:));
scatter3(diag_pts(1,1:2:end),diag_pts(2,1:2:end),diag_pts(3,1:2:end),'k.');
set(gca,'xlim',[4.955,4.958]*1e5,'ylim',[5.408, 5.416]*1e5)
set(gca,'fontsize',[14])
axis equal tight
lighting gouraud

end


    
function test_rrpa
    x = -150:1:150;
    y = -150:1:150;
    [X, Y] = meshgrid(x,y);

    corner = [-100; -100; -200];
    diagonal = [200; 200; 100];

    m = create_interaction_matrix([X(:)'; Y(:)'; 0 * X(:)'], corner, diagonal);

    calc_gz = reshape(m * 2000, [length(y), length(x)]) * 1E5;

    surf(X, Y, calc_gz, 'EdgeColor', 'none'); hold on;
    contour3(X, Y, calc_gz, 'k');
    axis equal
end