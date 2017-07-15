function forward_calculation_geo_tilt(n, LayerSlopeN, LayerSlopeE, X, Y, Elev, enable_plotting)

if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));
% Interpolate linearly since we're downsampling
                
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');
 %Xn = linspace(min(X(:)), max(X(:)), n)
 %Yn = linspace(min(Y(:)), max(Y(:)), n)
 
%Constants
rock_density=1850;
delta_rock_density=250;
delta_rock_density2=-680;
%delta rock density is used in the lower parts of the model
LayerElev=2183;
height=2130;
%lowest point of the top bottom interface
%LayerSlopeN=0;
%LayerSlopeE=0;
%Elevation/northing

dx = XI(1,2) - XI(1,1);
dy = YI(2,1) - YI(1,1);
assert(dx > 0 & dy > 0)

min_z = min(ElevI(:));
max_z = max(ElevI(:));

rho = repmat (rock_density, n*n,1);

rhoL = repmat(delta_rock_density, n*n, 1);

rhoLr = repmat(delta_rock_density2, n*n, 1);

% eval_pts = [Constants.base_station, Constants.tunnel_pts];

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elevations = eval_pts(3, :);
northing = eval_pts(2, :);
easting = eval_pts(1,:);

%Calculate Error in Elevation at station points using gradient

[gradX,gradY]=gradient(Elev,(max(max(X))-min(min(X)))/6266,(max(max(Y))-min(min(Y)))/4600);
grad=sqrt(gradX.*gradX+gradY.*gradY);
stationgrad=interp2(X, Y, grad, easting, northing);
ElevErr=.609*stationgrad;
freeairerr=.3086*2*ElevErr;

voxel_corners = [XI(:)'; YI(:)'; repmat(min_z, 1, n*n)];

pivotN=541000
pivotE=XI(floor((n^2)/2))
%Layer2=(YI-pivotN).*LayerSlopeN + (XI-pivotE).*LayerSlopeE + height;
Layer2= height*ones(size(ElevI));
Layer=(YI-pivotN).*LayerSlopeN + (XI-pivotE).*LayerSlopeE+LayerElev;
%comparelayer = [Layer; ElevI];

comparelayer = [Layer(:)'; ElevI(:)'];
comparelayer2= [Layer2(:)'; ElevI(:)'];
minElev=ones(1,n*n);
minElev2=ones(1,n*n);
minElev = min(comparelayer);
minElev2 = min(comparelayer2);

voxel_diag = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - min_z];
%voxel_diag_high = [repmat(dx, 1, n*n); repmat(dy, 1, n*n); ElevI(:)' - LayerElev];
voxel_diag_low= [repmat(dx, 1, n*n); repmat(dy, 1, n*n); minElev - minElev2 ];
voxel_diag_lower= [repmat(dx, 1, n*n); repmat(dy, 1, n*n); minElev2 - min_z ];
tic;
interaction_matrixLr = create_interaction_matrix(eval_pts, voxel_corners, voxel_diag_lower);
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

rho_oriented = repmat((-rock_density+delta_rock_density+delta_rock_density2), numel(tunnel_rooms), 1);

gz_vals = interaction_matrixL * rhoL + interaction_matrix * rho + interaction_matrixLr * rhoLr +tunnel_effect * rho_oriented;
%gz_vals = interaction_matrix * rho + tunnel_effect * rho_oriented;

% inverse = interaction_matrix \ gz_vals;
% diff = sum(abs(inverse - rho)./rho) / numel(rho)

offset_gz_vals = (gz_vals - gz_vals(strcmp(measured_points, 'BS_TN_1'))) * 1E5;

if ~enable_plotting
    return
end

% Increase the fineness of the default color map
set(0, 'DefaultFigureColormap', parula(1024*16));

%elevations = eval_pts(3, :);
%northing = eval_pts(2, :);
%easting = eval_pts(1,:);

measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

gz_avg_at_stations = cellfun(@mean, measured_values);
gz_error_at_stations = cellfun(@norm, measure_errors);

below_cutoff_height = elevations < 2150;

%     function do_plot(fig_num, fig_name, mask, N, gz_calc, gz_meas, error)
%         figure(fig_num); hold on;
%         scatter(N(mask), gz_calc(mask))
%         errorbar(N(mask), gz_meas(mask), error(mask), 'o');
%         
%         title([fig_name ' n = ' num2str(n) ', bottom density = ' num2str(rock_density) ', delta density = ' num2str(delta_rock_density) ', layer elevation = ' num2str(LayerElev)]);
%         legend(num2str(height), 'Observed');
%         xlabel('Northing (m)'); ylabel('gz (mgal)');
%         % saveas(gcf, ['figures/' fig_name ' stations_' num2str(n) '_' num2str(int64(Constants.rock_density))], 'png');
%     end
% 
% do_plot(10, 'Lower',  below_cutoff_height, northing, offset_gz_vals, gz_avg_at_stations, gz_error_at_stations);

%onslope= gz_calc(below_cutoff_height) > 0
%lowerrors=onslope(freeairerr)

% do_plot(11, 'Upper', ~below_cutoff_height, northing, offset_gz_vals, gz_avg_at_stations, freeairerr);
dataarray=[easting',northing',elevations',offset_gz_vals,gz_avg_at_stations,gz_error_at_stations, LayerSlopeN*ones(size(elevations')),LayerSlopeE*ones(size(elevations'))];
filename=['AllDataTopTilt_', num2str(n),'_', num2str(LayerSlopeN),'_', num2str(LayerSlopeE),'.dat'];
save(filename, 'dataarray', '-ascii', '-double');

UpAveObs=mean(gz_avg_at_stations(~below_cutoff_height));
LowAveObs=mean(gz_avg_at_stations(below_cutoff_height));
UpAveCalc=mean(offset_gz_vals(~below_cutoff_height));
LowAveCalc=mean(offset_gz_vals(~below_cutoff_height));

% rmsUp=(abs(UpAveObs^2-UpAveCalc^2))^(1/2);
% rmsLow=(abs(LowAveObs^2-LowAveCalc^2))^(1/2);
% rmsdata=[n, rock_density, delta_rock_density, LayerElev, rmsUp, rmsLow]
% filename2=['rms_', num2str(n),'_',num2str(rock_density),'_',num2str(delta_rock_density),'_',num2str(LayerElev),'.dat'];
% save(filename2, 'rmsdata', '-ascii', '-double')
%%
% figure(2); hold on;
% title('Elevation Data and Station Locations');
% xlabel('Easting (m)'); ylabel('Northing (m)');
% 
% surf(XI, YI, ElevI, 'EdgeAlpha', 0.15);
% scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:) + 2, 10, 'o', ...
%     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
% 
% diag_pts = voxel_corners + voxel_diag_low;
% scatter3(diag_pts(1,1:2:end),diag_pts(2,1:2:end),diag_pts(3,1:2:end),'r.');
% 
% axis equal tight
% lighting gouraud

% %%
% figure(3); hold on; axis equal;
% tunnel_pts = points_by_regexp(point_table, 'TS[0-9][0-9]');
% 
% scatter3(tunnel_pts(1,:), tunnel_pts(2,:), tunnel_pts(3,:));
% 
% arrayfun(@render, tunnel_rooms);
% 
% [p, X, Y, Z] = fill_plane(lc + [0;0;50], [0;1;0], [0;0;1], [300,300], 500);
% tic;
% 
% inc = 40;
% resolution = gravity_kernel_function(p, eval_pts(:,inc), gz_vals(inc) / Constants.rock_density);
% toc;
% surf(X, Y, Z, log10(reshape(abs(resolution), size(X))), 'EdgeAlpha', 0.1);
% scatter3(eval_pts(1,inc), eval_pts(2,inc), eval_pts(3,inc), 10, 'o', ...
%     'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
% colorbar;
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