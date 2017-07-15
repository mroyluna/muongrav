% make elevation profile around tunnel
% Mousumi Roy -- 12/5/16
clear all

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

n=400;
% Generate submesh on which to downsample
[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));

% Interpolate linearly since we're downsampling
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');

 
[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elevations = eval_pts(3, :);
northing = eval_pts(2, :);
easting = eval_pts(1, :);

figure(1); hold on;
title('Elevation Data and Station Locations');
xlabel('Easting (m)'); ylabel('Northing (m)');

surf(XI, YI, ElevI, 'EdgeAlpha', 0.15);
scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:) + 2, 10, 'o', ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');

set(gca,'fontname','Helvetica','fontsize',[14])

% pikc out a region within 40 m of   tunnel

tunx = 495740;
mask = XI >= tunx-20 & XI <= tunx+20;
Xm = XI(mask);
Ym = YI(mask);
Elevm = ElevI(mask);

figure(2); subplot(211); hold on
title('Lidar within \pm 20 m of Tunnel and Station Locations');
xlabel('Northing (m)');

plot(Ym, Elevm, 'k.');
plot(eval_pts(2,:)', eval_pts(3,:)', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
plot([5.405,5.418]*1e5, [2183 2183], 'r--');
set(gca,'xlim',[5.405,5.418]*1e5);
set(gca,'fontname','Helvetica','fontsize',[14])
%print 'elev20mprofile.png' -dpng

% pikc out a region within 50 m of   tunnel
tunx = 495740;
mask = XI >= tunx-50 & XI <= tunx+50;
Xm = XI(mask);
Ym = YI(mask);
Elevm = ElevI(mask);

subplot(212); hold on;
title('Lidar within \pm 50 m of Tunnel and Station Locations');
xlabel('Northing (m)');

plot(Ym, Elevm, 'k.');
plot(eval_pts(2,:)', eval_pts(3,:)', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
plot([5.405,5.418]*1e5, [2183 2183], 'r--');
set(gca,'xlim',[5.405,5.418]*1e5);
set(gca,'fontname','Helvetica','fontsize',[14])
%print 'elev20mprofile.png' -dpng


% pikc out a region within 100 m of   tunnel
tunx = 495740;
mask = XI >= tunx-100 & XI <= tunx+100;
Xm = XI(mask);
Ym = YI(mask);
Elevm = ElevI(mask);

% subplot(313); hold on;
% title('Lidar within \pm 50 m of Tunnel and Station Locations');
% xlabel('Northing (m)');
% 
% plot(Ym, Elevm, 'k.');
% plot(eval_pts(2,:)', eval_pts(3,:)', 'o','MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
% set(gca,'xlim',[5.405,5.418]*1e5);
% set(gca,'fontname','Helvetica','fontsize',[14])
print 'elevprofile.png' -dpng

