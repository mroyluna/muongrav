%megan jan 2017
%make figures of cross sections of the cliff to identify the upper layer boundary by eye 
n=400;
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

[XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));
% Interpolate linearly since we're downsampling
                
ElevI = interp2(X, Y, Elev, XI, YI, 'linear');
 %Xn = linspace(min(X(:)), max(X(:)), n)
 %Yn = linspace(min(Y(:)), max(Y(:)), n)
 
%Constants

%delta rock density is used in the lower parts of the model

%lowest point of the top bottom interface
%LayerSlopeN=0;
%LayerSlopeE=0;
%Elevation/northing

dx = XI(1,2) - XI(1,1);
dy = YI(2,1) - YI(1,1);
assert(dx > 0 & dy > 0)

min_z = min(ElevI(:));
max_z = max(ElevI(:));



% eval_pts = [Constants.base_station, Constants.tunnel_pts];

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elevations = eval_pts(3, :);
northing = eval_pts(2, :);
easting = eval_pts(1,:);

%Calculate Error in Elevation at station points using gradient



voxel_corners = [XI(:)'; YI(:)'; repmat(min_z, 1, n*n)];


%comparelayer = [Layer; ElevI];



%number 1 through n determines how far east the cut displayed in the
%next figure is
clf
figure(1)
%eastign of 495230
eastingcut=50
hold on

plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
set(gca, 'XTickLabel', num2str(get(gca, 'XTick')', '%.8f'));
%plot(northing, elevations,'o')

figure(2)
%easting of 495400
eastingcut=100
hold on
plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
%plot(northing, elevations,'o')


figure(3)
%easting of 495580
eastingcut=150
hold on
plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
%plot(northing, elevations,'o')




figure(5)
%easting of 495930
eastingcut=250
hold on
plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
%plot(northing, elevations,'o')

figure(6)
%easting 496110 
eastingcut=300
hold on
plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
%plot(northing, elevations,'o')

figure(7)
%easting of 496280
eastingcut=350
hold on
plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
%plot(northing, elevations,'o')


figure(8)
%easting of 496460
eastingcut=400
hold on
plot (YI(:,eastingcut),min_z*(ones(n,1)))
plot (YI(:,eastingcut), ElevI(:,eastingcut))
%plot(northing, elevations,'o')