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
eastingcut=1:3:400;
eastingsofcuts=ones(1,length(eastingcut));for i=1:length(eastingcut)
%number 1 through n determines how far east the cut displayed in the
%next figure is
% clf
% hold on
% plot (YI(:,eastingcut(i)), ElevI(:,eastingcut(i)))
% hold on
% plot(northing, elevations,'o')
% plot(YI(:,eastingcut(i)), Layer(:,eastingcut(i)))
% plot(YI(:,eastingcut(i)), Layer2(:,eastingcut(i)))
currenteasting=XI(1,eastingcut(i));
eastingsofcuts(i)=currenteasting;
%set(gca, 'XTickLabel', num2str(get(gca, 'XTick')', '%.1f'));
%h=figure;

%saveas(h,sprintf('hillsection%d.fig',eastingcut(i)));
end
