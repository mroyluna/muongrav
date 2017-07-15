%gridsearchgeo.m
%grid  search to find the height of the lower interface in a three layer
%model based on known geology (from the well logs)
close all;

if ~exist('enable_plotting', 'var')
    enable_plotting = true;
end

% [X, Y, Elev] = read_binary_file('TA41_Tunnel_LIDAR_NAVD88.bin', 2422420, 1540, 1573);
%[X, Y, Elev] = BinaryTerrain.read_file('Tunnel_points_20160715.bin', 23363400, 4600, 5079);

% [X, Y, Elev] = read_binary_file('TA41_Tunnel_LIDAR_NAVD88.bin', 2422420, 1540, 1573);
%[X, Y, Elev] = read_binary_file('TA41_Tunnel_LIDAR_NAVD88.bin', 2422420, 1540, 1573);
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

%smaller box 1
% X = X(1095:6266,:);
% Y = Y(1095:6266, :);
% Elev = Elev(1095:6266, :);

%smaller box 2
%X=X(300:6266,:);
%Y=Y(300:6266,:);
%Elev=Elev(300:6266,:);


%LayerSlopeN= -10
%LayerSlopeE=1

LayerSlopeNU= -4:1:0;
LayerSlopeEU= -4:1:0;
LayerSlopeNL= -10:1:-6;
LayerSlopeEL=  6:1:10;
LayerSlopeNU=tand(LayerSlopeNU);
LayerSlopeEU=tand(LayerSlopeEU);
LayerSlopeNL=tand(LayerSlopeNL);
LayerSlopeEL=tand(LayerSlopeEL);
n=400;

for  i=1:length(LayerSlopeNU)
    for j=1:length(LayerSlopeEU)
        for k=1:length(LayerSlopeNL)
            for l=1:length(LayerSlopeEL)
                [XI, YI] = meshgrid(linspace(min(X(:)), max(X(:)), n), ...
                    linspace(min(Y(:)), max(Y(:)), n));
                % Interpolate linearly since we're downsampling
                
                ElevI = interp2(X, Y, Elev, XI, YI, 'linear');
                LayerElev=2183;
                height=2130;
                pivotN=541000;
                pivotE=XI(floor((n^2)/2));
                Layer2=(YI-pivotN).*LayerSlopeNL(k) + (XI-pivotE).*LayerSlopeEL(l) + height;
                Layer=(YI-pivotN).*LayerSlopeNU(i) + (XI-pivotE).*LayerSlopeEU(j)+LayerElev;
                if isempty(find(Layer2 > Layer))
                    disp(['i, j, k, l', num2str(i),num2str(j),num2str(k),num2str(l)])
                    forward_calculation_geo_septilt(400, LayerSlopeNU(i),...
                       LayerSlopeEU(j), LayerSlopeNL(k), LayerSlopeEL(l),XI, YI, ElevI, X, Y, Elev, Layer, Layer2)
                    i
                end
		[i j k l]
            end
        end
    end
end
