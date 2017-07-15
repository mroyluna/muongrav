%gridsearch.m
% modified by MR on 12/1/16

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


rd  =[500:50:1600];
%drd =[-1800:200:1800];
drd =[0];
%le  =[2130:20:2210];
le  =[2130];

% rd  =1900
% drd=200
% le=2163

 for k=(1:length(le))
   count = 1;
   for i=(1:length(rd))
    for j=(1:length(drd))
        %if (rd(i)+drd(j) <= 2800) && (rd(i)+drd(j) >= 1200)
            forward_calculation_mod(400, rd(i), drd(j), le(k), X, Y, Elev);
        %else
        %    disp('out of range')
        %end
    end
   end
end
