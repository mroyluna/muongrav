%CODE TO EXPLORE THE DISCRETIZATION NEEDED TO ACHIEVE SIMILAR PREDICTED
%GRAVITY VALUES TO reference case in forward calcs of n=400 (400*400 voxels)
% but discretizing the volume in 3D -- note, Elena's muon inversion uses
% cubes that are 4 m on a side
% MR -- 6/1/2017
%

    close all;
    clear all;
    
    %% %Constants
    n=400;
    num_voxels   = n;
    num_stations = 60; % we will make two groups of stations at two different average elevations
    boxlength    = 10000; % box length in m
    range        = [0 boxlength]; % box length in m
    extrarange   = [0 boxlength];
    correlation_length = boxlength*0.1;
    rho_base     = 1800;
        
%% first load in and work with LIDAR elevation data as in all the gridsearch codes 
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

    Xrange = max(X(:)) - min(X(:)); % in m
    Yrange = max(Y(:)) - min(Y(:));
    
%% Now load in and work with station locations in the study area, as in other codes
% eval_pts = [Constants.base_station, Constants.tunnel_pts];

[point_table, measured_points] = build_table();

eval_pts = point_table{measured_points, Constants.xyz_index}';
elevations = eval_pts(3, :);
northing = eval_pts(2, :);
easting = eval_pts(1,:);

measured_values = point_table{measured_points, 'Measurements'};
measure_errors = point_table{measured_points, 'Errors'};

gz_avg_at_stations = cellfun(@mean, measured_values);
gz_error_at_stations = cellfun(@norm, measure_errors);
    
measured_gz = gz_avg_at_stations;

%% Now we use the discretization we have tested out in the code 
% optimize_discretization.m -- to make sure that the predicted gravity is
% nearly identical to the forward calculation using 400*400 voxels.
%  Reference discretization, n=400, total number is 400*400
% Divide the study area into an inner and outer region, then use a different
% discretization in each region.  test out different discretrizations - using 
% a coarse outer spacing, DX, DY, and a finer inner spaging, dx, dy, and a 
% division between inner and outer regions, X1in, Y1in, X2in, Y2in

    dx = 20;  % use 15 m spacing of voxels in inner region
    dy = dx;
    DX = 200; % use 200m spacing of voxels in outer region
    DY = DX;
    dz = dx;
    
    X1in = min(eval_pts(1,:)) - 20; % pad inner region containing all stations by 20 m
    X2in = max(eval_pts(1,:)) + 20;
    Y1in = min(eval_pts(2,:)) - 250; % pad by 250 m to south, to get canyon
    Y2in = max(eval_pts(2,:)) + 20;
    
     xnodes = [min(X(:)):DX:X1in-20,X1in:dx:X2in, X2in+20:DX:max(X(:))]; %x-positions of voxels
     ynodes = [min(Y(:)):DY:Y1in-20,Y1in:dy:Y2in, Y2in+20:DY:max(Y(:))]; %y-positions of voxels
    %xnodes = [X1in-200:dx:X2in+200]; %x-positions of voxels
    %ynodes = [Y1in-200:dy:Y2in+200]; %y-positions of voxels
    
    xshifted = wshift('1D', xnodes, 1);
    yshifted = wshift('1D', ynodes, 1);
    xshifted(end) = xshifted(end-1) + DX; %fix last x-value
    yshifted(end) = yshifted(end-1) + DY; %fix last y-value

    [XI, YI] = meshgrid(xnodes, ynodes);
    [XIshifted, YIshifted] = meshgrid(xshifted, yshifted); % will be used for top corner of voxel
    
% Interpolate linearly since we're downsampling
    ElevI = interp2(X, Y, Elev, XI, YI, 'linear');
    
    min_z = min(ElevI(:));
    max_z = max(ElevI(:));
 
% in 3D discretization number of layers we need at each x,y location:
    nlayers    = floor((ElevI-min_z)/dz);
    maxlayers  = max(nlayers(:))
    num_voxels = sum(nlayers(:))

%% %% Geometry: make voxels based on the optimized discretization above

for ilayer = 1:maxlayers
    layerbot     = (min_z + dz*(ilayer - 1))*ones(size(ElevI));
    layertop     = (min_z + dz*ilayer)*ones(size(ElevI));
    %layertop(ElevI < layertop) = ElevI(ElevI < layertop);
    allow_x       = XI(ElevI >= layertop);
    allow_x_shift = XIshifted(ElevI >= layertop);
    allow_y       = YI(ElevI >= layertop);
    allow_y_shift = YIshifted(ElevI >= layertop);
    allow_bot  = layerbot(ElevI >= layertop);
    allow_top  = layertop(ElevI >= layertop);
    bot_corner = [allow_x(:)'; allow_y(:)'; allow_bot(:)'];
    top_corner = [allow_x_shift(:)'; allow_y_shift(:)'; allow_top(:)'];
    if ilayer==1
        voxel_corner = bot_corner;
        voxel_top_corner = top_corner;
    else
        voxel_corner = [voxel_corner bot_corner];
        voxel_top_corner = [voxel_top_corner top_corner];
    end
end
voxel_diag = voxel_top_corner - voxel_corner;
 
voxel_sep = zeros(num_voxels);
for i=1:num_voxels
    for j=1:num_voxels
        voxel_sep(i,j) = norm(voxel_corner(:,i) + voxel_diag(:,i) / 2 - ...
                                  (voxel_corner(:,j) + voxel_diag(:,j) / 2)) ;
    end
end
%%
%Make G
tic
    G = create_interaction_matrix(eval_pts, voxel_corner, voxel_diag);
toc;

%% test out the forward calc 
    true_rho    = rho_base*ones([length(voxel_corner(1,:)) 1]);
    rho_avg     = mean(true_rho);
    true_gz     = G * true_rho;
    rho_std     = sqrt(cov(true_rho));
    
    %predicted g_z relative to base station
    gz_vals        = true_gz;
    offset_gz_vals = (gz_vals - gz_vals(strcmp(measured_points, 'BS_TN_1'))) * 1E5;
 
    %%
    %plot stations and voxels
    figure(1); clf; hold on;
    
    scatter3(voxel_corner(1,:), voxel_corner(2,:), voxel_corner(3,:), 'ro');
    scatter3(voxel_top_corner(1,:), voxel_top_corner(2,:), voxel_top_corner(3,:), 'b*');
    scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:));
    xlabel('Easting, m')
    ylabel('Northing, m')
    

    figure(4); clf
    scatter3(eval_pts(1,:), eval_pts(2,:), offset_gz_vals, 'ro'); hold on
    scatter3(eval_pts(1,:), eval_pts(2,:), measured_gz, 'b*');
    legend('true\_gz', 'measured\_gz')
    xlabel('Easting, m')
    ylabel('Northing, m')
    
