% INVERT GRAVITY OBSERVATIONS IN LANL MUONGRAVITY PROJECT - MAY 2017
%Modification made by MR -- 5/15-17/17; I tested that the fwd kernel is being
%correctly calculated -- test_single_station_voxel.m -- this shows that the
%gravity effects are equal and opposite when the station is above or
%beneath a voxel... in the middle of the voxel, by symmetry, the predicted
%gz is near zero -- tested the gz function etc, and all is good!
% 
% Now we will use the actual topo of the study area to drive the
% gravity-only inversion for the study area.  Got this working and tested
% the resolution using resolution_test.m

% 5-31-17: so far, only lateral rho variations are allowed.  Need to
% include vertical variations also.
%
% 6-1: got vertical variations working:
% lanl_gravity_inversion_jun_1_uniformguess.m
% 6-4: got layered starting guess working

    close all;
    clear all;
    
%% %Constants
    
    correlation_fac = 0.1;
    xycorr_fac   = 1e2;
    rho_base     = 1800;
    %rho_init     = 1800;
    %set up parameters for an initial rho_guess that is a geology-based
    %3-layer model
    rho_Qbt1     = 1420;
    rho_Qbt2     = 2100;
    rho_Qbt3     = 1850;
    up_elev      = 2183; %upper layer interface elevation
    low_elev     = 2150; %must be higher than min_z+dz, 2150 or so
    
    rho_ht       = rho_base*0.1;
    rho_std      = rho_ht*0.1;
       
    niter=50;
    mu = 0.09; % step size in the iterations
    
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
% DISCOVERED THAT you need to pad by around 465 m around the box enclosing
% all stations such that with a dx=dy=15 to 25m, we get same results as for the 
% ref. discretization used in forward models of 400*400 nodes.

    dx = 16;  % use 15 m spacing of voxels in inner region
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
%     xnodes = [X1in-465:dx:X2in+465]; %x-positions of voxels
%     ynodes = [Y1in-465:dy:Y2in+465]; %y-positions of voxels
%     num_voxels = length(xnodes)*length(ynodes)
    
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
    
%% Geometry: make voxels based on the optimized discretization above
for ilayer = 1:maxlayers
    layerbot     = (min_z + dz*(ilayer - 1))*ones(size(ElevI));
    layertop     = (min_z + dz*ilayer)*ones(size(ElevI));
    layerelev    = (min_z + dz*ilayer);
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
        rho_0 = [rho_Qbt1.*ones(1,length(allow_x))];
    else
        voxel_corner = [voxel_corner bot_corner];
        voxel_top_corner = [voxel_top_corner top_corner];
        if layerelev < low_elev
            rho_0 = [rho_0 rho_Qbt1.*ones(1,length(allow_x))];
        elseif layerelev < up_elev
            rho_0 = [rho_0 rho_Qbt2.*ones(1,length(allow_x))];
        else
            rho_0 = [rho_0 rho_Qbt3.*ones(1,length(allow_x))];
        end
    end
    
    
end
voxel_diag = voxel_top_corner - voxel_corner;
voxel_cen  = voxel_corner+voxel_diag*0.5;
voxel_sep = zeros(num_voxels);
for i=1:num_voxels
    for j=1:num_voxels
        voxel_sep(i,j)    = norm( voxel_cen(:,i)   - voxel_cen(:,j)   );
        voxel_sep_xy(i,j) = norm( voxel_cen(1:2,i) - voxel_cen(1:2,j) );
        voxel_sep_z(i,j)  =  abs( voxel_cen(3,i)   - voxel_cen(3,j)   );
    end
end
    
%% Make G
tic
    G = create_interaction_matrix(eval_pts, voxel_corner, voxel_diag);
toc;

%% test out the forward calc first
    true_rho    = rho_base*ones([length(voxel_corner(1,:)) 1]);
    rho_avg     = mean(true_rho);
    true_gz     = G * true_rho;
    
    %predicted g_z relative to base station
    gz_vals        = true_gz;
    offset_gz_vals = (gz_vals - gz_vals(strcmp(measured_points, 'BS_TN_1'))) * 1E5;

    
    %% make plot of stations and voxels
    figure(1); clf; hold on;
    
    scatter3(voxel_corner(1,:), voxel_corner(2,:), voxel_corner(3,:), 'ro');
    scatter3(voxel_top_corner(1,:), voxel_top_corner(2,:), voxel_top_corner(3,:), 'b*');
    scatter3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:));
    xlabel('Easting, m')
    ylabel('Northing, m')
 
    figure(2);clf
    scatter3(eval_pts(1,:), eval_pts(2,:), offset_gz_vals, 'ro'); hold on
    scatter3(eval_pts(1,:), eval_pts(2,:), measured_gz, 'b*');
    legend('true\_gz', 'measured\_gz')
    xlabel('Easting, m')
    ylabel('Northing, m')
        
%%    %Make data covar
    covar_data     = 1e-10*diag(gz_error_at_stations).^2;
    covar_data_inv = inv(covar_data);
    
    %Make prior model/rho covar, includes SMOOTHING   
    boxlength = Xrange;
    correlation_length = correlation_fac.*boxlength;
    xycorrelation_length = xycorr_fac*correlation_length;
    covar_rho      = (rho_std^2)*exp(-voxel_sep_z/correlation_length-voxel_sep_xy/xycorrelation_length);
    inv_covar_rho  = inv(covar_rho);
    
    %Make posterior model covariance
    inv_C_post     = G.' * (covar_data\G) + inv_covar_rho;
    %C_post         = inv(inv_C_post);
    operator1      = (inv_C_post\G');
    %operator2      = C_post * inv_covar_rho;
    
    %Iterate with rho_inv as new input model
    
    %rho_0 = repmat(rho_init, num_voxels, 1);
    rho_0  = rho_0';
    rho_guess = rho_0;                 
    for i=1:niter 
            pred_gz   = G * rho_guess;
            residual  = pred_gz - measured_gz*1e-5;  %note that the measured_gz is in mgals
            
            figure(2); scatter3(eval_pts(1,:), eval_pts(2,:), pred_gz*1e5, 'gd'); hold on
            term1     = operator1*(covar_data\residual);
            term2     = inv_C_post\(covar_rho\(rho_guess - rho_0));
            term      = term1 + term2;
            rho_inv   = rho_guess - mu*term;
            rho_old   = rho_guess;
            rho_guess = rho_inv;
            
            rmsrho     = sqrt(sum(((rho_guess - true_rho).^2 ))/num_voxels);
            rmsgz      = rms(residual);
            rmslist(i) = rmsgz;
            

%will not work in 3D 
%             figure(3); clf
%             subplot(211); hold on;
%             surf(XI, YI, reshape(rho_old, size(XI)), 'FaceAlpha', 0.5, 'FaceColor', [0;1;0]);
%             surf(XI, YI, reshape(rho_inv, size(XI)), 'FaceAlpha', 0.5, 'FaceColor', [1;0;0]);
%             %surf(x, y, reshape(rho_old, size(x)), 'FaceAlpha', 0.5, 'FaceColor', [0;0;1]);
%             view([30, 30])
%             legend('true\_rho', 'rho\_inv');
%             lighting gouraud
          
            figure(4)
            subplot(211);plot([1:i], rmslist,'o');title('gz_rms')
            subplot(212);plot(i, rmsrho,'o');hold on;title('rho_rms')

%will not work in 3D 
%             figure(5); 
% %             station_true_rho = interp2(XI, YI, reshape(true_rho, size(XI)), eval_pts(1,:), eval_pts(2,:), 'linear');
%             station_inv_rho  = interp2(XI, YI, reshape(rho_inv, size(XI)), eval_pts(1,:), eval_pts(2,:), 'linear');
%             %plot(station_true_rho,'ro');hold on
%             hold on; plot(station_inv_rho,'k+');plot(station_inv_rho,'k--'); 
%             
            disp('paused in iteration loop ')
            pause(2)
            
    end

disp('Inversion Complete - Now Run inv_post_proc.m')