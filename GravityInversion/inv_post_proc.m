%Post_processing of Inversion results -- run after lanl_gravity_inversion_junX.m
%Outputs
      
%will not work in 3D 
%     figure(6); clf;hold on;
%     for ind=1:num_voxels
%         color = (rho_inv(ind) / max(true_rho(:)))^4 * [1;0;0];
%         if max(color(:))>1
%             color(color>1) = 1;
%         end
%        
%        render_prism_orig(voxel_corner(:, ind), voxel_diag(:, ind), [1;0;0], [0;1;0], [0;0;1], color);
%     end
%     plot3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:)+2,'k.','markersize',[18]);
%     xlabel('Easting, m')
%     ylabel('Northing, m')
%     title(['Inverted density structure, L_{corr}= ' num2str(correlation_length) ' m, rho_init= ', num2str(rho_init)])
%     %for only inner region use
%     set(gca,'xlim',[X1in, X2in],'ylim',[Y1in, Y2in])
%--------------------------------------------ADAPTED OUTPUT FOR 3D------------------------   
% Need to come up with a clever way to visualize the density structure in 3D
% making slices, etc.

% use existing XI YI outputs from meshgrid to make arrays that fill a 3D
% space:

zreg = [min_z:dz:max_z];
%xreg = xnodes;
%yreg = ynodes;
xreg = [min(xnodes):dx:max(xnodes)];
yreg = [min(ynodes):dy:max(ynodes)];

[X3, Y3, Z3] = ndgrid(xreg', yreg', zreg');
FRHO    = scatteredInterpolant(voxel_cen(1,:)', voxel_cen(2,:)', voxel_cen(3,:)', rho_inv);
Intrho  = FRHO({xreg, yreg, zreg});
%comment the following out if xreg and yreg are not xnodes and ynodes
%Elevreg = ElevI';
%uncomment the following only if xred and yred are not xnodes and ynodes
FElev   = scatteredInterpolant(XI(:), YI(:), ElevI(:));
Elevreg = FElev({xreg,yreg});

%make a 3D topo mask that is 0 at pts outside the rock volume and 1 inside
Elev3d = repmat(Elevreg,1,1,length(zreg));
topomask = (Z3 <= Elev3d+dz);
%set density outside rock volume to be zero
Intrho(~topomask) = NaN;

%permute matrices to wrap in y first before x
XB = permute(X3,[2 1 3]);
YB = permute(Y3,[2 1 3]);
ZB = permute(Z3,[2 1 3]);
RhoB = permute(Intrho,[2 1 3]);
Elevreg = permute(Elevreg, [2 1 3]);

%now can make output figures
figure(6); clf
%ax1 = subplot(211);
%colormap(ax1, winter);
colormap(winter);
mesh(squeeze(XB(:,:,1)), squeeze(YB(:,:,1)), Elevreg);hold on
plot3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:)+2,'k.','markersize',[18]);
xlabel('Easting, m')
ylabel('Northing, m')
zlabel('Elevation, m')

figure(7);clf; hold on
%ax2 = subplot(212);
%colormap(ax2, flipud(pink));
cvals=[1300:200:2200];
%define slice x, y , z point
sx1 = [0.5*(X2in + X1in)-30];
sx2 = [0.5*(X2in + X1in)-200];
sx3 = [0.5*(X2in + X1in)+110];
sy = [0.5*(Y2in + Y1in)+150];
sz = [0.5*(max_z + min_z)];
sz = min_z;

hx1 = slice(XB,YB,ZB,RhoB,sx1,[],[]);
hx1.FaceColor = 'interp';
%hx1.DiffuseStrength = 0.8;
hx1.EdgeColor = [0.5, 0.5 0.5];
hx1.EdgeAlpha = [0.5];

%hx2 = slice(XB,YB,ZB,RhoB,sx2,[],[]);
%hx2.FaceColor = 'interp';
%hx2.DiffuseStrength = 0.8;

hx3 = slice(XB,YB,ZB,RhoB,sx3,[],[]);
hx3.FaceColor = 'interp';
%hx3.DiffuseStrength = 0.8;
hx3.EdgeColor = [0.5, 0.5 0.5];
hx3.EdgeAlpha = [0.5];

%hy = contourslice(XB,YB,ZB,RhoB,[],sy,[],cvals);
%hy.FaceColor = 'interp';
hy = slice(XB,YB,ZB,RhoB,[],sy,[]);
hy.FaceColor = 'interp';
%hy.DiffuseStrength = 0.8;
hy.EdgeColor = [0.5, 0.5 0.5];
hy.EdgeAlpha = [0.5];

hz = slice(XB,YB,ZB,RhoB,[],[],sz);
hz.FaceColor = 'interp';
%hz.DiffuseStrength = 0.8;
hz.EdgeColor = [0.5, 0.5 0.5];
hz.EdgeAlpha = [0.5];

stations = plot3(eval_pts(1,:), eval_pts(2,:), eval_pts(3,:),'o');
stations.MarkerSize=[6];
stations.MarkerEdgeColor='k';
stations.MarkerFaceColor='y';

lightangle(-45,25);
ndivs = 7;
crange=linspace(1900,2400,ndivs);
caxis([crange(1) crange(end)]);
colormap(flipud(pink(ndivs)));
set(gca,'xlim',[X1in X2in],'ylim',[Y1in,Y2in])
colorbar('vert')
xlabel('Easting, m')
ylabel('Northing, m')
zlabel('Elevation, m')
%title(['Inverted Density Structure with Uniform Initial Guess = ', num2str(rho_init)]);
title(['Inverted Density Structure with Layered Initial Guess']);
%filename = ['UniformInitRho' num2str(rho_init),'xycorrfac', num2str(xycorr_fac), '.fig'];
filename = ['LayeredInitRho_xycorrfac', num2str(xycorr_fac),'upelev',num2str(up_elev), '.fig'];
savefig(filename)    