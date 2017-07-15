% analyze_well.m -- MR -- 12/28/16
% used graphclick (a free software package) to digitize the density log
% from the Broxton and Vaniman paper... 
% this code plots the data and analyzes to find the mean densities in the
% intervals corresponding to Qbt3 and Qbt2 and lower units

load BroxtonVanimanDensityLog.dat
rho = BroxtonVanimanDensityLog(:,1);
dep = BroxtonVanimanDensityLog(:,2);
subplot(132);
plot(rho,dep)
set(gca,'ydir','reverse')

% depth of Qbt3/Qbt2 interface is about 175+25/2 ft = 187 ft
% depth of bottom of Qbt2 is about 325 ft
% bottom of low density package is about 815 ft

%rho1 = find(