%feature tilts
%find the tilt of the visable break in the slope of the cliff using data
%taken from the easting cut figures

%this section uses very course measurments, probably not useful anymore
% easting=[495230 495400 495580 495760 495930 496110 496280 496460];
% elev=[2187 2192 2189 2186 2181 2178 2172 2166];
% figure(1)
% plot(easting, elev, 'o')
% hold on
% title('Position of Visible Break in the Cliff')
% xlabel('Easting')
% ylabel('Elevation')
% s=polyfit(easting, elev, 1);
% y=s(1).*easting+s(2);
% plot(easting, y)
% format long
% tilt=atand(s(1))
%%
%do it again with many more measurements
%includes data for the slope of the top of the mesa and the bottom of the
%creek
%MUST run plotmanyeastingcuts before this to load in the lidar and find the
%eastings of each of the measured figures
%excel sheet is set up for n=400 eastingcuts=1:3:400
filename = 'elevations.xlsx';
range = 'A1:C134';
xl = xlsread(filename, range);
mesaelev=xl(:,1);
breakelev=xl(:,2);
creekelev=xl(:,3);
%this comes from the lidar data, must have it loaded in before this
figure(1)
boxx=[495568 495988];
boxy=[2250 2250];
studyarea=area(boxx, boxy); hold on
studyarea.FaceColor=[.9 .9 .9];
studyarea.EdgeColor=[.9 .9 .9];
plot(eastingsofcuts',mesaelev,'k.','markersize',[14])
hold on
plot(eastingsofcuts',breakelev,'k.','markersize',[14])
plot(eastingsofcuts',creekelev,'k.','markersize',[14])
a=polyfit(eastingsofcuts', mesaelev, 1);
b=polyfit(eastingsofcuts', breakelev, 1);
c=polyfit(eastingsofcuts', creekelev, 1);
tiltmesa=atand(a(1))
tiltbreak=atand(b(1))
tiltcreek=atand(c(1))
ya=a(2)+a(1)*eastingsofcuts;
yb=b(2)+b(1)*eastingsofcuts;
yc=c(2)+c(1)*eastingsofcuts;
%plot (eastingsofcuts,ya,'k--'); %MR commented out the polyfit lines
%plot (eastingsofcuts,yb, 'k--');
%plot (eastingsofcuts,yc, 'k--');
%title('Profiles of Mesa, Slope Break, and Creek Bottom Elevations')
ylabel('Elevation, m')
xlabel('Easting, m')
%box on

t1 = text(496100, 2225, {'Mesa',['Slope = -1.211' char(176)]})
t2 = text(496100, 2190, {'Break',['Slope = -1.207' char(176)]})
t3 = text(496100, 2110, {'Creek Bottom',['Slope = -1.656' char(176)]})
t1.FontSize = 16;t2.FontSize = 16;t3.FontSize = 16;
set(gca, 'fontsize', [16])
%%
%add the slope of the interface found in the foward calc

%Layer and ElevI come from forward_calc with the correct tilts inserted
pivotN=541000;
pivotE=XI(floor(n^2/2));
LayerSlopeNU=tand(0);
LayerSlopeEU=tand(5);
LayerSlopeNL=tand(-2);
LayerSlopeEL=tand(4);
Layer=(YI-pivotN).*LayerSlopeNU + (XI-pivotE).*LayerSlopeEU + LayerElev;
Layer2=(YI-pivotN).*LayerSlopeNL + (XI-pivotE).*LayerSlopeEL+height;


[smalldiff, pos]= min(abs(Layer(200:400,:)-ElevI(200:400,:)));
[smalldiff2, pos2]= min(abs(Layer(:,200:400)-ElevI(:,200:400)));


 Layer_on_cliff=ones(1,n);
 
 for i=1:n
     LayerCols=Layer(:,i);
     Layer_on_cliff(i)= LayerCols(pos(i)+199);
     
 end
 figure(1)
 hold on
 plot(XI(1,:), Layer_on_cliff,'r-.','linewidth',[1])
 %plot(XI(1,:), Layer(pos,YI(:,pos)),'b.')
 %plot(XI(1,:), Layer(1,:),'m.')
%% 
%now add the slope of the features in the inverse solution
%requires that you have the mesa, break, creek figure open as fig 1
% Upper plane - use two points: 
 
 x1 = 4.962e+5
 y1 = 5.41e+05
 z1 = 2155 + (2171-2155)/3;
 
 x2 = 4.957e+5
 y2 = 5.41e+05
 z2 = 2171 - (2171-2155)/3;
 x=495100:496500;
p=polyfit([x1 x2], [z1 z2], 1);
line=p(1)*x+p(2) ;
 figure(1)
 %plot(x,line,'c--') %MR - commented out 7/13/17
 
% Lower plane - use two points: 
 
 x1l = 4.959e+5
 y1l = 5.41e+05
 z1l = 2139 + (2171-2155)*2/3;
 
 x2l = 4.957e+5
 y2l = 5.41e+05
 z2l = 2139 + (2171-2155)/2;
 
 
pl=polyfit([x1l x2l], [z1l z2l], 1);
linel=pl(1)*x+pl(2) ;
 figure(1)
 plot(x,linel,'b-','linewidth',[2])
 
 %Upper pane using a higher interface
  x1h = 4.956e+5
 %this is wrong y1h = 5.41e+05
 z1h = 2187 - (2187-2171)*.1;
 
 x2h = 4.96e+5
 %wrong y2h = 5.41e+05
 z2h = 2187 - (2187-2171)/2;
 
 
ph=polyfit([x1h x2h], [z1h z2h], 1);
lineh=ph(1)*x+ph(2) ;
 figure(1)
 plot(x,lineh,'r-','linewidth',[2])
 
 
 axis([495000, 496500, 2050, 2250])


