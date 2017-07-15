comparelayer = [Layer(:)'; ElevI(:)'];
comparelayer2= [Layer2(:)'; ElevI(:)'];
minElev=ones(1,n*n);
minElev2=ones(1,n*n);
minElev = min(comparelayer);
minElev2 = min(comparelayer2);

top=area(YI(80001:80400)', ElevI(80001:80400)');
top.FaceColor=[.7 .7 0.9];
hold on
mid=area(YI(80001:80400)', minElev(80001:80400)');
mid.FaceColor=[1 .7 .7];
bot=area(YI(80001:80400)', minElev2(80001:80400)');
bot.FaceColor=[.9 0.7 0];
plot(northing,elevations,'ko','linewidth', 2)
set(gca,'xlim',[5.4060e5,5.412e5],'ylim',[2050,2260])
set(gca, 'fontsize',[16])
xlabel('Northing, m'); ylabel('Elevation, m')