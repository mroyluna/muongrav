%readmany.m -- 12//16 -- by Megan Lewis; modified -- MR
clear all
close all

cd results_11_30_16/RMS

files=dir('rms_*.dat') 
for i=1:length(files)
    eval(['load ' files(i).name]);
    nAmes=num2str(files(i).name);
    naMes=strrep(nAmes, '.dat','');
    names=strrep(naMes, '-','_');
    if strfind(names, '2130')
        allrms_2130z(i,:)=eval([names]);
        allrms_2130=allrms_2130z(any(allrms_2130z,2),:);
    elseif strfind(names, '2150')
        allrms_2150z(i,:)=eval([names]);
        allrms_2150=allrms_2150z(any(allrms_2150z,2),:);
    elseif strfind(names, '2170')
        allrms_2170z(i,:)=eval([names]);
        allrms_2170=allrms_2170z(any(allrms_2170z,2),:);
    else
        allrms_2190z(i,:)=eval([names]);
        allrms_2190=allrms_2190z(any(allrms_2190z,2),:);
    end
    %rmsl(i)=eval([names '(6)']);
end

cd ../../

layerelev = [2130:20:2190];
%layerelev = 2130;
for i = 1: numel(layerelev)
    lstr      = num2str(layerelev(i));
    eval(['UpDensity    = allrms_' lstr '(:,2);'])
    eval(['DeltaDensity = allrms_' lstr '(:,3);'])
    eval(['rmsU = allrms_' lstr '(:,5);'])
    eval(['rmsL = allrms_' lstr '(:,6);'])
    
    LowDensity   = DeltaDensity + UpDensity;
    rhospaceU=[min(UpDensity):50:max(UpDensity)];
    rhospaceL=[min(LowDensity):50:max(LowDensity)];
    npts = numel(UpDensity);
    % %make array for contouring
    xr = rhospaceL;
    yr = rhospaceU';
    rmsUint=griddata(LowDensity,UpDensity,rmsU,xr,yr);
    rmsLint=griddata(LowDensity,UpDensity,rmsL,xr,yr);
    figure;
%     plot3(LowDensity,UpDensity,rmsU,'o');hold on
%     mesh(xr,yr,rmsUint);
    C1 = contour(xr,yr,rmsUint,'k-');hold on
    clabel(C1,'FontSize',13,'Color','k','FontWeight','bold')
    C2 = contour(xr,yr,rmsLint,'r-');
    clabel(C2,'FontSize',13,'Color','r','FontWeight','bold')
    title(['LayerElev=' lstr ', Black = Upper Stations, Red = Lower Stations'])
    xlabel('Lower Layer Density')
    ylabel('Upper Layer Density')
end    

