function []=PlotFileSDDr(folderpath)
%PlotFileSDDr Plots part thin features.
%   PlotFileSDDr(FOLDERPATH)) plots the found thin features of part 
%   belonging to FOLDERPATH/partrecord.mat and saves plot as 
%   SDDrThinFeature.fig in Output\FOLDERPATH

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load variables for thin feature plotting
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'tri','coord','features','SDDr','MinD','colorv');

% Generate level 2 screen comment
ScreenComment('',['Generating ' partname ' thin features plot']);

colorv2=colorv;
% colorv2=colorv/max(colorv);
ind=find(features==0);
x1=coord([SDDr(1,:) SDDr(2,:)],1);
y1=coord([SDDr(1,:) SDDr(2,:)],2);
z1=coord([SDDr(1,:) SDDr(2,:)],3);

x2=coord([MinD(1,:) MinD(2,:)],1);
y2=coord([MinD(1,:) MinD(2,:)],2);
z2=coord([MinD(1,:) MinD(2,:)],3);

coord(ind,:)=NaN;
    
% Plot thin feature figure and save result
h=figure;
%trimesh(tri,coord(:,1),coord(:,2),coord(:,3),[colorv2 zeros(size(coord,1),2)])
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),colorv2)
axis equal
hold on
plot3(x1,y1,z1,'o','MarkerFaceColor',[1 0 0],'MarkerSize',10);
plot3(x2,y2,z2,'o','MarkerFaceColor',[0 0 1],'MarkerSize',10);
hold off
title([partname ' SDDR thin feature'])
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'SDDrThinFeature.fig']);    
close(h);
end