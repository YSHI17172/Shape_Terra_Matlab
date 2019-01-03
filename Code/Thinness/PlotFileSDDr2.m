function []=PlotFileSDDr2(filename)
%PlotFileSDDr2 Plots part thin features.
%   PlotFileSDDr2(FILENAME)) plots the found thin features of FILENAME and 
%   saves plot as  SDDrThinFeature.fig in ShapeTerra\Output\FILENAME

% Load variables for thin feature plotting
file=strcat('ShapeTerra\Output\',filename,'\partrecord.mat');
load(file,'tri','coord','features','SDDr','MinD','colorv');

colorv2=colorv;
% colorv2=colorv/max(colorv);
ind=find(features==0);
x1=coord([SDDr(1,:) SDDr(2,:)],1);
y1=coord([SDDr(1,:) SDDr(2,:)],2);
z1=coord([SDDr(1,:) SDDr(2,:)],3);

x2=coord([MinD(1,:) MinD(2,:)],1);
y2=coord([MinD(1,:) MinD(2,:)],2);
z2=coord([MinD(1,:) MinD(2,:)],3);
coord2=coord;
coord2(ind,:)=[];
colorv2(ind)=[];

% Plot thin feature figure and save result
h=figure;
hold on
%trimesh(tri,coord(:,1),coord(:,2),coord(:,3),[colorv2 zeros(size(coord,1),2)])
scatter3(coord2(:,1),coord2(:,2),coord2(:,3),100,colorv2,'filled');
axis equal
plot3(x1,y1,z1,'o','MarkerFaceColor',[1 0 0],'MarkerSize',10);
plot3(x2,y2,z2,'o','MarkerFaceColor',[0 0 1],'MarkerSize',10);
hold off
title(filename)
saveas(h,['ShapeTerra\Output\' filename '\SDDrThinFeature.fig']); 
close(h);
end