function []=PlotFileLO(folderpath)
%PlotFileLO Plots the left-overs of the part.
%   PlotFileLO(FOLDERPATH) Plots the left-overs of part with data in
%   FOLDERPATH/partrecord.mat and saves plot as LeftOvers.fig in FOLDERPATH

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load clusters and other variables for feature identification
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'tri','coord','clusters','features');

% Generate level 2 screen comment
ScreenComment('',['Generating ' partname ' left-over plots']);

% Create left-over plot
h=figure;
i=size(clusters,1)/2;
coord(features>0,:)=NaN;
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:));
axis equal
title([partname ' left-overs'])
colorbar

% Save plot and close figure
saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'LeftOvers.fig']);    
close(h);
end