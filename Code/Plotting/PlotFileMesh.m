function PlotFileMesh(folderpath)
%PlotFileMesh Plot the trimetric CATIA part mesh
%   PlotFileMesh(FOLDERPATH) will plot trimetric CATIA part mesh
%   and save it in folderpath\mesh.fig

% Load global file separator fs
global fs

% Determine partname folder
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load data from the file part
file=(['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat']);
load(file,'tri','coord');

% Create plot
h = figure();
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),0.5*ones(size(coord,1),1));
title([partname ' trimetric mesh'])
colorbar
% Save plot in Output folder and close figure
saveas(h,(['ShapeTerra' fs 'Output' fs folderpath fs 'mesh.fig']));
close(h);

% Display extensive comments only
ScreenComment('',['Plot and save ' partname ' mesh']);
end