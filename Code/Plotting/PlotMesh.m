function PlotMesh(filename,coord,tri)
%PlotMesh Plot the trimetric CATIA part mesh
%   PlotMesh(FILENAME,COORD,TRI) will plot trimetric CATIA part mesh
%   and save it in Output\Checking\filename-mesh.fig

% Load global file separator fs
global fs

% Create plot
h = figure();
trisurf(tri,coord(:,1),coord(:,2),coord(:,3),0.5*ones(size(coord,1),1));
title([filename ' trimetric mesh'])
colorbar
% Save plot in Output folder and close figure
saveas(h,['..' fs '..' fs 'Output' fs 'Checking' fs filename '-mesh.fig']);
close(h);

% Display extensive comments only
ScreenComment('',['Plot and save ' filename ' mesh']);
end