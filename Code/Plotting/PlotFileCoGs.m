function []=PlotFileCoGs(folderpath)
%PlotFileCoGs Plots persistence clusters with their CoGs
%   PlotFileCoGs(FOLDERPATH) Plots persistence clusters with their CoGs
%   for part belonging to to Output\FOLDERPATH\CoGClusteredPersistance-%.fig

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load files from output folder for folderpath
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'tri','coord','clusters','similold','COGs','COGPs');

% Define x, y and z coordinates
x=COGs(:,1,:);
y=COGs(:,2,:);
z=COGs(:,3,:);

% Plot and save CoG figures
ScreenComment('Saving CoGs plots','')
h=figure;
for i=1:length(similold)
    trimesh(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:),'edgecolor',...
        [0 0 0],'marker','o','markerfacecolor','flat','markersize',7);
    axis equal
    alpha 0;
    hold on
    plot3(x(:,i),y(:,i),z(:,i),'o','MarkerFaceColor',[0 0 0],'MarkerSize',18);
    hold off
    title(strcat(partname,' ',num2str(similold(i)*100),'%'))
    colorbar 
    ScreenComment('',['Saving CoG ' num2str(similold(i)*100) ...
        '%-peristence similarity cluster plot'])
    saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'CoG-ClusteredPersistence-'...
        num2str(similold(i)*100) '.fig']);  
end

close(h);
end