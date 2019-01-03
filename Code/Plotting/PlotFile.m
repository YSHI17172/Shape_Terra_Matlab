function []=PlotFile(folderpath)
%PlotFile Plots and save part clusters with same similarity persistence levels 
%   PlotFile(FOLDERPATH) plots clustered persistence similarity levels for 
%   specific part and saves the generated similarity figures as
%   folderpath\ClusteredPersistence-%.fig

% Load global file separator fs
global fs

% Retreive filename from folderpath
indsep = strfind(folderpath,fs);
%tradename = folderpath(1:indsep(1)-1);
partname = folderpath(indsep(1)+1:indsep(2)-1);
%meshname = folderpath(indsep(2)+1:indsep(3)-1);

% Load files from output folder for folderpath
file=['ShapeTerra' fs 'Output' fs folderpath fs 'partrecord.mat'];
load(file,'tri','coord','clusters','similold');

% Plot and save figures
ScreenComment('Saving Cluster plots','')
h=figure();
for i=1:length(similold)
    trisurf(tri,coord(:,1),coord(:,2),coord(:,3),clusters(2*i,:));
    axis equal
    title(strcat(partname,' ',num2str(similold(i)*100),'%'))
    colorbar
    ScreenComment('',['Saving ' num2str(similold(i)*100) '%' ...
        ' peristence similarity cluster plot'])
    saveas(h,['ShapeTerra' fs 'Output' fs folderpath fs 'ClusteredPersistence-'...
        num2str(similold(i)*100) '.fig']);    
end

close(h);
end